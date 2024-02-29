import pandas as pd
import sqlite3
import argparse
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate a multi-sheet XLSX document from distill sheets and a SQLite database.')
    parser.add_argument('--target_id_counts', type=str, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--db_name', type=str, help='Name of the SQLite database file.')
    parser.add_argument('--distill_sheets', nargs='+', help='List of paths to distill sheets.')
    parser.add_argument('--rrna_file', type=str, help='Path to the rrna_sheet.tsv file.', default=None)
    parser.add_argument('--trna_file', type=str, help='Path to the trna_sheet.tsv file.', default=None)
    parser.add_argument('--output_file', type=str, help='Path to the output XLSX file.')
    return parser.parse_args()

def compile_target_id_counts(target_id_counts):
    return pd.read_csv(target_id_counts, sep='\t')

def read_distill_sheets(distill_sheets):
    sheets_data = {}
    for sheet_path in distill_sheets:
        if file_contains_data(sheet_path):
            df = pd.read_csv(sheet_path, sep='\t')
            topic = df['topic_ecosystem'].unique().tolist()
            sheets_data.update({sheet_path: {'dataframe': df, 'topics': topic}})
        else:
            print(f"Skipping {sheet_path} as it contains 'NULL'.")
    return sheets_data

def file_contains_data(file_path):
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            return first_line != "NULL"
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return False

def compile_genome_stats(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    # Check for the presence of optional columns in the database
    cursor.execute("PRAGMA table_info(annotations)")
    columns_info = cursor.fetchall()
    column_names = [info[1] for info in columns_info]
    
    # Start with a query to select distinct samples
    query_base = "SELECT DISTINCT sample FROM annotations"
    df_samples = pd.read_sql_query(query_base, conn)
    
    # Initialize the genome_stats DataFrame with sample names
    df_genome_stats = pd.DataFrame({"sample": df_samples['sample']})
    
    # Dynamically build a query to include optional columns if they are present
    optional_columns = ['taxonomy', 'Completeness', 'Contamination']
    select_clauses = []
    for col in optional_columns:
        if col in column_names:
            # Prepare a SELECT clause to calculate average for numeric columns and group_concat for taxonomy
            if col in ['Completeness', 'Contamination']:
                select_clauses.append(f"AVG({col}) AS {col}")
            else:  # Assuming 'taxonomy' is not a numeric column
                select_clauses.append(f"GROUP_CONCAT(DISTINCT {col}) AS {col}")
    
    # If there are optional columns to include, modify the base query
    if select_clauses:
        query = f"SELECT sample, {', '.join(select_clauses)} FROM annotations GROUP BY sample"
        df_stats = pd.read_sql_query(query, conn)
        df_genome_stats = pd.merge(df_genome_stats, df_stats, on="sample", how="left")
    else:
        # If no optional columns, the genome_stats will only contain sample names at this point
        pass
    
    conn.close()
    return df_genome_stats

def add_sheet_from_dataframe(wb, df, sheet_name):
    ws = wb.create_sheet(title=sheet_name)
    for r in dataframe_to_rows(df, index=False, header=True):
        ws.append(r)

def query_annotations_for_gene_ids(db_name, gene_ids):
    """Fetch annotations for given gene IDs if they exist in the annotations database."""
    conn = sqlite3.connect(db_name)
    placeholders = ', '.join('?' for _ in gene_ids)
    query = f"SELECT DISTINCT gene_id FROM annotations WHERE gene_id IN ({placeholders})"
    df = pd.read_sql_query(query, conn, params=gene_ids)
    conn.close()
    return df

def add_rrna_trna_sheets(wb, rrna_file, trna_file):
    for tsv_file, sheet_name in [(rrna_file, 'rRNA'), (trna_file, 'tRNA')]:
        if tsv_file and file_contains_data(tsv_file):
            df = pd.read_csv(tsv_file, sep='\t')
            add_sheet_from_dataframe(wb, df, sheet_name)

def main():
    args = parse_arguments()

    wb = Workbook()
    wb.remove(wb.active)  # Remove the default sheet

    genome_stats_df = compile_genome_stats(args.db_name)
    add_sheet_from_dataframe(wb, genome_stats_df, "genome_stats")

    target_id_counts_df = compile_target_id_counts(args.target_id_counts)
    distill_data = read_distill_sheets(args.distill_sheets)
    for sheet_path, info in distill_data.items():
        df_distill = info['dataframe']
        for topic in info['topics']:
            df_topic = df_distill[df_distill['topic_ecosystem'] == topic]
            gene_ids = df_topic['gene_id'].unique().tolist()
            
            # Query the database to get a list of gene_ids that exist in annotations.db
            df_valid_gene_ids = query_annotations_for_gene_ids(args.db_name, gene_ids)
            
            # Filter the df_topic to keep only rows where gene_id exists in annotations.db
            df_topic_filtered = df_topic[df_topic['gene_id'].isin(df_valid_gene_ids['gene_id'])]
            
            # Merge the filtered df_topic with target_id_counts data
            df_merged = pd.merge(df_topic_filtered, target_id_counts_df, on="gene_id", how="left")
            df_final = df_merged.drop(columns=['query_id', 'sample', 'taxonomy', 'Completeness', 'Contamination'], errors='ignore')
            
            sheet_name = topic[:31]  # Excel sheet name character limit
            add_sheet_from_dataframe(wb, df_final, sheet_name)

    # Add rrna and trna sheets if not NULL
    add_rrna_trna_sheets(wb, args.rrna_file, args.trna_file)

    wb.save(args.output_file)

if __name__ == '__main__':
    main()
