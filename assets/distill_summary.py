import pandas as pd
import argparse
import glob
import logging

logging.basicConfig(level=logging.DEBUG)

def is_null_content(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
    return content == "NULL"

def distill_summary(combined_annotations_path, target_id_counts_df, output_path):
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]
    distill_sheets = glob.glob('*_distill_sheet.tsv')
    distill_summary_df = pd.DataFrame()
    additional_columns = set()
    has_target_id_column = False

    for distill_sheet in distill_sheets:
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue
        logging.info(f"Processing distill sheet: {distill_sheet}")
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        for index, row in distill_df.iterrows():
            gene_id = row['gene_id']
            gene_description = row['gene_description']
            pathway = row.get('pathway', None)
            topic_ecosystem = row.get('topic_ecosystem', None)
            category = row.get('category', None)
            subcategory = row.get('subcategory', None)

            # Check potential_gene_id_columns first
            for col in combined_annotations_df.columns:
                if col.endswith('_id') and col != "query_id":  # Exclude query_id
                    matched_indices = combined_annotations_df[col].str.contains(gene_id, na=False)
                    if matched_indices.any():
                        combined_ids = gene_id  # Use the matched gene_id directly
                        row_data = {'gene_id': combined_ids,
                                    'gene_description': gene_description,
                                    'pathway': pathway,
                                    'topic_ecosystem': topic_ecosystem,
                                    'category': category,
                                    'subcategory': subcategory}
                        # Add additional columns from distill sheet
                        for additional_col in set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'}:
                            if additional_col == 'target_id':
                                has_target_id_column = True
                            row_data[additional_col] = row.get(additional_col, None)
                        distill_summary_df = distill_summary_df.append(row_data, ignore_index=True)
                        break
            else:
                # Check potential_ec_columns
                for col in combined_annotations_df.columns:
                    if col.endswith('_EC'):
                        matched_indices = combined_annotations_df[col].str.contains(gene_id, na=False)
                        if matched_indices.any():
                            # Create a new row for each unique gene_id found
                            for combined_id in combined_annotations_df.loc[matched_indices, col.replace('_EC', '_id')]:
                                if gene_id in combined_annotations_df[col.replace('_EC', '_id')].values:
                                    gene_description += f'; {gene_id}'
                                    row_data = {'gene_id': combined_id,
                                                'gene_description': gene_description,  # Append EC to gene_description
                                                'pathway': pathway,
                                                'topic_ecosystem': topic_ecosystem,
                                                'category': category,
                                                'subcategory': subcategory}
                                    # Add additional columns from distill sheet
                                    for additional_col in set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'}:
                                        if additional_col == 'target_id':
                                            has_target_id_column = True
                                        row_data[additional_col] = row.get(additional_col, None)
                                    distill_summary_df = distill_summary_df.append(row_data, ignore_index=True)
                            break
                else:
                    logging.info(f"No match found for gene ID {gene_id}.")

    # Merge distill_summary_df with target_id_counts_df
    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'],
                                  how='left')
    
    # Remove the 'target_id' column if it wasn't added as an additional column
    if not has_target_id_column:
        distill_summary_df.drop('target_id', axis=1, inplace=True, errors='ignore')

    # Define columns to output
    required_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    additional_columns = [col for col in distill_summary_df.columns if col not in required_columns and col not in target_id_counts_df.columns]

    columns_to_output = required_columns + additional_columns + list(target_id_counts_df.columns)

    # Ensure all required columns are present
    for col in columns_to_output:
        if col not in distill_summary_df.columns:
            distill_summary_df[col] = None

    # Drop duplicates based on subset of columns
    deduplicated_df = distill_summary_df.drop_duplicates(subset=required_columns, ignore_index=True).copy()

    # Write output to file
    deduplicated_df.to_csv(output_path, sep='\t', index=False, columns=columns_to_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True,
                        help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True,
                        help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True,
                        help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
