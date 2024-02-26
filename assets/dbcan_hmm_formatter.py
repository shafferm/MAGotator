import pandas as pd
import argparse
import sqlite3

def get_sig_row(row):
    return row['full_evalue'] < 1e-18

def calculate_bit_score(row):
    return row['full_score'] / row['domain_number']

def fetch_descriptions_from_db(target_ids, db_file):
    conn = sqlite3.connect(db_file)
    descriptions = {}
    for target_id in target_ids:
        cursor = conn.execute("SELECT description, ec FROM dbcan_description WHERE id=?", (target_id,))
        row = cursor.fetchone()
        if row:
            description, ec = row
            descriptions[target_id] = {'description': description, 'ec': ec}
        else:
            descriptions[target_id] = {'description': "", 'ec': ""}  # Handle case where description is not found
    conn.close()
    return descriptions


def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    parser.add_argument("--db_file", type=str, help="Path to the SQLite database file.")

    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)
    print(f"Loaded HMM search results from: {args.hits_csv}")

    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)

    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)

    # Fetch descriptions from the database
    target_ids = hits_df['target_id'].unique()
    descriptions = fetch_descriptions_from_db(target_ids, args.db_file)

    # Assign descriptions and ECs to hits
    hits_df['dbcan_description'] = hits_df['target_id'].map(lambda x: descriptions[x]['description'])
    hits_df['dbcan_EC'] = hits_df['target_id'].map(lambda x: descriptions[x]['ec'])

    print("Saving the formatted output to CSV...")
    selected_columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'target_id', 'score_rank', 'bitScore', 'dbcan_description', 'dbcan_EC']
    modified_columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'dbcan_id', 'dbcan_score_rank', 'dbcan_bitScore', 'dbcan_description', 'dbcan_EC']

    # Ensure the columns exist in the DataFrame before renaming
    if set(selected_columns).issubset(hits_df.columns):
        # Rename the selected columns
        hits_df.rename(columns=dict(zip(selected_columns, modified_columns)), inplace=True)

        # Save the formatted output to CSV
        try:
            hits_df[modified_columns].to_csv(args.output, index=False)
            print(f"Formatted output saved to: {args.output}")
        except Exception as e:
            print(f"Error occurred while saving the formatted output: {e}")

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
