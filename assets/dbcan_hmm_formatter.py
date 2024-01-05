import pandas as pd
import argparse

def get_sig_row(row):
    return row['full_evalue'] < 1e-5

def calculate_bit_score(row):
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def extract_subfamily(row, ch_dbcan_subfam, ch_dbcan_fam):
    target_id = row['target_id'].replace('.hmm', '')

    matching_rows_subfam = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]
    matching_rows_fam = ch_dbcan_fam[ch_dbcan_fam['target_id'] == target_id]

    if not matching_rows_subfam.empty:
        return matching_rows_subfam.iloc[0]['subfamily']
    elif not matching_rows_fam.empty:
        return matching_rows_fam.iloc[0]['subfamily']
    else:
        return ""


def extract_subfam_ec(row, ch_dbcan_subfam):
    target_id = row['target_id'].replace('.hmm', '')
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]

    if not matching_rows.empty:
        # Filter out rows with NaN values in 'subfam-EC'
        matching_rows = matching_rows.dropna(subset=['subfam-EC'])

        if not matching_rows.empty:
            # Filter values that match the EC pattern
            ec_values = matching_rows['subfam-EC'].str.extract(r'(\d+\.\d+\.\d+\.\d+)', expand=False).dropna()
            if not ec_values.empty:
                # Concatenate all matching EC values with "; "
                return "; ".join(ec_values.unique())

    return ""

def extract_subfam_genbank(row, ch_dbcan_subfam):
    target_id = row['target_id'].replace('.hmm', '')
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]

    if not matching_rows.empty:
        # Filter out rows with NaN values in 'subfam-GenBank'
        matching_rows = matching_rows.dropna(subset=['subfam-GenBank'])

        if not matching_rows.empty:
            # Concatenate all GenBank values with "; "
            return "; ".join(matching_rows['subfam-GenBank'].astype(str).unique())

    return ""



def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the fam file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfam file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)
    print("Loading subfam file...")
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", comment='#', header=None, names=['target_id', 'subfamily', 'subfam-GenBank', 'subfam-EC', 'score'], engine='python')
    print("Loading fam file...")
    ch_dbcan_fam = pd.read_csv(args.fam, comment='#', header=None, names=['target_id', 'subfamily'], engine='python', on_bad_lines='skip', delimiter='\t', usecols=[0, 1], quoting=3)

    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)

    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    hits_df['subfamily'] = hits_df.apply(lambda row: extract_subfamily(row, ch_dbcan_subfam, ch_dbcan_fam), axis=1)

    hits_df['subfam-GenBank'] = hits_df.apply(lambda row: extract_subfam_genbank(row, ch_dbcan_subfam), axis=1)
    hits_df['subfam-EC'] = hits_df.apply(lambda row: extract_subfam_ec(row, ch_dbcan_subfam), axis=1)

    sig_hits_df = hits_df[hits_df.apply(get_sig_row, axis=1)]
    sig_hits_df = sig_hits_df.sort_values(by='score_rank')

    print("Saving the formatted output to CSV...")
    selected_columns = ['query_id', 'target_id', 'score_rank', 'bitScore', 'subfamily', 'subfam-GenBank', 'subfam-EC']
    sig_hits_df[selected_columns].to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
