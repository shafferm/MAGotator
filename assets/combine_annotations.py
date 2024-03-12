import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def extract_samples_and_paths(annotation_files):
    samples_and_paths = []
    for i in range(0, len(annotation_files), 2):
        sample = annotation_files[i].strip('[], ')
        path = annotation_files[i + 1].strip('[], ')
        samples_and_paths.append((sample, path))
    return samples_and_paths

def assign_rank(row):
    # Adjust this function based on your database bit score columns
    rank = 'E'  # Default rank
    if 'kegg_bitScore' in row and row['kegg_bitScore'] > 350:
        rank = 'A'
    elif 'uniref_bitScore' in row and row['uniref_bitScore'] > 350:
        rank = 'B'
    elif ('kegg_bitScore' in row and row['kegg_bitScore'] > 60) or ('uniref_bitScore' in row and row['uniref_bitScore'] > 60):
        rank = 'C'
    elif any(db in row for db in ['pfam_id', 'dbcan_id', 'merops_id']) and not (row.get('kegg_bitScore', 0) > 60 or row.get('uniref_bitScore', 0) > 60):
        rank = 'D'
    # Note: Rank 'E' is already the default
    return rank

def organize_columns(df):
    # Ensure the 'rank' column is included after 'strandedness'
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    other_columns = [col for col in df.columns if col not in base_columns]
    final_columns_order = base_columns + other_columns
    return df[final_columns_order]

def combine_annotations(annotation_files, output_file):
    samples_and_paths = extract_samples_and_paths(annotation_files)
    combined_data = pd.DataFrame()

    for sample, path in samples_and_paths:
        try:
            annotation_df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            logging.warning(f"Empty DataFrame for sample: {sample}, skipping.")
            continue
        except Exception as e:
            logging.error(f"Error loading DataFrame for sample {sample}: {str(e)}")
            continue

        annotation_df.insert(0, 'sample', sample)
        combined_data = pd.concat([combined_data, annotation_df], ignore_index=True, sort=False)

    # Drop duplicates
    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'])

    # Assign ranks based on defined criteria
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)

    # Organize and sort columns
    combined_data = organize_columns(combined_data)
    combined_data = combined_data.sort_values(by='query_id', ascending=True)

    combined_data.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine annotation files with ranks and avoid duplicating specific columns.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names.")
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output)
        logging.info(f"Combined annotations saved to {args.output}, with ranks assigned and sorted by 'query_id'.")
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
