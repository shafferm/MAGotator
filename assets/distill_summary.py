import pandas as pd
import argparse

def distill_summary(combined_annotations, genome_summary_form, target_id_counts, output, add_modules):
    # Read input files into pandas dataframes
    combined_data = pd.read_csv(combined_annotations, sep='\t')
    genome_summary_data = pd.read_csv(genome_summary_form, sep='\t')
    target_id_counts_data = pd.read_csv(target_id_counts, sep='\t', header=None, names=['gene_id'])

    # Create a set of all gene_id values from target_id_counts
    valid_gene_ids = set(target_id_counts_data['gene_id'])

    # Process additional module files
    for i, add_module_file in enumerate(add_modules, start=1):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')

            # Print debugging information
            print(f"Additional Module {i} Columns: {additional_module_data.columns}")

            # Add gene_id values from additional modules to the set
            valid_gene_ids.update(additional_module_data['gene_id'])

            # Dynamically identify columns ending with "_id" in combined_data
            id_columns = [col for col in combined_data.columns if col.endswith('_id')]

            # Update merging logic to consider dynamically identified gene_id columns
            for gene_id_col in id_columns:
                if gene_id_col in additional_module_data.columns:
                    combined_data = pd.merge(combined_data, additional_module_data, left_on=gene_id_col, right_on=gene_id_col, how='left')

    # Filter rows with valid gene_id values
    combined_data = combined_data[combined_data['gene_id'].isin(valid_gene_ids)]

    # Merge with genome_summary_data on gene_id
    result = pd.merge(combined_data, genome_summary_data, on='gene_id', how='left')

    # Write the result to the output file
    result.to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts file')
    parser.add_argument('--output', required=True, help='Path to the output file')
    parser.add_argument('--add_module1', required=False, help='Path to the additional module1 file')
    parser.add_argument('--add_module2', required=False, help='Path to the additional module2 file')
    parser.add_argument('--add_module3', required=False, help='Path to the additional module3 file')
    parser.add_argument('--add_module4', required=False, help='Path to the additional module4 file')
    parser.add_argument('--add_module5', required=False, help='Path to the additional module5 file')

    args = parser.parse_args()
    add_modules = [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]

    distill_summary(args.combined_annotations, args.genome_summary_form, args.target_id_counts, args.output, add_modules)
