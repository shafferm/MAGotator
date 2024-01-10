import argparse
import pandas as pd

def distill_summary(combined_annotations_file, genome_summary_form_file, output_file):
    try:
        # Read input files into dataframes
        combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
        genome_summary_form = pd.read_csv(genome_summary_form_file, sep='\t')

        # Check if 'gene_id' column exists in genome_summary_form
        if 'gene_id' not in genome_summary_form.columns:
            raise KeyError("'gene_id' column not found in genome_summary_form file.")

        # Merge data based on gene_id
        merged_data = pd.merge(genome_summary_form, combined_annotations, how='left',left_on='gene_id', right_on='gene_id')

        # Write output to TSV file
        merged_data.to_csv(output_file, sep='\t', index=False)

        print("Distill summary completed successfully.")

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    # Call distill_summary function
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
