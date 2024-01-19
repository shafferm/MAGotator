process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_rrnas

    output:
    path("collected_rrnas.tsv"), emit: rrna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    from collections import Counter

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in tsv_files]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + samples)

    # Create a dictionary to store counts for each sample
    sample_counts = {sample: Counter() for sample in samples}

    # Iterate through each input file
    for file in tsv_files:
        sample_name = os.path.basename(file).replace("_processed_rrnas.tsv", "")
        input_df = pd.read_csv(file, sep='\t')

        # Populate gene_id column with collective values from the "type" column
        gene_ids = input_df['type'].tolist()
        gene_counts = Counter(gene_ids)

        # Populate gene_description column
        collected_data = pd.concat([collected_data, pd.DataFrame({'gene_id': gene_counts.keys()})], ignore_index=True)

        # Set module column values to "rRNA"
        collected_data['module'] = 'rRNA'

        # Update counts for each sample
        sample_counts[sample_name] += gene_counts

    # Add counts to the output DataFrame
    for sample in samples:
        collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample][x] if x in sample_counts[sample] else 0)

    # Populate other columns based on the given rules
    collected_data['gene_description'] = collected_data['gene_id'] + " gene"
    collected_data['header'] = ""
    collected_data['subheader'] = ""

    # Remove the first row that includes extra column names
    collected_data = collected_data[collected_data['gene_id'] != 'gene_id']

    # Sort the whole table by the gene_id column
    collected_data.sort_values(by='gene_id', inplace=True)

    # Write the collected data to the output file
    collected_data.to_csv("collected_rrnas.tsv", sep="\t", index=False)

    """
}
