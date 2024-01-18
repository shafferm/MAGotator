process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_trnas

    output:
    path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + samples)

    # Iterate through each TSV file
    for file in tsv_files:
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file, sep='\t', skiprows=[1])

        # Construct the gene_id column
        df['gene_id'] = df['type'] + ' (' + df['codon'] + ')'

        # Construct the gene_description column
        df['gene_description'] = df['type'] + ' tRNA with ' + df['codon'] + ' Codon'

        # Construct the module column
        df['module'] = df['type'] + ' tRNA'

        # Add constant values to header and subheader columns
        df['header'] = 'tRNA'
        df['subheader'] = ''

        # Extract sample name from the file name
        sample_name = os.path.basename(file).replace("_processed_trnas.tsv", "")

        # Deduplicate the gene_id column
        df_deduplicated = df.drop_duplicates(subset=['gene_id'])

        # Update the corresponding columns in the collected_data DataFrame
        collected_data.loc[:, 'gene_id'] = df_deduplicated['gene_id']
        collected_data.loc[:, 'gene_description'] = df_deduplicated['gene_description']
        collected_data.loc[:, 'module'] = df_deduplicated['module']
        collected_data.loc[:, 'header'] = df_deduplicated['header']
        collected_data.loc[:, 'subheader'] = df_deduplicated['subheader']
        
        # Count occurrences of each unique gene_id value
        gene_id_counts = df['gene_id'].value_counts()
        
        # Populate counts in the corresponding sample-named columns
        for unique_gene_id, count in gene_id_counts.items():
            collected_data.loc[collected_data['gene_id'] == unique_gene_id, sample_name] = count

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
