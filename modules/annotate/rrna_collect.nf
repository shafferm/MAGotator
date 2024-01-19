process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    val combined_rrnas

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
    sample_counts = {sample: [] for sample in samples}

    # Iterate through each input file
    for file in tsv_files:
        sample_name = os.path.basename(file).replace("_processed_rrnas.tsv", "")
        input_df = pd.read_csv(file, sep='\t')
        
        # Populate gene_id column with collective values from the "type" column
        gene_ids = input_df['type'].tolist()
        collected_data['gene_id'] = gene_ids
        
        # De-replicate gene_id column to ensure no duplicates
        collected_data['gene_id'] = collected_data['gene_id'].astype(str).apply(lambda x: list(set(x.split(','))))

        # Populate gene_description column
        collected_data['gene_description'] = [f"{gene} gene" for gene in collected_data['gene_id']]
        
        # Set module column values to "rRNA"
        collected_data['module'] = 'rRNA'

        # Count occurrences of each type value for each sample
        for gene_id in gene_ids:
            count_dict = Counter(gene_id.split(','))
            for gene, count in count_dict.items():
                sample_counts[sample_name].append(count)

    # Add sample count values to the output DataFrame
    for sample, counts in sample_counts.items():
        collected_data[sample] = counts

    # Save collected data to collected_rrnas.tsv
    collected_data.to_csv("collected_rrnas.tsv", sep='\t', index=False)
    """
}
