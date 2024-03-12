process ADD_BIN_QUALITY {

    errorStrategy 'finish'

    input:
    path( combined_annotations, stageAs: "input-raw-annotations.tsv" )
    file( ch_bin_quality )

    output:
    path("raw-annotations.tsv"), emit: annots_bin_quality_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Load combined_annotations.tsv
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')

    # Load checkm TSV
    checkm_path = "${ch_bin_quality}"
    checkm_data = pd.read_csv(checkm_path, sep='\t')

    # Replace "." with "-" in the sample column for comparison
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    checkm_data.iloc[:, 0] = checkm_data.iloc[:, 0].str.replace(".", "-")  # Assume first column matches "sample"

    # Merge checkm data with combined_annotations
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on=checkm_data.columns[0], how="left")

    # After merging, Completeness and Contamination are updated only if they are NaN in the combined_annotations
    merged_data['Completeness'] = merged_data['Completeness_x'].fillna(merged_data['Completeness_y'])
    merged_data['Contamination'] = merged_data['Contamination_x'].fillna(merged_data['Contamination_y'])

    # Drop temporary columns generated by the merge
    merged_data.drop(columns=[col for col in merged_data.columns if '_x' in col or '_y' in col], inplace=True)

    # Ensure "rank" column is correctly preserved
    # This step might be redundant if the "rank" column is not modified above, but it ensures clarity.
    if 'rank' not in merged_data.columns and 'rank' in combined_annotations.columns:
        merged_data['rank'] = combined_annotations['rank']

    # Drop the column from checkm_data that was used for merging if it's not 'sample'
    if checkm_data.columns[0] != 'sample':
        merged_data.drop(columns=[checkm_data.columns[0]], inplace=True)

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    """
}