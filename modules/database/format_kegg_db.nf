process FORMAT_KEGG_DB {

    input:
    path( ch_kegg_pep, stageAs: "kegg_pep_loc" )
    path( gene_ko_link_loc, stageAs: "genes_ko_link.gz")
    path( ch_format_kegg_db_script )

    script:
    """

    python ${ch_format_kegg_db_script} --kegg_loc ${ch_kegg_pep} --gene_ko_link_loc ${gene_ko_link_loc} --download_date ${params.kegg_download_date} --threads ${params.threads}

    """
}
