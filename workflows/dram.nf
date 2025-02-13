/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dram_pipeline'
include { getFastaChannel        } from '../subworkflows/local/utils_pipeline_setup.nf'

// Pipeline steps
include { RENAME_FASTA           } from "${projectDir}/modules/local/rename/rename_fasta.nf"
include { CALL                   } from "${projectDir}/subworkflows/local/call.nf"
include { COLLECT_RNA            } from "${projectDir}/subworkflows/local/collect_rna.nf"
include { MERGE                  } from "${projectDir}/subworkflows/local/merge.nf"
include { ANNOTATE               } from "${projectDir}/subworkflows/local/annotate.nf"
include { ADD_AND_COMBINE        } from "${projectDir}/subworkflows/local/add_and_combine.nf"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DRAM {

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'dram_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }



    default_channel = Channel.fromPath(params.distill_dummy_sheet)
    ch_fasta = getFastaChannel(params.input_fasta, params.fasta_fmt)

    if( params.rename ) {
        RENAME_FASTA( ch_fasta )
        ch_fasta = RENAME_FASTA.out.renamed_fasta
    }
    if (params.merge_annotations){
        MERGE()
    }else {
        ch_quast_stats = default_channel
        ch_gene_locs = default_channel
        ch_called_proteins = default_channel
        ch_collected_fna = default_channel

        if (params.call){
            CALL( ch_fasta )
            ch_quast_stats = CALL.out.ch_quast_stats
            ch_gene_locs = CALL.out.ch_gene_locs
            ch_called_proteins = CALL.out.ch_called_proteins
            ch_collected_fna = CALL.out.ch_collected_fna

        } else{
            
        }

        // TODO: Simplify this when we implement distill
        if (params.call || params.distill_topic || params.distill_ecosystem || params.distill_custom){
            COLLECT_RNA( ch_fasta, default_channel )
        }
        if (params.annotate){
            ANNOTATE( ch_fasta, ch_gene_locs, ch_called_proteins, default_channel )
            
        }
        if (params.annotate || params.distill_topic || params.distill_ecosystem || params.distill_custom){
            if (params.annotate){
                ch_combined_annotations = ANNOTATE.out.ch_combined_annotations
            } else {
                ch_combined_annotations = Channel
                    .fromPath(params.annotations, checkIfExists: true)
                    .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --annotate, you must provide an annotations TSV file (--annotations <path>) with approprite formatting. Cannot find any called gene files matching: ${params.annotations}\nNB: Path needs to follow pattern: path/to/directory/" }
            }
            ADD_AND_COMBINE( ch_combined_annotations )
        }
    }

    // //
    // // MODULE: MultiQC
    // //
    // ch_multiqc_config        = Channel.fromPath(
    //     "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ?
    //     Channel.fromPath(params.multiqc_config, checkIfExists: true) :
    //     Channel.empty()
    // ch_multiqc_logo          = params.multiqc_logo ?
    //     Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
    //     Channel.empty()

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     []
    // )

    emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
