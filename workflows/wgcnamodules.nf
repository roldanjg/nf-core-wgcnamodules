/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_wgcnamodules_pipeline'
include { DIFFERENTIALABUNDANCE } from '../subworkflows/local/differentialabundancemod'
include { CUSTOM_FILES_TRANSFORMATION } from '../modules/local/customtranformation/main'
include { WGCNA } from '../modules/local/wgcna/main'
include { TDTHUBMODULES } from '../modules/local/tdthub/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow WGCNAMODULES {

    main:

    if (params.genes && params.diffgenes == true) {
        exit 1, 'introuducing a list of genes in --genes is not compatible with setting --diffgenes true.' 
    }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW:  DIFFERENTIAL EXPRESSION ANALISYS
    //
    // if (params.diffgenes == true) {
    //     DIFFERENTIALABUNDANCE()
    //     ch_versions = ch_versions.mix(DIFFERENTIALABUNDANCE.out.versions)
    // }

    input_wgcna = file(params.input, checkIfExists: true)
    if (params.contrast) { contrast_wgcna = file(params.contrast, checkIfExists: true) } else { contrast_wgcna = null }
    if (params.genes) { 
        genes = file(params.genes, checkIfExists: true) 
        } else if ( params.diffgenes == true) {
             genes = "${params.outdir}/tables/differential" 
             } else (
                genes = "/home/joaquingr/nf-core-wgcnamodules/.nextflow.log"
             )
    tpms_wgcna = file(params.salmon_dir+"/salmon.merged.gene_tpm.tsv", checkIfExists: true)

    CUSTOM_FILES_TRANSFORMATION(
            input_wgcna,
            contrast_wgcna,
            tpms_wgcna,
            genes
        )
    ch_versions = ch_versions.mix(CUSTOM_FILES_TRANSFORMATION.out.versions)
    WGCNA(
        CUSTOM_FILES_TRANSFORMATION.out.csv
    )
    ch_versions = ch_versions.mix(WGCNA.out.versions)
    TDTHUBMODULES(
        WGCNA.out.modulesextraction
    )
    ch_versions = ch_versions.mix(TDTHUBMODULES.out.versions)
    // ch_versions = ch_versions.mix(DIFFERENTIALABUNDANCE.out.versions)

// --------------------------------------------------------QUITAR TODO ESTO
    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     ch_samplesheet
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
// --------------------------------------------------------QUITAR TODO ESTO

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

// --------------------------------------------------------QUITAR TODO ESTO
    //
    // MODULE: MultiQC
    //
    // ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    // ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    // summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )

    // emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    // versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
// --------------------------------------------------------QUITAR TODO ESTO

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
