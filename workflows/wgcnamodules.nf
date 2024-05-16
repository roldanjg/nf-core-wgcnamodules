/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
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

    if (params.genes && params.diff_exp_genes == true) {
        exit 1, 'Introuducing a list of genes in --genes is not compatible with setting --diff_exp_genes true.' 
    }

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW:  DIFFERENTIAL EXPRESSION ANALISYS
    //

    if (params.diff_exp_genes == true) {
        DIFFERENTIALABUNDANCE()
        ch_versions = ch_versions.mix(DIFFERENTIALABUNDANCE.out.versions)
    }

    input_wgcna = file(params.input, checkIfExists: true)
    if (params.contrast) { contrast_wgcna = file(params.contrast, checkIfExists: true) } else { contrast_wgcna = null }
    if (params.genes) {
        genes = file(params.genes, checkIfExists: true)
        } else if ( params.diff_exp_genes == true) {
             genes = DIFFERENTIALABUNDANCE.out.tables.collect{it[1]}.map{it.join("+c+")}
             } else (
                genes = "all"
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

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_wgcnamodules_software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
