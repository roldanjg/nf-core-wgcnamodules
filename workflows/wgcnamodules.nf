/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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
    if (params.diffgenes == true) {
        DIFFERENTIALABUNDANCE()
        ch_versions = ch_versions.mix(DIFFERENTIALABUNDANCE.out.versions)
    }

    input_wgcna = file(params.input, checkIfExists: true)
    if (params.contrast) { contrast_wgcna = file(params.contrast, checkIfExists: true) } else { contrast_wgcna = null }
    if (params.genes) {
        genes = file(params.genes, checkIfExists: true)
        } else if ( params.diffgenes == true) {
             genes = DIFFERENTIALABUNDANCE.out.tables.collect{it[1]}.map{it.join("+c+")}
             } else (
                genes = "all"
             )
    // println "$projectDir"
    // DIFFERENTIALABUNDANCE.out.tables.collect{it[1]}.join(";")
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
    ch_versions = ch_versions.mix(DIFFERENTIALABUNDANCE.out.versions)


    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
