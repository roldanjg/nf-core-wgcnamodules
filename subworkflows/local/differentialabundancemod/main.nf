
def exp_meta = [ "id": params.study_name  ]
matrix_file = file(params.salmon_dir+"/salmon.merged.gene_counts.tsv", checkIfExists: true)

ch_in_raw = Channel.of([ exp_meta, matrix_file])

ch_input = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ])

ch_transcript_lengths = Channel.of([ exp_meta, file(params.salmon_dir+"/salmon.merged.gene_lengths.tsv", checkIfExists: true)]).first()
ch_control_features = [[],[]]

report_file = file(params.report_file, checkIfExists: true)
logo_file = file(params.logo_file, checkIfExists: true)
css_file = file(params.css_file, checkIfExists: true)
citations_file = file(params.citations_file, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM} from '../../../modules/nf-core/deseq2/differential/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../../modules/nf-core/custom/dumpsoftwareversions/main'
include { CUSTOM_MATRIXFILTER } from '../../../modules/nf-core/custom/matrixfilter/main'
include { RMARKDOWNNOTEBOOK } from '../../../modules/nf-core/rmarkdownnotebook/main'
include { SHINYNGS_STATICDIFFERENTIAL as PLOT_DIFFERENTIAL  } from '../../../modules/nf-core/shinyngs/staticdifferential/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY    } from '../../../modules/nf-core/shinyngs/staticexploratory/main'
include { SHINYNGS_VALIDATEFOMCOMPONENTS  as VALIDATOR} from '../../../modules/nf-core/shinyngs/validatefomcomponents/main'
include { ZIP as MAKE_REPORT_BUNDLE                         } from '../../../modules/nf-core/zip/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// process checalgo (
//     input: estose
//     output:
//     estose
// )
workflow DIFFERENTIALABUNDANCE {

    // Set up some basic variables
    ch_versions = Channel.empty()
    // Channel for the contrasts file
    ch_contrasts_file = Channel.from([[exp_meta, file(params.contrast)]])

    

    // // Otherwise we can just use the matrix input; save it to the workdir so that it does not
    // // just appear wherever the user runs the pipeline
    matrix_as_anno_filename = "${workflow.workDir}/matrix_as_anno.${matrix_file.getExtension()}"

    ch_features_matrix = ch_in_raw

    ch_features = ch_features_matrix
        .map{ meta, matrix ->
            // println meta
            // println matrix
            matrix.copyTo(matrix_as_anno_filename)
            [ meta, file(matrix_as_anno_filename) ]
        }
    ch_matrices_for_validation = ch_in_raw
    
    VALIDATOR(
        ch_input.join(ch_matrices_for_validation),
        ch_features,
        ch_contrasts_file
    )
    ch_raw = VALIDATOR.out.assays
    ch_matrix_for_differential = ch_raw
    ch_contrasts = VALIDATOR.out.contrasts
        .map{it[1]}
        .splitCsv ( header:true, sep:'\t' )
        .map{
            it.blocking = it.blocking.replace('NA', '')
            if (!it.id){
                it.id = it.values().join('_')
            }
            tuple(it, it.variable, it.reference, it.target)
        }
    CUSTOM_MATRIXFILTER(
        ch_matrix_for_differential,
        VALIDATOR.out.sample_meta
    )
    ch_samples_and_matrix = VALIDATOR.out.sample_meta
        .join(CUSTOM_MATRIXFILTER.out.filtered)     // -> meta, samplesheet, filtered matrix
        .first()
    
    DESEQ2_NORM (
            ch_contrasts.first(),
            ch_samples_and_matrix,
            ch_control_features,
            ch_transcript_lengths
        )

        // Run the DESeq differential module, which doesn't take the feature
        // annotations

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_samples_and_matrix,
        ch_control_features,
        ch_transcript_lengths
    )
    // DESEQ2_DIFFERENTIAL.out.results.flatten().view()
    // Let's make the simplifying assumption that the processed matrices from
    // the DESeq runs are the same across contrasts. We run the DESeq process
    // with matrices once for each contrast because DESeqDataSetFromMatrix()
    // takes the model, and the model can vary between contrasts if the
    // blocking factors included differ. But the normalised and
    // variance-stabilised matrices are not (IIUC) impacted by the model.

    ch_norm = DESEQ2_NORM.out.normalised_counts
    ch_differential = DESEQ2_DIFFERENTIAL.out.results
    ch_model = DESEQ2_DIFFERENTIAL.out.model

    ch_versions = ch_versions
        .mix(DESEQ2_NORM.out.versions)
    ch_versions = ch_versions
        .mix(DESEQ2_DIFFERENTIAL.out.versions)

    ch_processed_matrices = ch_norm
    if ('rlog' in params.deseq2_vs_method){
        ch_processed_matrices = ch_processed_matrices.join(DESEQ2_NORM.out.rlog_counts)
    }
    if ('vst' in params.deseq2_vs_method){
        ch_processed_matrices = ch_processed_matrices.join(DESEQ2_NORM.out.vst_counts)
    }
    ch_processed_matrices = ch_processed_matrices
        .map{ it.tail() }
    ch_contrast_variables = ch_contrasts
    .map{
        [ "id": it[1] ]
    }
    .unique()
    ch_mat = ch_raw.combine(ch_processed_matrices)

    ch_all_matrices = VALIDATOR.out.sample_meta                // meta, samples
    .join(VALIDATOR.out.feature_meta)                       // meta, samples, features
    .join(ch_mat)                                           // meta, samples, features, raw, norm (or just norm)
    .map{
        tuple(it[0], it[1], it[2], it[3..it.size()-1])
    }
    .first()

    // PLOT_EXPLORATORY(
    //     ch_contrast_variables
    //         .combine(ch_all_matrices.map{ it.tail() })
    // )

    // // Differential analysis using the results of DESeq2

    // PLOT_DIFFERENTIAL(
    //     ch_differential,
    //     ch_all_matrices
    // )
    
    // // Gather software versions

    // ch_versions = ch_versions
    //     .mix(VALIDATOR.out.versions)
    //     .mix(PLOT_EXPLORATORY.out.versions)
    //     .mix(PLOT_DIFFERENTIAL.out.versions)

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    // // Generate a list of files that will be used by the markdown report

    // ch_report_file = Channel.from(report_file)
    //     .map{ tuple(exp_meta, it) }

    // ch_logo_file = Channel.from(logo_file)
    // ch_css_file = Channel.from(css_file)
    // ch_citations_file = Channel.from(citations_file)

    // ch_report_input_files = ch_all_matrices
    //     .map{ it.tail() }
    //     .map{it.flatten()}
    //     .combine(VALIDATOR.out.contrasts.map{it.tail()})
    //     .combine(CUSTOM_DUMPSOFTWAREVERSIONS.out.yml)
    //     .combine(ch_logo_file)
    //     .combine(ch_css_file)
    //     .combine(ch_citations_file)
    //     .combine(ch_differential.map{it[1]}.toList())
    //     .combine(ch_model.map{it[1]}.toList())

    //  // Make a params list - starting with the input matrices and the relevant
    // // params to use in reporting

    // def report_file_names = [ 'observations', 'features' ] +
    //     params.exploratory_assay_names.split(',').collect { "${it}_matrix".toString() } +
    //     [ 'contrasts_file', 'versions_file', 'logo', 'css', 'citations' ]
    // def params_pattern = ~/^(report|study|observations|features|filtering|exploratory|differential|deseq2|gsea).*/
    // ch_report_params = ch_report_input_files
    //     .map{
    //         params.findAll{ k,v -> k.matches(params_pattern) } +
    //         [report_file_names, it.collect{ f -> f.name}].transpose().collectEntries()
    //     }

    // // Render the final report
    // RMARKDOWNNOTEBOOK(
    //     ch_report_file,
    //     ch_report_params,
    //     ch_report_input_files
    // )
    emit:
    tables         = ch_differential
    versions       = ch_versions                                // channel: [ versions.yml ]
}