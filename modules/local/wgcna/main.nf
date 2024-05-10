process WGCNA {
    label "process_low"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "conda-forge::r-base bioconda::r-wgcna bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap"
    
    input:
    path wgcnai

    output:
    output:
    path 'gene_info.csv' , emit: modulesextraction
    path 'modules_dendogram_before_after*'
    path 'modules_distance*'
    path 'module_eigengene_values*'
    path 'module_trait_relationship*'
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    wgcna.r $wgcnai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
