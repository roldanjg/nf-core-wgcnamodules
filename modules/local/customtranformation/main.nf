process CUSTOM_FILES_TRANSFORMATION {
    label 'process_single'

    conda "conda-forge::python=3.11.8 conda-forge::pandas=2.2.1"

    input:
    path input_wgcna
    path contrast_wgcna
    path tpms_wgcna
    path genes

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    if (params.diffgenes == true) {
        """
        customETL.py \\
            --input_wgcna $input_wgcna \\
            --contrast_wgcna $contrast_wgcna \\
            --tpms_wgcna $tpms_wgcna \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    } else {
        """
        customETL.py \\
            --input_wgcna $input_wgcna \\
            --contrast_wgcna $contrast_wgcna \\
            --tpms_wgcna $tpms_wgcna \\
            --genes $genes \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
}
