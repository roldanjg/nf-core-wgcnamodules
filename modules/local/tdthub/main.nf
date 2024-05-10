process TDTHUBMODULES {
    label 'process_single'

    conda "conda-forge::python=3.11.8 conda-forge::pandas=2.2.1"

    input:
    path genes_modules

    output:
    path "*.csv"
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:


    """
    customTDTHub.py --genes_modules $genes_modules

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """

}
