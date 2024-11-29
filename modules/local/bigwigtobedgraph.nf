process UCSC_BIGWIGTOBEDGRAPH {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bigwigtobedgraph:469--h9b8f530_0' :
        'quay.io/biocontainers/ucsc-bigwigtobedgraph:469--h9b8f530_0' }"

    input:
    tuple val(meta), path(bigwig)
    path  sizes

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '469' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bigWigToBedGraph \\
        $bigwig \\
        ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '469' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
