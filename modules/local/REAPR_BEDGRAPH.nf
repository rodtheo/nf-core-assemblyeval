process REAPER_BEDGRAPH {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/ubuntu:22.04':
        'quay.io/biocontainers/ubuntu:22.04' }"

    input:
    tuple val(meta), path(scoreerrorsperbase)

    output:
    tuple val(meta), path("*_reaper_score_per_base.bedGraph"), emit: reapr_per_base_bedgraph
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '22.04' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    zcat ${scoreerrorsperbase} \\
    | awk -F"\t" 'BEGIN{ OFS="\t" } {\$4=\$2-1} {print \$1,\$4,\$2,\$3}' > ${prefix}_reaper_score_per_base.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UBUNTU: $VERSION)
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '22.04' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_reaper_score_per_base.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UBUNTU: $VERSION)
    END_VERSIONS
    """
}
