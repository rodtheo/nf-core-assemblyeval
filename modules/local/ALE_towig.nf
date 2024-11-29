process ALE_TO_WIGGLE {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ale:20180904--py27ha92aebf_0':
        'rodtheo/genomics:eval_assem_ale_reapr' }"

    input:
    tuple val(meta), path(aleout)

    output:
    tuple val(meta), path("*_ALEoutput.txt.wig"), emit: ale_wig_base
    tuple val(meta), path("*_ALEoutput.txt-depth.wig"), emit: ale_wig_depth
    tuple val(meta), path("*_ALEoutput.txt-insert.wig"), emit: ale_wig_insert
    tuple val(meta), path("*_ALEoutput.txt-kmer.wig"), emit: ale_wig_kmer
    tuple val(meta), path("*_ALEoutput.txt-place.wig"), emit: ale_wig_place
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20180904' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    ale2wiggle.py \\
        $aleout

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ale: $VERSION
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20180904' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_ALEoutput.txt.wig
    touch ${prefix}_ALEoutput.txt-depth.wig
    touch ${prefix}_ALEoutput.txt-insert.wig
    touch ${prefix}_ALEoutput.txt-kmer.wig
    touch ${prefix}_ALEoutput.txt-place.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ale: $VERSION
    END_VERSIONS
    """
}
