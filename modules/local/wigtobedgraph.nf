process UCSC_WIGTOBEDGRAPH_GZ {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-wigtobigwig:447--h2a80c09_1' :
        'rodtheo/parse_assembly_eval:1.5.0' }"

    input:
    tuple val(meta), path(wig)
    path sizes

    output:
    // tuple val(meta), path("*.bedGraph.gz"), emit: bedgraph_gz
    tuple val(meta), path("*.bedGraph.gz"), path("*.bedGraph.gz.tbi"), emit: gz_tbi
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    // bedtools sort -i ${prefix}.bedGraph | bgzip  --threads ${task.cpus} -c $args3 > ${prefix}.bedGraph.gz
    """
    wigToBigWig \\
        $args \\
        $wig \\
        $sizes \\
        ${prefix}.bw

    bigWigToBedGraph \\
        $args2 \\
        ${prefix}.bw \\
        ${prefix}.bedGraph

    bgzip  --threads ${task.cpus} -c $args3 ${prefix}.bedGraph > ${prefix}.bedGraph.gz
    tabix --threads ${task.cpus} -0 -p bed $args4 ${prefix}.bedGraph.gz

    rm -rf ${prefix}.bw ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}