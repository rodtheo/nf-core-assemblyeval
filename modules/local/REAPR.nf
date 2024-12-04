process REAPR {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.5.0--pyhdfd78af_0':
        'rodtheo/genomics:eval_assem_ale_reapr' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/reapr:v1.0.18dfsg-4-deb_cv1':
    //     'biocontainers/reapr:v1.0.18dfsg-4-deb_cv1' }"
        

    input:
    tuple val(meta), path(asm), path(bam)       // Required:    One genome assembly to evaluate path(bam)            // Required:    Corresponding aligned file
    

    output:
    tuple val(meta), path("*-REAPR"), emit: reapr_dir
    tuple val(meta), path("*-REAPR/03.score.errors.gff.gz"), emit:reapr_score_errors
    tuple val(meta), path("*-REAPR/00.in.bam"), emit:reapr_bam
    tuple val(meta), path("*-REAPR/00.in.bam.bai"), emit:reapr_bai
    tuple val(meta), path("*-REAPR/05.summary.report.txt"), emit:reapr_summary
    tuple val(meta), path("*-REAPR/03.score.per_base.gz"), emit: reapr_score_per_base
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    def VERSION = '1.0.18' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    reapr pipeline \\
        $asm \\
        $bam \\
        ${prefix}-REAPR

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        REAPR: $VERSION)
    END_VERSIONS
    """
}
