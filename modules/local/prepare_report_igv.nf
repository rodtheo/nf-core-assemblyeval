process PREPARE_REPORT_IGV {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::mulled-v2-206e8d7a1e1dc26bc177258a09390bdd3a096228"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-206e8d7a1e1dc26bc177258a09390bdd3a096228' :
        'quay.io/biocontainers/mulled-v2-206e8d7a1e1dc26bc177258a09390bdd3a096228:c967269490aed6e0283bfea2e02d0d77bf259f8d-0' }"

    input:
    path(template_jinja)
    val(meta)
    path(reapr_bam)
    path(reapr_bai)
    path(ale_wig_base)
    path(ale_wig_depth)
    path(ale_wig_insert)
    path(ale_wig_kmer)
    path(ale_wig_place)
    path(reapr_score_per_base)


    output:
    tuple val(meta), path('*-track-config.json'), emit: track_config
    tuple val(meta), path(reapr_bam), path(reapr_bai), path(ale_wig_base), path(ale_wig_depth), path(ale_wig_insert), path(ale_wig_kmer), path(ale_wig_place), path(reapr_score_per_base), emit: track_files
    // path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_igv_report_prerunner.py \\
    $template_jinja \\
    $reapr_bam \\
    $reapr_bai \\
    $ale_wig_base \\
    $ale_wig_depth \\
    $ale_wig_insert \\
    $ale_wig_kmer \\
    $ale_wig_place \\
    $reapr_score_per_base \\
    ${prefix}-track-config.json
    """
}
