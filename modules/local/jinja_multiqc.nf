process JINJA_MULTIQC {
    // tag "$meta"
    label 'process_single'

    conda "bioconda::mulled-v2-206e8d7a1e1dc26bc177258a09390bdd3a096228"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-206e8d7a1e1dc26bc177258a09390bdd3a096228' :
        'quay.io/biocontainers/mulled-v2-206e8d7a1e1dc26bc177258a09390bdd3a096228:c967269490aed6e0283bfea2e02d0d77bf259f8d-0' }"

    input:
    path(template_jinja)
    val(meta_id)
    path(image_to_encode_left)
    path(image_to_encode_right)


    output:
    path('*_jinja_out_mqc.html'), emit: rendered
    // path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    jinja_multiqc.py \\
    $template_jinja \\
    $meta_id \\
    $image_to_encode_left \\
    $image_to_encode_right
    """
}
