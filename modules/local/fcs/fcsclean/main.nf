process FCS_CLEAN {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-gx.sif':
        'docker.io/ncbi/fcs-gx:0.4.0' }"

    input:
    tuple val(meta), path(assembly), path(fcs_gx_report)
    // path gxdb

    output:
    tuple val(meta), path("*.contam.fasta")         , emit: fcs_contam_fasta, optional: true
    tuple val(meta), path("*.contamfree.fasta")     , emit: fcs_clean_fasta , optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.4.0'

    """
    /app/bin/gx \\
        clean-genome \\
        --input $assembly \\
        --action-report  $fcs_gx_report \\
        --contam-fasta-out ${prefix}.contam.fasta \\
        --output ${prefix}.contamfree.fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.4.0'

    """
    mkdir -p out
    touch out/${prefix}.contam.fasta
    touch out/${prefix}.contamfree.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """
}
