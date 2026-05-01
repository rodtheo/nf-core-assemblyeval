process CRAQ {
    tag "$meta.id"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/craq:1.0.9--hdfd78af_0':
        'quay.io/biocontainers/craq:1.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(asm)
    tuple path(bam_shortreads), path(bai)
    path bam_longreads
    // val mode                              // Required:    One of genome, proteins, or transcriptome
    // val lineage                           // Required:    lineage to check against, "auto" enables --auto-lineage instead
    // path busco_lineages_path              // Recommended: path to busco lineages - downloads if not set
    // path config_file                      // Optional:    busco configuration file
    

    output:
    tuple val(meta), path("*-craq"), emit: craq_dir
    tuple val(meta), path("*-craq/runAQI_out/out_final.Report"), emit: summary_report
    tuple val(meta), path("*-craq/runAQI_out/locER_out/out_final.CRE.bed"), emit: locER_bed_CRE
    tuple val(meta), path("*-craq/runAQI_out/locER_out/out_final.CRH.bed"), emit: locER_bed_CRH
    tuple val(meta), path("*-craq/runAQI_out/locER_out/ambiguous.RE.RH"), emit: locER_bed_ambiguous
    tuple val(meta), path("*-craq/runAQI_out/strER_out/out_final.CSE.bed"), emit: strER_bed_CRE
    tuple val(meta), path("*-craq/runAQI_out/strER_out/out_final.CSH.bed"), emit: strER_bed_CRH
    tuple val(meta), path("*-craq/runAQI_out/strER_out/ambiguous.SE.SH"), emit: strER_bed_ambiguous
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    def craq_longreads = bam_longreads ? "-sms $bam_longreads" : ''
    // def busco_config = config_file ? "--config $config_file" : ''
    // def compleasm_lineage = lineage.equals('auto') ? '--autolineage' : "-l ${lineage}"
    // def busco_lineage_dir = busco_lineages_path ? "--download_path ${busco_lineages_path}" : ''
    """
    craq \\
        -g $asm \\
        --thread $task.cpus \\
        -sms $bam_shortreads \\
        $craq_longreads \\
        $args \\
        --output_dir ${prefix}-craq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        craq: 1.0.9
    END_VERSIONS
    """
}
