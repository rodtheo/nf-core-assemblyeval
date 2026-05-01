process SUBSET_REAPR {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/ubuntu:22.04':
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(meta), path(scoreerrors)
    tuple val(meta2), path(strER_CSE_bed), path(strER_CSH_bed), path(locER_CRE_bed), path(locER_CRH_bed)

    output:
    tuple val(meta), path("*_reaper_failure.bed"), emit: reaper_bed_failure
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '22.04' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    zcat $scoreerrors \\
    |sed 's/ /_/g' | grep "FCD_failure" | cut -f 1,4,5,9 > ${prefix}_reaper_base.bed

    cat $strER_CSE_bed $strER_CSH_bed $locER_CRE_bed $locER_CRH_bed | cut -f 1,2,3,5 > ${prefix}_craq.bed

    cat ${prefix}_reaper_base.bed ${prefix}_craq.bed | awk -F"\t" '{print \$1"\t"\$2-1"\t"\$3"\t"\$4}' > ${prefix}_reaper_failure_unsorted.bed

    bedtools sort -i ${prefix}_reaper_failure_unsorted.bed > ${prefix}_reaper_failure.bed

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
    touch ${prefix}_ALEoutput.txt.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UBUNTU: $VERSION)
    END_VERSIONS
    """
}
