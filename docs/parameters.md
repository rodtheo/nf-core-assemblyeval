# nf-core/assemblyeval pipeline parameters

A Nextflow pipeline to evaluate and compare genome assemblies using multiple QC tools, combining results into a unified HTML report.

## Input/output options

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `input` | Path to the input `assemblysheet.yaml` describing assemblies and sequencing reads. | `string` | | True | |
| `outdir` | The output directory where the results will be saved. Use absolute paths for cloud infrastructure. | `string` | | True | |
| `email` | Email address for pipeline completion summary. | `string` | | | |

## Assembly options

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `domain` | Organism domain. Used to configure domain-specific tools. (accepted: `euk`\|`prok`) | `string` | `prok` | | |
| `aligner` | Short-read aligner to use for mapping Illumina reads to assemblies. | `string` | `bwa` | | |
| `long_read_technology` | Minimap2 preset for long-read alignment. Use `map-ont` for Oxford Nanopore reads or `map-pb` for PacBio reads. (accepted: `map-ont`\|`map-pb`) | `string` | `map-ont` | | |

## Completeness options

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `busco` | Use BUSCO instead of COMPLEASM for gene-space completeness assessment. | `boolean` | `false` | | |
| `busco_augustus_species` | Augustus species model to use when running BUSCO in `genome` mode. | `string` | | | |

## Correctness options

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `skip_correctness` | Skip all correctness assessment steps (ALE, REAPR, CRAQ). | `boolean` | `false` | | |
| `reapr_by_chr` | Run REAPR independently per chromosome/contig instead of on the whole assembly. Produces per-contig error profiles at higher computational cost. | `boolean` | `false` | | |

## Contamination options

> [!WARNING]
> The contamination module (FCS-GX) requires a large reference database (~500 GB disk space) and at least **512 GB of RAM** for optimal performance. It is strongly recommended to run this module on an HPC system. Alternatively, pre-filter your assemblies using [NCBI FCS on Galaxy](https://usegalaxy.org/) before running this pipeline.

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `skip_contamination` | Skip the contamination screening module (FCS-adaptor + FCS-GX + FCS-clean). | `boolean` | `true` | | |

## Scoring options

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `yaml_weights_file` | Path to a YAML file specifying the weights used to compute the normalized weighted score combining all evaluation metrics. | `string` | `assets/score_weights.yaml` | | |

## MultiQC options

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` | | | |
| `multiqc_title` | Custom title for the MultiQC report. | `string` | | | |
| `multiqc_logo` | Custom logo file to supply to MultiQC. | `string` | | | |
| `multiqc_methods_description` | Custom MultiQC yaml file containing methods description. | `string` | | | |

## Max resource options

Use these to set upper limits on resources for any single job. These are not the resources used by the pipeline; they are the maximum allowed per process.

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. | `integer` | `16` | | |
| `max_memory` | Maximum amount of memory that can be requested for any single job. | `string` | `128.GB` | | |
| `max_time` | Maximum amount of time that can be requested for any single job. | `string` | `240.h` | | |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `custom_config_version` | Git commit id for Institutional configs. | `string` | `master` | | True |
| `custom_config_base` | Base directory for Institutional configs. | `string` | `https://raw.githubusercontent.com/nf-core/configs/master` | | True |
| `config_profile_name` | Institutional config name. | `string` | | | True |
| `config_profile_description` | Institutional config description. | `string` | | | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
| ---- | ---- | ---- | ---- | ---- | ---- |
| `version` | Display version and exit. | `boolean` | | | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. (accepted: `symlink`\|`rellink`\|`link`\|`copy`\|`copyNoFollow`\|`move`) | `string` | `symlink` | | True |
| `email_on_fail` | Email address for completion summary, only when the pipeline fails. | `string` | | | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` | | | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` | | | True |
| `hook_url` | Incoming hook URL for messaging service. | `string` | | | True |