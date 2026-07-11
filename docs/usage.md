# nf-core/assemblyeval: Usage<!-- omit in toc -->

* [Assemblysheet input](#assemblysheet-input)

* [External databases](#external-databases)

  * [NCBI FCS-GX database](#ncbi-fcs-gx-database)

  * [BUSCO lineages](#busco-lineages)

* [Other parameters](#other-parameters)

  * [Assembly options](#assembly-options)

  * [Completeness options](#completeness-options)

  * [Correctness options](#correctness-options)

  * [Scoring options](#scoring-options)

* [Minimum System Requirements](#minimum-system-requirements)

* [Running the pipeline](#running-the-pipeline)

  * [Updating the pipeline](#updating-the-pipeline)

  * [Reproducibility](#reproducibility)

* [Core Nextflow arguments](#core-nextflow-arguments)

  * `[-profile](#-profile)`

  * `[-resume](#-resume)`

  * `[-c](#-c)`

* [Custom configuration](#custom-configuration)

  * [Resource requests](#resource-requests)

  * [Custom Containers](#custom-containers)

  * [Custom Tool Arguments](#custom-tool-arguments)

  * [nf-core/configs](#nf-coreconfigs)

* [Running in the background](#running-in-the-background)

## Assemblysheet input

You will need to create an `assemblysheet.yaml` file with information about the assemblies you would like to evaluate before running the pipeline. Use the `--input` parameter to specify its location. It must be a valid YAML file. An [example assemblysheet](../assets/chla_test_input_channels.yaml) has been provided with the pipeline.

The file has the following top-level structure:

```yaml
samples:
  - metadata:
      ...
    assembly:
      ...
    illumina:
      ...
    ont:
      ...
```

Its fields are:

* **`metadata`** _(required)_: Organism-level metadata shared across all assemblies in the sheet:

  * `id`: A unique name for the sample or organism (e.g. `Ctrachomatis`).

  * `kmer_size`: K-mer length used by `Meryl` for k-mer counting (e.g. `21`).

  * `ploidy`: Ploidy level of the organism (e.g. `1` for haploid, `2` for diploid).

  * `organism_domain`: Either `euk` (eukaryote) or `prok` (prokaryote). Used to configure domain-specific tool parameters.

  * `taxid`: NCBI taxonomy ID for the organism as a quoted string (e.g. `"315277"`). A taxonomy ID can be obtained by searching a _Genus species_ at [https://www.ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

  * `busco_lineages`: A YAML list of one or more BUSCO lineage dataset names (e.g. `chlamydiae_odb10`). To select a lineage, see [https://busco.ezlab.org/list_of_lineages.html](https://busco.ezlab.org/list_of_lineages.html).

* **`assembly`** _(required)_: A list of assemblies to evaluate and compare. Each entry requires:

  * `id`: A unique identifier for the assembly (no spaces). This tag is used throughout the pipeline and in the final reports.

  * `pri_asm`: Path to the primary assembly FASTA file.

* **`illumina`** _(required)_: Paired-end Illumina short reads used for alignment-based correctness assessment:

  * `read1`: Path to the forward (R1) FASTQ file.

  * `read2`: Path to the reverse (R2) FASTQ file.

* **`ont`** _(required for CRAQ long-read correctness assessment)_: Long reads for alignment with Minimap2:

  * `reads`: Path to the long reads FASTQ file. If using PacBio reads instead of ONT, set `--long_read_technology map-pb` at runtime.

A complete example evaluating four assemblies of _Chlamydia trachomatis_:

```yaml
samples:
  - metadata:
      id: Ctrachomatis
      kmer_size: 21
      ploidy: 1
      organism_domain: prok
      taxid: "315277"
      busco_lineages:
        - "chlamydiae_odb10"
    assembly:
      - id: ref
        pri_asm: "./data_test/GCF_000012125.1_ASM1212v1_genomic.fasta"
      - id: ctracho_5inversions_default
        pri_asm: "./data_test/ctracho_5inversions_default.simseq.genome.fa"
      - id: ctracho_10persnp_default
        pri_asm: "./data_test/ctracho_10persnp_default.simseq.genome.fa"
      - id: ctracho_5inversions_split
        pri_asm: "./data_test/Ctrachomatis_ctracho_5inversions_default_cleaned_breaked.fasta"
    illumina:
      - read1: "./data_test/mason_R1_001.fastq.gz"
        read2: "./data_test/mason_R2_001.fastq.gz"
    ont:
      - reads: "./data_test/ont_reads_pbsim3.fq.gz"
```

## External databases

### NCBI FCS-GX database

> \[!WARNING\]\
> The FCS-GX contamination module requires a large reference database (\~500 GB of disk space) and at least **512 GB of RAM** for optimal performance. It is strongly recommended to run this module on an HPC system. Alternatively, pre-filter your assemblies using [NCBI FCS on Galaxy](https://usegalaxy.org/) before running this pipeline with `--skip_contamination` (the default).

If contamination screening is enabled by setting `--skip_contamination false`, the FCS-GX database must be available on the host filesystem. Setup instructions and download scripts are available at [https://github.com/ncbi/fcs/wiki/FCS-GX](https://github.com/ncbi/fcs/wiki/FCS-GX). The database directory must contain the following files:

```
all.assemblies.tsv
all.blast_div.tsv.gz
all.gxi
all.gxs
all.manifest
all.meta.jsonl
all.README.txt
all.seq_info.tsv.gz
all.taxa.tsv
```

The path to this directory must be made accessible to the pipeline's container runtime (Docker/Singularity) by mounting it appropriately.

### BUSCO lineages

BUSCO lineage datasets are downloaded automatically during pipeline execution. Lineages are specified per-sample in the `assemblysheet.yaml` under `busco_lineages`. Available lineages are listed at [https://busco.ezlab.org/list_of_lineages.html](https://busco.ezlab.org/list_of_lineages.html).

By default, the pipeline uses **COMPLEASM** as a faster BUSCO-compatible alternative, which also requires a BUSCO lineage to be specified. Note that **lineage version `odb12` must be used** with COMPLEASM. BUSCO is only used when `--busco true` is set. See [Completeness options](#completeness-options) below.

## Other parameters

This section provides additional information for selected parameters. For the full parameter list, see [parameters.md](./parameters.md).

### Assembly options

* **`--domain`**: Sets the organism domain (`euk` or `prok`). This affects tool configurations that behave differently for eukaryotes and prokaryotes. Default: `prok`.

* **`--aligner`**: Short-read aligner used to map Illumina reads to the assemblies prior to correctness assessment. Default: `bwa` (BWA-MEM2).

* **`--long_read_technology`**: Minimap2 preset for long-read alignment. Use `map-ont` for Oxford Nanopore Technology reads (default) or `map-pb` for PacBio reads. This affects both CRAQ correctness assessment and any long-read-based analyses.

### Completeness options

* **`--busco`**: When set to `true`, uses BUSCO instead of COMPLEASM for gene-space completeness assessment. COMPLEASM is used by default as it is faster and produces BUSCO-compatible results. Default: `false`.

* **`--busco_augustus_species`**: Augustus species model for BUSCO when running in `genome` mode. Only relevant when `--busco true` is set. Refer to the [Augustus species list](https://github.com/Gaius-Augustus/Augustus/tree/master/config/species) for available models.

### Correctness options

* **`--skip_correctness`**: Skips all three correctness assessment tools (ALE, REAPR, and CRAQ). Useful when short-read or long-read data are not available, or to reduce runtime. Default: `false`.

* **`--reapr_by_chr`**: When set to `true`, REAPR is run independently per chromosome/contig rather than on the whole assembly at once. This produces per-contig error profiles and is more sensitive for detecting localised errors, but requires significantly more compute time and memory for assemblies with many sequences. Default: `false` (whole-assembly mode).

### Scoring options

* **`--yaml_weights_file`**: Path to a YAML file specifying the weights used to compute a normalized weighted score combining all evaluation metrics. The default file is `assets/score_weights.yaml`. Users can supply a custom file to tune the relative importance of each metric.

  The YAML file has three metric group sections (`correctness`, `completeness`, `contiguity`), each containing per-metric weights that must sum to `1.0`. A final `score` section defines the contribution of each group to the overall assembly score. The default weights are:

  ```yaml
  correctness:
    - ale: 0.1
    - reapr_fcd: 0.4
    - R-AQI: 0.1
    - S-AQI: 0.1
    - pct_frameshift: 0.1
    - merfin_qv_ast: 0.2
  completeness:
    - pctcomplete: 0.5
    - pctmissing: 0.4
    - merfin_completeness: 0.1
  contiguity:
    - contigs: 0.2
    - n50: 0.4
    - auN: 0.2
    - largest: 0.2
  score:
    - correctness: 0.55
    - completeness: 0.40
    - contiguity: 0.05
  ```

  These defaults prioritise correctness (55%) and completeness (40%) over contiguity (5%), reflecting the expectation that modern long-read assemblers already produce highly contiguous assemblies.

## Minimum System Requirements

Resource requirements depend heavily on genome size and the number of assemblies being compared. The pipeline has been tested in the following scenarios (contamination module excluded):

* **Small prokaryotic/fungal genome (\~30 MB, e.g. _Aspergillus_ sp.)**: A run comparing **\[TODO: number of assemblies\]** assemblies was executed on a desktop workstation with an AMD Ryzen 9 7950X (16-core / 32-thread) and 92 GB of RAM. Peak memory usage reached **\[TODO: peak memory\]** and the total runtime was approximately **\[TODO: runtime\]**.

* **Medium eukaryotic genome (\~700 MB)**: A comparable run took approximately **12 hours** and consumed between **150–200 GB** of disk space for intermediate working files.

For the optional contamination module (FCS-GX), the minimum requirements are:

* **FCS-GX**: 1 CPU + **512 GB RAM**. Performance degrades substantially with less memory, as portions of the database must be read from disk. An HPC system is strongly recommended for this module.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run rodtheo/nf-assemblyeval \
    --input ./assemblysheet.yaml \
    --outdir ./results \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work/               # Directory containing the Nextflow working files
<OUTDIR>/           # Finished results in the specified output location (--outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, e.g. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> \[!WARNING\]\
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in YAML format:

```bash
nextflow run rodtheo/nf-assemblyeval -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: "./assemblysheet.yaml"
outdir: "./results/"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available — even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, regularly update the cached version:

```bash
nextflow pull rodtheo/nf-assemblyeval
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used. First, go to the [nf-core/assemblyeval releases page](https://github.com/rodtheo/nf-assemblyeval/releases) and find the latest pipeline version — numeric only (e.g. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) — e.g. `-r 1.0.0`.

<!-- TODO: update release page URL once the pipeline is published -->

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> \[!TIP\]\
> If you wish to share such a params file (e.g. as supplementary material for academic publications), make sure to NOT include cluster-specific paths to files or institution-specific profiles.

## Core Nextflow arguments

> \[!NOTE\]\
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) — see below.

> \[!IMPORTANT\]\
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility. When this is not possible, Conda is also supported.

Note that multiple profiles can be loaded, for example: `-profile test,docker` — the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines.

* `test`

  * A profile with a complete configuration for automated testing using the bundled _C. trachomatis_ dataset

  * Includes links to test data so needs no other parameters

* `docker`

  * A generic configuration profile to be used with [Docker](https://docker.com/)

* `singularity`

  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)

* `podman`

  * A generic configuration profile to be used with [Podman](https://podman.io/)

* `shifter`

  * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)

* `charliecloud`

  * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)

* `apptainer`

  * A generic configuration profile to be used with [Apptainer](https://apptainer.org/)

* `wave`

  * A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `≥24.03.0-edge`).

* `conda`

  * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most pipeline steps, if the job exits with an error it will automatically be resubmitted with higher resource requests (2× original, then 3× original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline step for a particular tool. By default, nf-core pipelines use containers and software from the [Biocontainers](https://biocontainers.pro/) or [Bioconda](https://bioconda.github.io/) projects. To use a different container, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default. To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

If you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings, it may be a good idea to request that your custom config file is uploaded to the `[nf-core/configs](https://github.com/nf-core/configs)` git repository. Before doing so, please test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in `[nf-core/configs/docs](https://github.com/nf-core/configs/tree/master/docs)`), and amend `[nfcore_custom.config](https://github.com/nf-core/configs/blob/master/nfcore_custom.config)` to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the `[#assemblyeval](https://nfcore.slack.com/channels/assemblyeval)`[ channel](https://nfcore.slack.com/channels/assemblyeval).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or a similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run Nextflow within a cluster job submitted to your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machine can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```