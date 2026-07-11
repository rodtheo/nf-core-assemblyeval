<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-assemblyeval_logo_dark.png">
    <img alt="nf-core/assemblyeval" src="docs/images/nf-core-assemblyeval_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/assemblyeval/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/assemblyeval/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/assemblyeval/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/assemblyeval/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/assemblyeval/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/assemblyeval)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23assemblyeval-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/assemblyeval)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**Assemblyeval** accepts genome assemblies (FASTA), paired-end Illumina reads, and long reads (ONT or PacBio) via a YAML samplesheet, optionally cleaning them of contamination before evaluation. The pipeline systematically assesses contiguity (QUAST), completeness (COMPLEASM/BUSCO, Merfin), and correctness (ALE, REAPR, CRAQ) — combining read-alignment-based and k-mer-based evidence into a single normalized weighted score. Final outputs include a comparative MultiQC report ranking assemblies across all metrics and an IGV-based interactive report for manual inspection of putative structural errors.


<!-- ![pipeline](assets/assembly-eval-pipeline-v5.png) -->

```mermaid
flowchart TD
    INPUT[/"YAML file describing inputs (relevant metadata about the sample, assemblies in FASTA, short and long reads in FASTQ)"/]:::yellow
    INPUT --> PREPARE["PREPARE_INPUT\nassemblies · illumina · long reads"]

    PREPARE -->|"Illumina reads"| FASTQC["🔘 FASTQC"]
    PREPARE --> SKIP_CONT{"--skip_contamination?"}:::white

    SKIP_CONT -->|"no"| FCS_A["🔴 FCS_FCSADAPTOR\nadapter screen"]
    FCS_A --> FCS_G["🔴 FCS_FCSGX\nforeign genome screen"]
    FCS_G --> FCS_C["🔴 FCS_CLEAN\nremove contaminants"]
    FCS_C --> ASMS(("Assemblies (ASMs)")):::white
    SKIP_CONT -->|"yes"| ASMS

    ASMS --> COMP_BRANCH{"--busco?"}:::white
    COMP_BRANCH -->|"true"| BUSCO["🟢 BUSCO_BUSCO"]
    COMP_BRANCH -->|"false (default)"| COMPLEASM["🟢🔵 COMPLEASM"]

    ASMS --> QUAST["🟡 QUAST\ncontiguity stats"]

    ASMS -->|"+ short reads"| BWA["BWAMEM2_INDEX\nBWAMEM2_MEM\nread alignment"]
	ASMS -->|"+ long reads"| MINIMAP2["MINIMAP2 alignment"]
    BWA --> SAMTOOLS["SAMTOOLS_SORT\nSAMTOOLS_INDEX"]

    SAMTOOLS --> REAPR_BRANCH{"--reapr_by_chr?"}:::white
    ASMS --> REAPR_BRANCH
    REAPR_BRANCH -->|"false (default)"| HDR["HEADER_FASTA_REAPR\nsanitize FASTA headers"]
    HDR --> REAPR["🔵 REAPR\nwhole-genome error detection"]
    REAPR_BRANCH -->|"true"| HDR_CHR["HEADER_FASTA_REAPR\nsanitize + split by chr"]
    HDR_CHR --> BED_INT["BEDTOOLS_INTERSECT\nsubset BAM per chromosome"]
    BED_INT --> REAPR_CHR["🔵 REAPR_BY_CHR\nper-chromosome error detection"]

    SAMTOOLS -->|"BAM + BAI + ASM"| ALE["🔵 ALE\nassembly likelihood score"]
    SAMTOOLS -->|"BAM + BAI + ASM"| CRAQ["🔵 CRAQ\nstructural error detection"]
	MINIMAP2 --> CRAQ

    ASMS -->|"+ Illumina reads"| MERYL["🟣 MERYL_COUNT\nMERYL_HISTOGRAM"]
    MERYL --> GS2["🟣 GENOMESCOPE2\ngenome profiling"]
    GS2 --> MERFIN["🟢🔵 MERFIN\nQV* + completeness"]

    ALE --> ALE_WIG["ALE_TO_WIGGLE\n5 WIG tracks: base · depth · insert · kmer · place"]
    ALE_WIG --> WIG_BG["UCSC_WIGTOBEDGRAPH + TABIX"]

    REAPR --> REAPR_BG["REAPER_BEDGRAPH + TABIX\nper-base REAPR score"]
    REAPR_CHR --> REAPR_BG

    CRAQ --> CRAQ_TAB["BEDTOOLS_SORT + TABIX\nregional AQI bedgraph"]

    WIG_BG --> SUB_REAPR["🔘 SUBSET_REAPR\nfilter REAPR failure regions"]
    REAPR_BG --> SUB_REAPR
    CRAQ_TAB --> SUB_REAPR

    SUB_REAPR --> PREP_IGV["🔘 PREPARE_REPORT_IGV\ntrack config JSON"]
    PREP_IGV --> IGVREPORTS["🔘 IGVREPORTS"]

    BUSCO --> PARSE["🔘 PARSE_RESULTS\nnormalize metrics · compute weighted score\nscore_weights.yaml"]
    COMPLEASM --> PARSE
    QUAST --> PARSE
    ALE --> PARSE
    REAPR --> PARSE
    REAPR_CHR --> PARSE
    CRAQ --> PARSE
    MERFIN --> PARSE

    FASTQC --> MULTIQC["🔘 MULTIQC"]
    PARSE --> MULTIQC

    MULTIQC --> OUT1[/"multiqc/\nMain HTML Report"/]:::orange
    IGVREPORTS --> OUT2[/"igvreports/\nPer-assembly IGV Report"/]:::orange

    subgraph Legenda
        direction LR
        L1["🟢 Completeness"]
        L2["🔴 Contamination"]
        L3["🟡 Contiguity"]
        L4["🔵 Correctness"]
        L5["🟣 KMER_MODULE"]
        L6["🔘 Output"]
    end

    classDef white fill:#ffffff,stroke:#000000,color:#000000
    classDef yellow fill:#FFFF00,stroke:#000000,color:#000000
    classDef orange fill:#EDC001,stroke:#000000,color:#000000
```

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Parse and prepare input samplesheet (`PREPARE_INPUT`: assemblies · Illumina reads · long reads)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Contamination screening and removal (_Optional_, `--skip_contamination`):
   1. Adapter/vector screen ([`FCS-adaptor`](https://github.com/ncbi/fcs))
   2. Foreign genome screen ([`FCS-GX`](https://github.com/ncbi/fcs))
   3. Generate cleaned assembly ([`FCS-CLEAN`](https://github.com/ncbi/fcs))
4. Contiguity assessment ([`QUAST`](https://github.com/ablab/quast))
5. Completeness assessment (choice of):
   1. [`COMPLEASM`](https://github.com/huangnengCSU/compleasm) _(default)_
   2. [`BUSCO`](https://busco.ezlab.org/) _(optional, `--busco`)_
6. K-mer profiling and QV estimation:
   1. K-mer counting ([`Meryl`](https://github.com/marbl/meryl))
   2. Genome profiling ([`GenomeScope2`](https://github.com/tbenavi1/genomescope2.0))
   3. K-mer completeness and QV* score ([`Merfin`](https://github.com/arangrhie/merfin))
7. Short-read alignment to assemblies ([`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2), [`SAMtools`](https://www.htslib.org/))
8. Long-read alignment to assemblies ([`Minimap2`](https://github.com/lh3/minimap2))
Here's the corrected item 9:

9. Correctness assessment:
   1. Assembly likelihood score ([`ALE`](https://github.com/sc932/ALE)) — evaluates short-read alignments using a Bayesian framework to compute a global log-probability score that the assembly is correct, integrating four per-nucleotide metrics: alignment concordance, insert size consistency, depth of coverage uniformity, and k-mer frequency
   2. Structural error detection ([`CRAQ`](https://github.com/bioinfo-biols/CRAQ)) — uses clipped (split-mapped) short- and long-read alignments to identify regional and structural errors, reporting a Regional Assembly Quality Index (R-AQI) and a Structural Assembly Quality Index (S-AQI) on a 0–100 scale
   3. Error region detection via Fragment Coverage Distribution ([`REAPR`](https://www.sanger.ac.uk/tool/reapr/)) — evaluates the discrepancy between observed paired-read mapping coverage and theoretical expectations (FCD error), scoring each base from 0 to 1 based on read orientation, insert size consistency, and soft-clipping evidence; two run modes are supported:
      1. Whole-assembly mode _(default)_ — runs REAPR across the entire assembly at once
      2. Per-contig mode _(optional, `--reapr_by_chr true`)_ — splits the assembly by contig/chromosome, subsets the BAM accordingly with [`BEDTools`](https://github.com/arq5x/bedtools2/), and runs REAPR independently per contig; recommended for large assemblies (>500 Mbp)
10. Convert correctness scores to indexed BEDGraph tracks ([`UCSC wigToBedGraph`](http://hgdownload.soe.ucsc.edu/admin/exe/), [`Tabix`](https://www.htslib.org/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
11. Filter and prepare candidate misassembly regions for visualization (`SUBSET_REAPR`, `PREPARE_REPORT_IGV`)
12. Generate interactive per-assembly error report ([`IGV-Reports`](https://github.com/igvteam/igv-reports))
13. Aggregate metrics, compute weighted ranking scores and present QC ([`MultiQC`](http://multiqc.info/))

> [!WARNING]
> The contamination module requires the pre-download of the **FCS-GX database**, which demands a substantial amount of disk space and a **minimum of 512 GB of RAM** for full-performance execution. Running with less memory is possible but may result in up to a **10,000× performance penalty**. For this reason, we strongly recommend running the contamination steps (`--skip_contamination false`) only on **HPC (High-Performance Computing)** environments.
>
> If you do not have access to HPC infrastructure, you can run **FCS-GX online via Galaxy** — [ncbi_fcs_gx on Galaxy](https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fncbi_fcs_gx%2Fncbi_fcs_gx%2F0.5.0%2Bgalaxy0&version=latest) — and then use the resulting decontaminated FASTA file as input to **nf-core/assemblyeval** with `--skip_contamination`.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare an `assemblysheet.yaml` with the following fields:

- **`metadata`**: organism-level metadata shared across assemblies:
  - `id`: sample or organism name
  - `kmer_size`: k-mer size for tools like `Meryl`
  - `ploidy`: ploidy level
  - `organism_domain`: `euk` or `prok`
  - `taxid`: NCBI taxonomy ID
  - `busco_lineages`: BUSCO lineage(s) for completeness evaluation
- **`assembly`**: one or more assemblies to evaluate. Each entry requires a unique `id` (no spaces) and `pri_asm` with the path to the FASTA file.
- **`illumina`**: paired-end Illumina reads (`read1` and `read2`).
- **`ont`**: long reads file path. For PacBio data, pass `--map-pb` when running the pipeline.


Below is an example [`assemblysheet.yaml`](assets/chla_test_input_channels.yaml) evaluating four assemblies of *C. trachomatis*:

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
    illumina:
      - read1: "./data_test/mason_R1_001.fastq.gz"
        read2: "./data_test/mason_R2_001.fastq.gz"
    ont:
      - reads: "./data_test/ont_reads_pbsim3.fq.gz"
```

To get started, clone the repository, enter the folder, and test the pipeline using the bundled C. trachomatis test dataset:

```bash
git clone https://github.com/rodtheo/nf-assemblyeval.git

cd nf-assemblyeval

export NXF_VER=25.04.3

nextflow run main.nf \
   -profile docker \
   --max_cpus <N> \
   --input assets/chla_test_input_channels.yaml \
   --outdir <OUTDIR>
```

Once the pipeline completes, the main outputs are:
- **General evaluation report**: `<OUTDIR>/multiqc/Report-for-General-Evaluation-of-Assemblies_multiqc_report.html`
- **Per-assembly error inspection**: `<OUTDIR>/igvreports/`

> [!WARNING]
> Provide pipeline parameters via the CLI or `--params-file`. The `-c` Nextflow option handles configuration only and **cannot** set pipeline parameters. See [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details, see the [usage documentation](https://nf-co.re/assemblyeval/usage) and [parameter documentation](https://nf-co.re/assemblyeval/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://rodtheo.github.io/simposio-biomol-2024/simulations/Report-for-General-Evaluation-of-Assemblies_multiqc_report.html) and the [IGV report](https://rodtheo.github.io/simposio-biomol-2024/simulations/Ctrachomatis_ctracho_5inversions_default_report_mqc.html).

After the pipeline finishes, the main outputs are self-contained HTML located at `<OUTDIR>/multiqc/` and `<OUTDIR>/igvreport/`.

For more details about the output files and reports, please refer to the
[output documentation]().

## Disclaimer

This is not an official nf-core community pipeline (at least, not yet), although it strives to follow the community's recommended standards as closely as possible.

## Credits

nf-core/assemblyeval was originally written by Rodrigo Theodoro Rocha.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#assemblyeval` channel](https://nfcore.slack.com/channels/assemblyeval) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/assemblyeval for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
