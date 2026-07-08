Here's a simple README for the script:

---

# parse_results.py

A Python script to parse and aggregate genome assembly evaluation metrics from multiple tools, normalize them, and output a summary table compatible with [MultiQC](https://multiqc.info/).

## Description

This script collects results from several assembly QC tools for multiple assemblies, computes per-metric normalizations, applies user-defined weights, and produces a scored summary table ready to be rendered as a MultiQC custom table.

## Supported Tools

| Tool | Metrics Extracted |
|------|------------------|
| **BUSCO** | Complete, single-copy, duplicated, fragmented, missing gene counts and percentages |
| **Compleasm** | Same as BUSCO + frameshift event rate |
| **QUAST** | Assembly length, contig count, N50, largest contig, auN |
| **ALE** | Assembly likelihood score |
| **REAPR** | Total errors, FCD errors, low fragment coverage errors |
| **Merfin** | QV* score, k-mer completeness |
| **CRAQ** | R-AQI and S-AQI scores |

## Requirements

```bash
pip install pandas numpy scikit-learn pyyaml
```

## Usage

```bash
python parse_results.py \
  -g "[sample_A, sample_B]" \
  -a "ale_A.txt.gz ale_B.txt.gz" \
  -r "reapr_A.txt reapr_B.txt" \
  -b "busco_summary_A.txt busco_summary_B.txt" \
  -q "quast_A.tsv quast_B.tsv" \
  --merfin_qv_res "merfin_qv_A.txt merfin_qv_B.txt" \
  --merfin_comp_res "merfin_comp_A.txt merfin_comp_B.txt" \
  --compleasm_table "compleasm_A.tsv compleasm_B.tsv" \
  --craq_aqi_bedgraph "craq_A.txt craq_B.txt" \
  --yaml-weights-file weights.yaml \
  -f results/assembly_stats.tsv \
  [--busco true|false]
```

### Key Arguments

| Argument | Description |
|----------|-------------|
| `-g` / `--genomes_ids` | Comma-separated list of sample names (e.g. `[sampleA, sampleB]`) |
| `-a` / `--ale_res` | Space-separated ALE output files (gzipped) |
| `-r` / `--reapr_res` | Space-separated REAPR summary report files |
| `-b` / `--busco_re_summary` | Space-separated BUSCO or Compleasm short summary files |
| `-q` / `--quast_res` | Space-separated QUAST TSV report files |
| `--merfin_qv_res` | Space-separated Merfin QV* output files |
| `--merfin_comp_res` | Space-separated Merfin completeness output files |
| `--compleasm_table` | Space-separated Compleasm full table files |
| `--craq_aqi_bedgraph` | Space-separated CRAQ report files |
| `-f` / `--file_out` | Output file path (TSV format) |
| `--busco` | Use BUSCO parser instead of Compleasm (default: `false`) |
| `--yaml-weights-file` | Path to YAML file defining metric weights for scoring |
| `-l` / `--log-level` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`) |

## Weights YAML File

The scoring system requires a YAML file defining weights for three quality categories — **correctness**, **contiguity**, and **completeness** — plus a final **score** that combines the three. Weights within each category must sum to `1.0`.

```yaml
correctness:
  - ale_norm: 0.5
  - reapr_total_errors_norm: 0.5
contiguity:
  - n50_norm: 0.4
  - auN_norm: 0.3
  - contigs_norm: 0.3
completeness:
  - pctcomplete_norm: 0.5
  - merfin_completeness_norm: 0.5
score:
  - score_correctness: 0.4
  - score_contiguity: 0.3
  - score_completeness: 0.3
```

## Output

The script produces a TSV file with:
- Raw and normalized metrics per assembly
- Per-category scores (`score_correctness`, `score_contiguity`, `score_completeness`)
- A final composite `Score`
- MultiQC-compatible headers prepended to the file for direct use as a [MultiQC custom data table](https://multiqc.info/docs/development/custom_content/)

## Normalization Strategy

| Strategy | Applied to |
|----------|-----------|
| **Higher is better** (min-max) | N50, auN, largest contig, BUSCO/Compleasm complete (%) |
| **Lower is better** (inverted min-max) | Contig count, REAPR errors |
| **Closer to target is better** | ALE score (→ 0), QV* (→ 100), CRAQ AQI (→ 100), Merfin completeness (→ 100) |