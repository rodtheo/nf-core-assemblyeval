#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import logging
import sys
import yaml
from pathlib import Path

logger = logging.getLogger()


# ─────────────────────────────────────────────
# Normalization functions (from parse_results.py)
# ─────────────────────────────────────────────

def normalization_max_best(input_value, input_array, decimal_places=4):
    res = (input_value - np.nanmin(input_array)) / (np.nanmax(input_array) - np.nanmin(input_array))
    return f"{res:.{decimal_places}f}"


def normalization_closer(input_value, input_array, target=0, decimal_places=4):
    res = 1 - (np.abs(input_value - target) / np.nanmax(np.abs(input_array - target)))
    return f"{res:.{decimal_places}f}"


def normalization_lower_best(input_value, input_array, decimal_places=4):
    res = 1 - float(normalization_max_best(input_value, input_array, decimal_places=10))
    return f"{res:.{decimal_places}f}"


# ─────────────────────────────────────────────
# Column name -> display name for final output
# ─────────────────────────────────────────────

DISPLAY_NAMES = {
    'name':                  'Assembly',
    'ale':                   'ALE score',
    'merfin_completeness':   'Merfin Completness',
    'R-AQI':                 'R-AQI',
    'S-AQI':                 'S-AQI',
    'merfin_qv_ast':         'Merfin QV*',
    'reapr_fcd':             'REAPR fcd',
    'genomesize':            'Assembly length',
    'contigs':               'contigs',
    'n50':                   'N50',
    'largest':               'Largest contig',
    'auN':                   'auN',
    'pctsingle':             'COMPLEASM single (%)',
    'nsingle':               'COMPLEASM single',
    'pctduplicated':         'COMPLEASM duplicated (%)',
    'nduplicated':           'COMPLEASM duplicated',
    'pctfragmented':         'COMPLEASM fragmented Class I (%)',
    'nfragmented':           'COMPLEASM fragmented Class I',
    'pctincomplete':         'COMPLEASM fragmented Class II (%)',
    'nincomplete':           'COMPLEASM fragmented Class II',
    'pctmissing':            'COMPLEASM missing (%)',
    'nmissing':              'COMPLEASM missing',
    'pctcomplete':           'COMPLEASM complete (%)',
    'ncomplete':             'COMPLEASM complete',
    'pct_frameshift':        'Genes with Frameshift (%)',
    'n_frameshift':          'Genes with Frameshift (N)',
    'score_correctness':     'score_correctness',
    'score_contiguity':      'score_contiguity',
    'score_completeness':    'score_completeness',
    'Score':                 'Score',
}

# Columns to drop from final output (internal/intermediate)
DROP_COLS = {
    'neglike_ale', 'reapr_percentage_errors', 'reapr_low', 'reapr_total_errors',
}

NORM_SUFFIX = '_norm'


def normalize_items(items: list) -> list:
    """Re-compute all _norm fields jointly across all items."""

    def _arr(key):
        vals = []
        for it in items:
            v = it.get(key)
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                vals.append(np.nan)
        return np.array(vals, dtype=float)

    ale_scores                 = _arr('neglike_ale')
    reapr_total_errors_scores  = _arr('reapr_total_errors')
    reapr_fcd_scores           = _arr('reapr_fcd')
    r_aqi_scores               = _arr('R-AQI')
    s_aqi_scores               = _arr('S-AQI')
    pct_frameshift_scores      = _arr('pct_frameshift')
    pctcomplete_score          = _arr('pctcomplete')
    pctmissing_score           = _arr('pctmissing')
    merfin_completeness_scores = _arr('merfin_completeness')
    contigs_scores             = _arr('contigs')
    n50_scores                 = _arr('n50')
    auN_scores                 = _arr('auN')
    largest_scores             = _arr('largest')
    merfin_qv_ast_scores       = _arr('merfin_qv_ast')

    for it in items:
        nale = it.get('neglike_ale')
        valid_ale = ale_scores[~np.isnan(ale_scores)]
        if nale is not None and not (isinstance(nale, float) and np.isnan(nale)) and len(valid_ale) > 0:
            it['ale_norm'] = normalization_closer(float(nale), ale_scores, target=0, decimal_places=4)
        else:
            it['ale_norm'] = np.nan

        it['R-AQI_norm']               = normalization_closer(float(it['R-AQI']),              r_aqi_scores,                target=100, decimal_places=4)
        it['S-AQI_norm']               = normalization_closer(float(it['S-AQI']),              s_aqi_scores,                target=100, decimal_places=4)
        it['pct_frameshift_norm']      = normalization_closer(float(it['pct_frameshift']),      pct_frameshift_scores,       target=100, decimal_places=4)
        it['pctmissing_norm']          = normalization_closer(float(it['pctmissing']),          pctmissing_score,            target=0,   decimal_places=4)
        it['merfin_completeness_norm'] = normalization_closer(float(it['merfin_completeness']), merfin_completeness_scores,  target=100, decimal_places=4)
        it['merfin_qv_ast_norm']       = normalization_closer(float(it['merfin_qv_ast']),       merfin_qv_ast_scores,        target=100, decimal_places=4)

        it['reapr_total_errors_norm']  = normalization_lower_best(float(it['reapr_total_errors']), reapr_total_errors_scores, decimal_places=4)
        it['reapr_fcd_norm']           = normalization_lower_best(float(it['reapr_fcd']),           reapr_fcd_scores,          decimal_places=4)
        it['contigs_norm']             = normalization_lower_best(float(it['contigs']),             contigs_scores,            decimal_places=4)

        it['pctcomplete_norm']         = normalization_max_best(float(it['pctcomplete']), pctcomplete_score, decimal_places=4)
        it['n50_norm']                 = normalization_max_best(float(it['n50']),         n50_scores,        decimal_places=4)
        it['auN_norm']                 = normalization_max_best(float(it['auN']),         auN_scores,        decimal_places=4)
        it['largest_norm']             = normalization_max_best(float(it['largest']),     largest_scores,    decimal_places=4)

    return items


def compute_scores(data: pd.DataFrame, data_weights: dict) -> pd.DataFrame:
    """Apply YAML-defined weighted scores."""
    column_types = {
        'correctness':  {'cols': [], 'weights': []},
        'contiguity':   {'cols': [], 'weights': []},
        'completeness': {'cols': [], 'weights': []},
        'score':        {'cols': [], 'weights': []},
    }

    for type_error, weight in data_weights.items():
        if type_error != 'score':
            for metric in weight:
                metric_key, metric_value = list(metric.items())[0]
                column_types[type_error]['cols'].append(metric_key + '_norm')
                column_types[type_error]['weights'].append(float(metric_value))
            if round(np.sum(column_types[type_error]['weights']), 6) != 1.0:
                logger.error(f"Weights for '{type_error}' do not sum to 1.")
                sys.exit(1)
        else:
            for metric in weight:
                metric_key, metric_value = list(metric.items())[0]
                column_types['score']['cols'].append('score_' + metric_key)
                column_types['score']['weights'].append(float(metric_value))
            if round(np.sum(column_types['score']['weights']), 6) != 1.0:
                logger.error(f"Weights for 'score' do not sum to 1.")
                sys.exit(1)

    data['score_correctness']  = np.average(data[column_types['correctness']['cols']].astype(float),  weights=column_types['correctness']['weights'],  axis=1) * 100
    data['score_contiguity']   = np.average(data[column_types['contiguity']['cols']].astype(float),   weights=column_types['contiguity']['weights'],   axis=1) * 100
    data['score_completeness'] = np.average(data[column_types['completeness']['cols']].astype(float), weights=column_types['completeness']['weights'], axis=1) * 100
    data['Score']              = np.average(data[column_types['score']['cols']].astype(float),         weights=column_types['score']['weights'],         axis=1)

    norm_cols = (
        column_types['correctness']['cols'] +
        column_types['contiguity']['cols'] +
        column_types['completeness']['cols']
    )
    data = data.drop(columns=[c for c in norm_cols if c in data.columns])
    return data


def write_multiqc_table(data: pd.DataFrame, file_out: Path):
    data.to_csv(file_out, index=False, sep="\t")
    with open(file_out, 'r') as f:
        content = f.read()
    with open(file_out, 'w') as f:
        f.write("# id: 'assemblies_table'\n")
        f.write("# section_anchor: 'assemblies_table'\n")
        f.write("# description: 'a custom text introduction (a few sentences) for this section'\n")
        f.write("# section_name: 'Assembly Stats'\n")
        f.write("# plot_type: 'table'\n")
        f.write(content)
    logger.info(f"Output written to: {file_out}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Read *.base.tsv outputs from parse_results.py, re-normalize jointly "
            "across all inputs, and output a single MultiQC-compatible TSV."
        ),
        epilog="Example: python concatenate_results_tables.py -i run1.base.tsv run2.base.tsv --yaml-weights-file weights.yaml -f merged.tsv",
    )
    parser.add_argument("-i", "--input_tables", nargs="+", required=True, metavar="BASE_TSV", type=Path,
                        help="One or more *.base.tsv files produced by parse_results.py.")
    parser.add_argument("--yaml-weights-file", required=True, type=Path,
                        help="YAML file with score weights.")
    parser.add_argument("-f", "--file_out", required=True, type=Path,
                        help="Output TSV file path.")
    parser.add_argument("-l", "--log-level", choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
                        default="WARNING", help="Logging verbosity (default: WARNING).")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # 1. Read and concatenate
    dfs = []
    for tsv_path in args.input_tables:
        if not tsv_path.exists():
            logger.error(f"File not found: {tsv_path}")
            sys.exit(1)
        logger.info(f"Reading: {tsv_path}")
        dfs.append(pd.read_csv(tsv_path, sep="\t", header=0, index_col=False))

    combined_df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Total rows after concatenation: {len(combined_df)}")

    # 2. Drop existing _norm columns (will be recomputed jointly)
    old_norm_cols = [c for c in combined_df.columns if c.endswith(NORM_SUFFIX)]
    combined_df = combined_df.drop(columns=old_norm_cols)
    logger.info(f"Dropped {len(old_norm_cols)} pre-existing _norm column(s).")

    # 3. Convert to list of dicts for normalization
    items = combined_df.to_dict(orient='records')

    # 4. Re-normalize jointly
    items = normalize_items(items)

    # 5. Merge _norm columns back into DataFrame
    norm_keys = [k for k in items[0].keys() if k.endswith(NORM_SUFFIX)]
    for i, row_dict in enumerate(items):
        for nk in norm_keys:
            combined_df.at[i, nk] = row_dict[nk]

    # 6. Drop internal columns
    combined_df = combined_df.drop(columns=[c for c in DROP_COLS if c in combined_df.columns])

    # 7. Rename to display names
    combined_df = combined_df.rename(columns=DISPLAY_NAMES)

    # 8. Load weights and compute scores
    with open(args.yaml_weights_file, 'r') as f:
        data_weights = yaml.safe_load(f)

    combined_df = compute_scores(combined_df, data_weights)

    # 9. Write output
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    write_multiqc_table(combined_df, args.file_out)
    print(f"Done. Results written to: {args.file_out}")


if __name__ == "__main__":
    sys.exit(main())