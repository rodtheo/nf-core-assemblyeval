#!/usr/bin/env python
from jinja2 import Environment, FileSystemLoader
import sys
from pathlib import Path
import base64
import pandas as pd
import csv

def main():
    # python jinja_multiqc.py templa.html $meta.id $img_left $img_right
    meta_id = sys.argv[2]
    lines_to_write = []
    with open(sys.argv[3], "r") as file_fcs_report:
        for line in file_fcs_report.readlines():
            if not line.startswith("#"):
                elems = ["_".join(x.split(" ")) for x in line.rstrip().strip("\n").split("\t")]
                lines_to_write.append(elems)
        # fcs_report = pd.read_csv(file_fcs_report, skiprows=2, names=['seq_id', 'start_pos', 'end_pos',
        #  'seq_len', 'action', 'div', 'agg_cont_cov', 'top_tax_name'])

    mqc_id = "# id: 'contaminant_{}'\n# parent_id: 'contaminant'".format(meta_id)
    mqc_section = "# section_name: 'Contaminant Screening - {}'".format(meta_id)
    header_mqc = """# plot_type: 'table'
# pconfig:
#     namespace: 'FCS'
#     scale: False
# headers:
#     seq_id:
#         description: 'This is a longer hover text for my column'
#     start_pos:
#         description: 'Start coordinate for the identified contamination. If only a portion of the sequence is identified as contaminant, these values indicate the range that should be removed.'
#         format: '{:,.0f}'
#     end_pos:
#         description: 'End coordinate for the identified contamination. If only a portion of the sequence is identified as contaminant, these values indicate the range that should be removed.'
#         format: '{:,.0f}'
#     seq_len:
#         description: 'Length of the entire sequence in Column 1. Only a portion may be identified as contaminant, according to the start_pos and end_pos columns.'
#         format: '{:,.0f}'
#     action:
#         description: 'The recommended action: EXCLUDE, TRIM, FIX, REVIEW, REVIEW_RARE or INFO'
#     div:
#         description: 'The taxonomic division assigned to the contaminant sequence by FCS-GX.'
#     agg_cont_cov:
#         description: 'The percentage alignment coverage for the contaminant in the range indicated by the start_pos and end_pos columns. If the range is composed of multiple contigs separated by gaps, these gaps are ignored for computing coverage. Low coverage values often indicate contamination with novel organisms (e.g., novel genera of bacteria). FCS-GX is tuned for high specificity, even in cases of low reported coverage values.'
#         format: '{:,.0f}'
#     top_tax_name:
#         description: 'The taxonomic name of the top contaminant organism identified by FCS-GX.'"""

        # encoded_string = base64.b64encode(image_file_left.read())
        # data_base64_left = encoded_string.decode()

    # with open(sys.argv[4], "rb") as image_file_right:
    #     encoded_string = base64.b64encode(image_file_right.read())
    #     data_base64_right = encoded_string.decode()

    # asms = [
    #     {"id": meta_id, "img_left": "data:image/png;base64, "+data_base64_left, "img_right": "data:image/png;base64, "+data_base64_right}
    # ]

    # asm = {"id": meta_id, "img_left": "data:image/png;base64, "+data_base64_left, "img_right": "data:image/png;base64, "+data_base64_right}
    row_count = len(lines_to_write)
    if row_count != 0:
        path_template = Path(sys.argv[1])
        path_base = path_template.parent
        # name_template = path_template.name
        # env = Environment(loader=FileSystemLoader(path_base))
        # template = env.get_template("template_mqc.html")
        # template = env.get_template(name_template)

        # filename = path_base/"{}_jinja_out_mqc.html".format(meta_id)
        filename = path_base/"{}_contaminant_out_mqc.txt".format(meta_id)
        # fcs_report.to_csv(filename, sep="\t", index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar="\\")

        with open(filename, mode="w", encoding="utf-8") as output:
            # content = output.read()
            # output.seek(0,0)
            # output.write("# id: 'table_{}'\n".format(meta_id))
            # output.write("# parent_id: 'contaminant'\n")
            output.write(mqc_id+"\n")
            output.write(mqc_section+"\n")
            output.write(header_mqc+"\n")
            # output.write("# plot_type: 'table'\n")
            # output.write("# section_name: 'Test my table'\n")
            output.write(' '.join(['seq_id', 'start_pos', 'end_pos','seq_len', 'action', 'div', 'agg_cont_cov', 'top_tax_name']))
            for l in lines_to_write:
                output.write('\n')
                output.write(' '.join(l))
            


if __name__ == "__main__":
    sys.exit(main())