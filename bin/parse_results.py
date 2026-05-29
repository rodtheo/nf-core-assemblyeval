#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import argparse
import logging
import sys
from pathlib import Path
import re
from collections import OrderedDict
from collections import defaultdict
import yaml
import gzip

logger = logging.getLogger()

def normalization_max_best(input_value, input_array, decimal_places=4):
    # min-max normalization
    # print(f"{num:.4f}")  # Output: 3.1416
    res = ((input_value-np.nanmin(input_array))/(np.nanmax(input_array)-np.nanmin(input_array)))
    return f"{res:.{decimal_places}f}"

def normalization_closer(input_value, input_array, target=0, decimal_places=4):
    # normalization for when values are closer to a target (in this case 0) are best
    res = 1 - ((np.abs(input_value - target)) / (np.nanmax(np.abs(input_array - target))))
    return f"{res:.{decimal_places}f}"

def normalization_lower_best(input_value, input_array, decimal_places=4):
    # normalization for when lower values are better
    res = 1 - float(normalization_max_best(input_value, input_array, decimal_places=10))
    return f"{res:.{decimal_places}f}"

def parse_busco(out_short_summary):
    dict_busco = {}
    with open(out_short_summary, 'r') as infbusco:
        for idxline, line in enumerate(infbusco.readlines()):
            # line = line.strip("\n")
            if line.startswith("#"):
                match_db_obj = re.match(r'.+\(.+number of genomes\:\s(\d+),\snumber of BUSCOs\:\s(\d+)\)', line)
                #if match_db_obj:
                #    print(match_db_obj.group(1), match_db_obj.group(2))
            else:
                match_summary_obj = re.match(r'\s+C\:(\d+.\d+)\%\[S\:(\d+.\d)\%,D\:(\d+.\d+)\%\],F\:(\d+.\d+)\%,M:(\d+.\d+)%,n:(\d+)', line)
                if match_summary_obj:
                    dict_busco['pctcomplete'] = match_summary_obj.group(1)
                    dict_busco['pctsingle'] = match_summary_obj.group(2)
                    dict_busco['pctduplicated'] = match_summary_obj.group(3)
                    dict_busco['pctfragmented'] = match_summary_obj.group(4)
                    dict_busco['pctmissing'] = match_summary_obj.group(5)
                    dict_busco['total_busco_searched_genes'] = match_summary_obj.group(6)
                #    print(idxline)
                elif idxline == 9:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['ncomplete'] = n
                elif idxline == 10:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nsingle'] = n
                elif idxline == 11:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nduplicated'] = n
                elif idxline == 12:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nfragmented'] = n
                elif idxline == 13:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nmissing'] = n

        dict_busco['pct_frameshift'] = 'X'
    return(dict_busco)

def parse_compleasm(out_short_summary, compleasm_table):
    dict_busco = {}
    with open(out_short_summary, 'r') as infbusco:
        for idxline, line in enumerate(infbusco.readlines()):
            # line = line.strip("\n")
            if line.startswith("#"):
                match_db_obj = re.match(r'## lineage:\s(\S+)$', line)
                    # S + D
                    
                    # S (Single Copy Complete Genes): The BUSCO genes that can be entirely aligned in the assembly, with only one copy present.
                    
                    # D (Duplicated Complete Genes): The BUSCO genes that can be completely aligned in the assembly, with more than one copy present.
                    
                    # F (Fragmented Genes, subclass 1): The BUSCO genes which only a portion of the gene is present in the assembly, and the rest of the gene cannot be aligned.
                    # PLUS
                    # I (Fragmented Genes, subclass 2): The BUSCO genes in which a section of the gene aligns to one position in the assembly, while the remaining part aligns to another position.
                    
                    # M (Missing Genes): The BUSCO genes with no alignment present in the assembly.
                    
            elif line.startswith("S"):
                match_n = re.match(r'S:(\S+)%,\s(\d+)$', line)
                dict_busco['pctsingle'] = match_n.group(1)
                dict_busco['nsingle'] = match_n.group(2)
            elif line.startswith("D"):
                match_n = re.match(r'D:(\S+)%,\s(\d+)$', line)
                dict_busco['pctduplicated'] = match_n.group(1)
                dict_busco['nduplicated'] = match_n.group(2)
            elif line.startswith("F"):
                match_n = re.match(r'F:(\S+)%,\s(\d+)$', line)
                dict_busco['pctfragmented'] = match_n.group(1)
                dict_busco['nfragmented'] = match_n.group(2)
            elif line.startswith("I"):
                match_n = re.match(r'I:(\S+)%,\s(\d+)$', line)
                dict_busco['pctincomplete'] = match_n.group(1)
                dict_busco['nincomplete'] = match_n.group(2)
            elif line.startswith("M"):
                match_n = re.match(r'M:(\S+)%,\s(\d+)$', line)
                dict_busco['pctmissing'] = match_n.group(1)
                dict_busco['nmissing'] = match_n.group(2)
            else:
                continue
        
    # COMPLETE = S + D
    dict_busco['pctcomplete'] = float(dict_busco['pctsingle']) + float(dict_busco['pctduplicated'])
    dict_busco['ncomplete'] = int(dict_busco['nsingle']) + int(dict_busco['nduplicated'])

    full_table = pd.read_csv(compleasm_table, sep="\t", header=0)
    
    frameshifts = full_table[full_table['Frameshift events'].astype('Int64') != 0]
    
    pct_frameshifts = (frameshifts.shape[0]/full_table.shape[0])*100

    dict_busco['pct_frameshift'] = f"{pct_frameshifts:.2f}"
    dict_busco['n_frameshift'] = frameshifts.shape[0]

    # FRAGMENTED GENES IN BUSCO = 
    # F (Fragmented Genes, subclass 1): The BUSCO genes which only a portion of the gene is present in the assembly, and the rest of the gene cannot be aligned.
    # PLUS
    # I (Fragmented Genes, subclass 2): The BUSCO genes in which a section of the gene aligns to one position in the assembly, while the remaining part aligns to another position.
    # dict_busco['pctfragmented'] = dict_busco['pctfragmented'] + dict_busco['pctincomplete']
    # dict_busco['nfragmented'] = dict_busco['nfragmented'] + dict_busco['nincomplete']

    return(dict_busco)

def parse_results_to_table(genomes_ids, ale_res, reapr_res, busco_re_summary, quast_res, file_out, busco, merfin_qv_res, merfin_comp_res, compleasm_table, craq_aqi_bedgraph, yaml_weights_file):
                # [ ...
                # [ sample_id, [ ale_res ], [reapr_res], [busco_re_summary], [quast_res]],
                # [ sample_id_B, [ ale_res_B ], [reapr_res_B], [busco_re_summary_B], [quast_res_B]]
                #  ... ]
                # ale_res=expand('evaluate_assembly/{sample}/ALEScore_{sample}.finished',sample=samples['sample']),
                # reapr_res=expand('evaluate_assembly/{sample}/REAPR_{sample}.finished',sample=samples['sample']),
                # busco_res=expand('evaluate_assembly/{sample}/busco/BUSCO_{sample}.finished', sample=samples['sample']),
                # quast_res='evaluate_assembly/quast_results/QUAST.OK'
        output = "evaluate_assembly/results.html"
        items = []
        items_classes = []
        my_list = [item.strip() for item in genomes_ids[1:-1].split(',')]
        my_list_ale = [item.strip() for item in ale_res.split(' ')]
        my_list_reapr = [item.strip() for item in reapr_res.split(' ')]
        my_list_busco = [item.strip() for item in busco_re_summary.split(' ')]
        my_list_table_compleasm = [item.strip() for item in compleasm_table.split(' ')]
        my_list_quast = [item.strip() for item in quast_res.split(' ')]
        my_merfin_qv_file = [item.strip() for item in merfin_qv_res.split(' ')]
        my_merfin_completeness_file = [item.strip() for item in merfin_comp_res.split(' ')]
        my_craq_file = [item.strip() for item in craq_aqi_bedgraph.split(' ')]
        for idx, pseudo_samp in enumerate(my_list):
                # samp = pseudo_samp.split("/")[-1]
                samp = pseudo_samp
                dict_sample = OrderedDict()
                dict_classes = OrderedDict()
                dict_sample['name'] = samp
#                        print(pseudo_samp, samp)
                # busco_summary_file = "evaluate_assembly/{}/run_{}/short_summary_{}.txt".format(pseudo_samp, samp, samp)
                if busco:
                    busco_summary_file = my_list_busco[idx]
                    if Path(busco_summary_file).exists():
                            dict_samp = OrderedDict(parse_busco(busco_summary_file))
                    else:
                            dict_samp['ncomplete'] = 'X'
                            dict_samp['pctcomplete'] = 'X'
                            dict_samp['nduplicated'] = 'X'
                            dict_samp['pctduplicated'] = 'X'
                            dict_samp['nfragmented'] = 'X'
                            dict_samp['pctfragmented'] = 'X'
                else:
                    compleasm_summary_file = my_list_busco[idx]
                    compleasm_file = my_list_table_compleasm[idx]
                    if Path(compleasm_summary_file).exists():
                            dict_samp = OrderedDict(parse_compleasm(compleasm_summary_file, compleasm_file))
                # genome_path = samples.loc[samp, "assembly"]
                # genome_name  = genome_path.split('/')[-1]
                # ale_file = "evaluate_assembly/{}/ALEoutput.txt".format(pseudo_samp, genome_name)
                if len(my_list_ale) == len(my_list):
                    ale_file = my_list_ale[idx]
                    if Path(ale_file).exists():
                        # with open(ale_file) as infreapr:
                        with gzip.open(ale_file, 'rt') as infreapr:
                                # for line in infreapr.readlines():
                                for line in infreapr:
                                        match_score = re.match(r'#\sALE_score:\s(-\d+.\d+)', line)
                                        if match_score:
                                            ale_score = float(match_score.group(1))
                                            dict_sample['ale'] = float(ale_score)
                                            dict_sample['neglike_ale'] = -1*float(ale_score)


                # MATCH MERFIN RESULTS - COMPLETENESS
                merfin_file_completeness = my_merfin_completeness_file[idx]
                if Path(merfin_file_completeness).exists():
                    with open(merfin_file_completeness) as infmc:
                            # for line in infmc.readlines():
                            for line in infmc:
                                    match_comp = re.match(r'^COMPLETENESS:\s+(\S+)$', line)
                                    if match_comp:
                                            comp_score = float(match_comp.group(1))
                                            dict_sample['merfin_completeness'] = float(comp_score)*100

                # PARSE CRAQ RESULTS
                craq_report = my_craq_file[idx]
                if Path(craq_report).exists():
                    with open(f"{craq_report}", "r") as craq_report_content:
                        # for line in craq_report_content.readlines():
                        for line in craq_report_content:
                            if line.startswith("Genome\t"):
                                fields = line.split("\t")
                                cre = re.search(r"\(([^)]+)\)", fields[5]).group(1)
                                cse = re.search(r"\(([^)]+)\)", fields[6]).group(1)
                                dict_sample['R-AQI'] = float(cre)
                                dict_sample['S-AQI'] = float(cse)

                # MATCH MERFIN RESULTS - QV*
                merfin_file_qv = my_merfin_qv_file[idx]
                if Path(merfin_file_qv).exists():
                    with open(merfin_file_qv) as infmc:
                            # for line in infmc.readlines():
                            for line in infmc:
                                    match_comp = re.match(r'^Merfin\sQV\*:\s(\S+)$', line)
                                    if match_comp:
                                            qv_score = float(match_comp.group(1))
                    dict_sample['merfin_qv_ast'] = float(qv_score)

#                        ale_file = "evaluate_assembly/{}/ALEoutput.txt".format(pseudo_samp, genome_name)
#                        if path.exists(ale_file):
#                                with open(ale_file, 'r') as inale:
#                                        nline = 0
#                                        for line in inale.readlines():
#                                                if nline == 0:
#                                                        line_list = line.split(" ")
#                                                        dict_sample['ale'] = float(line_list[-1].strip('\n'))
                # BEGIN - PARSING REAPR RESULTS
                # if len(my_list_reapr) == len(my_list):
                reapr_by_genome = defaultdict(list)
                for gen in my_list:
                    for reapr_file in my_list_reapr:
                        if reapr_file.startswith(gen):
                            reapr_by_genome[gen].append(reapr_file)

                
                # for gen, reapr_files in reapr_by_genome.items():
                reapr_percentage_errors_total = 0
                reapr_errors_total = 0
                fcd_errors_total = 0
                low_frag_errors_total = 0
                for reapr_file in reapr_by_genome[my_list[idx]]:
                    print("Parsing REAPR file: ", reapr_file, "for genome: ", my_list[idx])
                # reapr_file = my_list_reapr[idx]
                    if Path(reapr_file).exists():
                        with open(reapr_file) as infreapr:
                            # for line in infreapr.readlines():
                                for line in infreapr:
                                        match_percentage_errors = re.match(r'Error free bases:\s(\d+\.\d+)%', line)
                                        if match_percentage_errors:
                                                reapr_percentage_errors = match_percentage_errors.group(1)
                                        match_errors = re.match(r'^(\d+)\serrors.$', line)
                                        if match_errors:
                                                reapr_errors = match_errors.group(1)
                                        match_fcd_errors = re.match(r'FCD errors within a contig:\s(\d+)', line)
                                        if match_fcd_errors:
                                                fcd_errors = match_fcd_errors.group(1)
                                        match_low_frag_cov = re.match(r'Low fragment coverage within a contig:\s(\d+)', line)
                                        if match_low_frag_cov:
                                                low_frag_errors = match_low_frag_cov.group(1)
                        reapr_percentage_errors_total += float(reapr_percentage_errors)
                        reapr_errors_total += int(reapr_errors)
                        fcd_errors_total += int(fcd_errors)
                        low_frag_errors_total += int(low_frag_errors)
                dict_sample['reapr_percentage_errors'] = reapr_percentage_errors_total
                dict_sample['reapr_total_errors'] = reapr_errors_total
                dict_sample['reapr_fcd'] = fcd_errors_total
                dict_sample['reapr_low'] = low_frag_errors_total
                # END - PARSING REAPR RESULTS
#                        reapr_file = "evaluate_assembly/{}/reapr_results/05.summary.report.txt".format(pseudo_samp)
#                        if path.exists(reapr_file):
#                                with open(reapr_file, 'r') as inreapr:
#                                        for line in inreapr.readlines():
#                                                if line.endswith('errors:\n'):
#                                                        line_list = line.split(" ")
                #                                 dict_sample['reapr'] = line_list[0]
                # parsing quast results
                # df_quast = pd.read_table('evaluate_assembly/quast_results/report.tsv', sep='\t')
                if Path(my_list_quast[idx]).exists():
                    df_quast = pd.read_table(my_list_quast[idx], sep='\t')
                    array_sample = df_quast.iloc[:,1].values
                    print(array_sample)
                    dict_sample['genomesize'] = '{:,}'.format(int(array_sample[6]))
                    dict_sample['contigs']    = '{:,}'.format(int(array_sample[12]))
                    dict_sample['n50']        = '{:,}'.format(int(array_sample[16]))
                    dict_sample['largest']    = '{:,}'.format(int(array_sample[13]))
                    dict_sample['auN']    = '{:,}'.format(float(array_sample[18]))

                    dict_sample['genomesize'] = int(array_sample[6])
                    dict_sample['contigs']    = int(array_sample[12])
                    dict_sample['n50']        = int(array_sample[16])
                    dict_sample['largest']    = int(array_sample[13])
                    dict_sample['auN']    = float(array_sample[18])
                    # concatenating results
                    # dict_sample['genomesize_class'] = "tg-lboi"
                    dict_appended = {**dict_sample, **dict_samp}
                    for k, v in dict_appended.items():
                            key_class = '{}_class'.format(k)
                            dict_classes[key_class] = "tg-lboi"
                    items.append(dict_appended)
                    items_classes.append(dict_classes)
        # print(items)
        # Normalize ALE scores if they are present for all samples
        if len(my_list_ale) == len(my_list):
            ale_scores = []
            reapr_total_errors_scores = []
            reapr_fcd_scores = []
            r_aqi_scores = []
            s_aqi_scores = []
            pct_frameshift_scores = []
            pctcomplete_score = []
            pctmissing_score = []
            merfin_completeness_scores = []
            contigs_scores = []
            n50_scores = []
            auN_scores = []
            largest_scores = []
            merfin_qv_ast_scores = []
            for it in items:
                ale_scores.append(it['neglike_ale'])
                reapr_total_errors_scores.append(it['reapr_total_errors'])
                reapr_fcd_scores.append(it['reapr_fcd'])
                r_aqi_scores.append(it['R-AQI'])
                s_aqi_scores.append(it['S-AQI'])
                pct_frameshift_scores.append(float(it['pct_frameshift']))
                pctcomplete_score.append(float(it['pctcomplete']))
                pctmissing_score.append(float(it['pctmissing']))
                merfin_completeness_scores.append(it['merfin_completeness'])
                contigs_scores.append(it['contigs'])
                n50_scores.append(it['n50'])
                auN_scores.append(it['auN'])
                largest_scores.append(it['largest'])
                merfin_qv_ast_scores.append(it['merfin_qv_ast'])
            ale_scores = np.array(ale_scores)
            reapr_total_errors_scores_array = np.array(reapr_total_errors_scores)
            reapr_fcd_scores_array = np.array(reapr_fcd_scores)
            r_aqi_scores_array = np.array(r_aqi_scores)
            s_aqi_scores_array = np.array(s_aqi_scores)
            pct_frameshift_scores_array = np.array(pct_frameshift_scores)
            pctcomplete_score_array = np.array(pctcomplete_score)
            pctmissing_score_array = np.array(pctmissing_score)
            merfin_completeness_scores_array = np.array(merfin_completeness_scores)
            contigs_scores_array = np.array(contigs_scores)
            n50_scores_array = np.array(n50_scores)
            auN_scores_array = np.array(auN_scores)
            largest_scores_array = np.array(largest_scores)
            merfin_qv_ast_scores_array = np.array(merfin_qv_ast_scores)
            for it in items:
                # min-max normalization
                # it['ale_norm'] = '{:.4f}'.format(((it['neglike_ale']-np.nanmin(ale_scores))/(np.nanmax(ale_scores)-np.nanmin(ale_scores))))
                # normalization for when values are closer to a target (in this case 0) are best
                # it['ale_norm'] = '{:.4f}'.format( 1 - ((np.abs(it['neglike_ale'] - 0)) / (np.abs(np.nanmax(ale_scores) - 0))))

                # normalization for when closer values are better
                it['ale_norm'] = normalization_closer(it['neglike_ale'], ale_scores, target=0, decimal_places=4)
                it['R-AQI_norm'] = normalization_closer(it['R-AQI'], r_aqi_scores_array, target=100, decimal_places=4)
                it['S-AQI_norm'] = normalization_closer(it['S-AQI'], s_aqi_scores_array, target=100, decimal_places=4)
                it['pct_frameshift_norm'] = normalization_closer(float(it['pct_frameshift']), pct_frameshift_scores_array, target=100, decimal_places=4)
                it['pctcomplete_norm'] = normalization_closer(float(it['pctcomplete']), pctcomplete_score_array, target=100, decimal_places=4)
                it['pctmissing_norm'] = normalization_closer(float(it['pctmissing']), pctmissing_score_array, target=0, decimal_places=4)
                it['merfin_completeness_norm'] = normalization_closer(it['merfin_completeness'], merfin_completeness_scores_array, target=100, decimal_places=4)   
                it['merfin_qv_ast_norm'] = normalization_closer(it['merfin_qv_ast'], merfin_qv_ast_scores_array, target=100, decimal_places=4)

                # normalization for when lower values are better
                it['reapr_total_errors_norm'] = normalization_lower_best(it['reapr_total_errors'], reapr_total_errors_scores_array, decimal_places=4)
                it['reapr_fcd_norm'] = normalization_lower_best(it['reapr_fcd'], reapr_fcd_scores_array, decimal_places=4)
                it['contigs_norm'] = normalization_lower_best(it['contigs'], contigs_scores_array, decimal_places=4)
                

                # normalization for when higher values are better
                it['pctcomplete_norm'] = normalization_max_best(it['pctcomplete'], pctcomplete_score_array, decimal_places=4)
                it['n50_norm'] = normalization_max_best(it['n50'], n50_scores_array, decimal_places=4)
                it['auN_norm'] = normalization_max_best(it['auN'], auN_scores_array, decimal_places=4)
                it['largest_norm'] = normalization_max_best(it['largest'], largest_scores_array, decimal_places=4)
            # print(it['ale'])
        #myList = [list(col) for col in zip(*[d.values() for d in items])]
        #myList_argmax = np.argmax(myList, axis=1)
        #print("argmax",myList_argmax)

        #for idx, arg in enumerate(myList_argmax):
#               key_checked = list(dict_appended.keys())[idx]
#                       key_class = "{}_class".format(key_checked)
#                       items_classes[arg][key_class] = "mark"
                #print("K",items[arg][key_checked])
        # res_items = []
        # df = pd.DataFrame(items)
        # df_unnameA = df.drop('name', axis=1)
        # df_unnameA = df_unnameA.drop('ale_norm', axis=1)
        # print(df_unnameA)
        # df_unname = df_unnameA.apply(pd.to_numeric)
        # print("HERE")
        # df_unname['pct_nonduplicated'] = 100. - df_unname['pctduplicated']
        # df_unname['pct_integral'] = 100. - df_unname['pctfragmented']
        # df_unname['pct_found'] = 100. - df_unname['pctmissing']
        # df_sub = df_unname[['genomesize', 'contigs', 'n50', 'largest', 'pctcomplete', 'pct_nonduplicated', 'pct_integral', 'pct_found']]
        # df_sub['name'] = df['name']
        # df_sub_sorted = df_sub.sort_values('name')
        # df_sub_sorted.columns = ['Genome Size (bp)', 'Number of Contigs', 'N50', 'Largest Contig (bp)', 'BUSCO Complete Genes (%)', 'BUSCO Single-Copy Genes (%)', 'BUSCO Non-fragmented Genes (%)', 'BUSCO Found Genes (%)', 'Assembly']
        # df_sub_sorted = df_sub_sorted[['Assembly', 'Genome Size (bp)', 'Number of Contigs', 'N50', 'Largest Contig (bp)', 'BUSCO Complete Genes (%)', 'BUSCO Single-Copy Genes (%)', 'BUSCO Non-fragmented Genes (%)', 'BUSCO Found Genes (%)']]
        # k = df_sub_sorted.style.hide_index().background_gradient('viridis', axis=0, subset=['Genome Size (bp)', 'Number of Contigs', 'N50', 'Largest Contig (bp)', 'BUSCO Complete Genes (%)', 'BUSCO Single-Copy Genes (%)', 'BUSCO Non-fragmented Genes (%)', 'BUSCO Found Genes (%)'])
        # with open('assets/results_heat.html', 'w') as fheat:
        # 	fheat.write(k.render())

# 		for t in list(zip(items, items_classes)):
# 				nd = {**t[0], **t[1]}
# 				res_items.append(nd)
# 		loader = jinja2.FileSystemLoader('template.html')
# 		env = jinja2.Environment(loader=loader)
# 		template = env.get_template('')
# 		output_jinja2 = template.render(items=res_items)
# #               print(output_jinja2)
# 		with open(output[0], 'w') as outfile:
# 				outfile.write(output_jinja2)
        # c = dict([(k,[a[k],b[k]]) for k in items])
        
        with open(yaml_weights_file, 'r') as file:
            data_weights = yaml.safe_load(file)
            print(data_weights)        

        c = pd.DataFrame(items)
        if busco:
            if len(my_list_ale) == len(my_list):
                # {'name': 'Sample_AssemblerA', 'ale': -100937528.909207, 'reapr_total_errors': '402', 'reapr_fcd': '196', 'reapr_low': '206', 'genomesize': 4090859, 
                # 'contigs': 2, 'n50': 4065161, 'largest': 4065161, 'pctcomplete': '6.2', 'pctsingle': '5.7', 'pctduplicated': '0.5', 'pctfragmented': '3.0', 'pctmissing': '90.8', 'total_busco_searched_genes': '758', 'ncomplete': '47', 'nsingle': '43', 'nduplicated': '4', 'nfragmented': '23', 'nmissing': '688'}, {'name': 'Sample_AssemblerB', 'ale': -100937528.909207, 'reapr_total_errors': '402', 'reapr_fcd': '196', 'reapr_low': '206', 'genomesize': 4090859, 'contigs': 2, 'n50': 4065161, 'largest': 4065161, 'pctcomplete': '6.2', 'pctsingle': '5.7', 'pctduplicated': '0.5', 'pctfragmented': '3.0', 'pctmissing': '90.8', 'total_busco_searched_genes': '758', 'ncomplete': '47', 'nsingle': '43', 'nduplicated': '4', 'nfragmented': '23', 'nmissing': '688'}
                dict_names = {'name': 'Assembly', 'ale': 'ALE score (neglog)', 'reapr_total_errors': 'REAPR erros', 'reapr_fcd': 'REAPR fcd', 'reapr_low': 'REAPR low',
                'genomesize': 'Assembly length', 'contigs': 'contigs', 'n50': 'N50', 'largest': 'Largest contig', 'auN': 'auN',
                'pctcomplete': 'BUSCO complete (%)', 'pctsingle': 'BUSCO single (%)', 'pctduplicated': 'BUSCO duplicated (%)', 'pctfragmented':'BUSCO fragmented (%)',
                'pctmissing': 'BUSCO missing (%)',
                'ncomplete': 'BUSCO complete', 'nsingle': 'BUSCO single', 'nduplicated': 'BUSCO duplicated',
                    'nfragmented': 'BUSCO fragmented', 'nmissing': 'BUSCO missing','total_busco_searched_genes': 'BUSCO SEARCHED GENES',
                    'merfin_completeness': 'Merfin Completness','merfin_qv_ast': 'Merfin QV*', 'reapr_percentage_errors': 'Error-free bases (%)'}
            else:
                dict_names = {'name': 'Assembly', 'genomesize': 'Assembly length', 'contigs': 'contigs', 'n50': 'N50', 'largest': 'Largest contig', 'auN': 'auN',
                'pctcomplete': 'BUSCO complete (%)', 'pctsingle': 'BUSCO single (%)', 'pctduplicated': 'BUSCO duplicated (%)', 'pctfragmented':'BUSCO fragmented (%)',
                'pctmissing': 'BUSCO missing (%)',
                'ncomplete': 'BUSCO complete', 'nsingle': 'BUSCO single', 'nduplicated': 'BUSCO duplicated',
                    'nfragmented': 'BUSCO fragmented', 'nmissing': 'BUSCO missing','total_busco_searched_genes': 'BUSCO SEARCHED GENES',
                    'merfin_completeness': 'Merfin Completness','merfin_qv_ast': 'Merfin QV*'}

            c.rename(columns=dict_names, inplace=True)

            # c.columns = ['Assembly', 'ALE score (neglog)', 'REAPR erros', 'REAPR fcd', 'REAPR low',
            #  'Assembly length', 'contigs', 'N50', 'Largest contig', 'BUSCO complete (%)', 'BUSCO single (%)', 'BUSCO duplicated (%)', 'BUSCO fragmented (%)', 'BUSCO missing (%)', 'BUSCO complete', 'BUSCO single', 'BUSCO duplicated', 'BUSCO fragmented', 'BUSCO missing', 'ALE normalized']
        else:
            if len(my_list_ale) == len(my_list):
                dict_names = {'name': 'Assembly', 'ale': 'ALE score (neglog)', 'reapr_total_errors': 'REAPR erros', 'reapr_fcd': 'REAPR fcd', 'reapr_low': 'REAPR low',
             'genomesize': 'Assembly length', 'contigs': 'contigs', 'n50': 'N50', 'largest': 'Largest contig', 'auN': 'auN',
             'pctcomplete': 'COMPLEASM complete (%)', 'pctsingle': 'COMPLEASM single (%)', 'pctduplicated': 'COMPLEASM duplicated (%)', 'pctfragmented':'COMPLEASM fragmented Class I (%)', 'pctincomplete':'COMPLEASM fragmented Class II (%)',
              'pctmissing': 'COMPLEASM missing (%)',
               'ncomplete': 'COMPLEASM complete', 'nsingle': 'COMPLEASM single', 'nduplicated': 'COMPLEASM duplicated', 'nfragmented': 'COMPLEASM fragmented Class I', 
               'nincomplete': 'COMPLEASM fragmented Class II', 'nmissing': 'COMPLEASM missing','total_busco_searched_genes': 'COMPLEASM SEARCHED GENES',
                'merfin_completeness': 'Merfin Completness','merfin_qv_ast': 'Merfin QV*', 'reapr_percentage_errors': 'Error-free bases (%)', 'pct_frameshift': 'Genes with Frameshift (%)', 'n_frameshift': 'Genes with Frameshift (N)'}
            else:
                dict_names = {'name': 'Assembly','genomesize': 'Assembly length', 'contigs': 'contigs', 'n50': 'N50', 'largest': 'Largest contig', 'auN': 'auN',
             'pctcomplete': 'COMPLEASM complete (%)', 'pctsingle': 'COMPLEASM single (%)', 'pctduplicated': 'COMPLEASM duplicated (%)', 'pctfragmented':'COMPLEASM fragmented Class I (%)', 'pctincomplete':'COMPLEASM fragmented Class II (%)',
              'pctmissing': 'COMPLEASM missing (%)',
               'ncomplete': 'COMPLEASM complete', 'nsingle': 'COMPLEASM single', 'nduplicated': 'COMPLEASM duplicated', 'nfragmented': 'COMPLEASM fragmented Class I', 
               'nincomplete': 'COMPLEASM fragmented Class II', 'nmissing': 'COMPLEASM missing','total_busco_searched_genes': 'COMPLEASM SEARCHED GENES',
                'merfin_completeness': 'Merfin Completness','merfin_qv_ast': 'Merfin QV*', 'pct_frameshift': 'Genes with Frameshift (%)', 'n_frameshift': 'Genes with Frameshift (N)'}
            c.rename(columns=dict_names, inplace=True)
        # c.to_excel("evaluate_assembly/results.xlsx", index=False)
        

        c.to_csv(file_out, index=False, sep="\t")

        data_read = pd.read_csv(file_out, sep="\t", header=0, index_col=False)
        if busco:
            if len(my_list_ale) == len(my_list) == len(my_list_reapr):
                greater_better = ['Merfin QV*', 'Error-free bases (%)', 'BUSCO complete (%)', 'BUSCO single (%)', 'auN', 'Merfin Completness']
                smaller_better = ['REAPR erros', 'BUSCO duplicated (%)', 'BUSCO fragmented (%)', 'BUSCO missing (%)']
            else:
                greater_better = ['Merfin QV*', 'BUSCO complete (%)', 'BUSCO single (%)', 'auN', 'Merfin Completness']
                smaller_better = ['BUSCO duplicated (%)', 'BUSCO fragmented (%)', 'BUSCO missing (%)']
        else:
            if len(my_list_ale) == len(my_list) == len(my_list_reapr):
                greater_better = ['Merfin QV*', 'Error-free bases (%)', 'COMPLEASM complete (%)', 'COMPLEASM single (%)', 'auN', 'Merfin Completness']
                smaller_better = ['REAPR erros', 'COMPLEASM duplicated (%)', 'COMPLEASM fragmented Class I (%)', 'COMPLEASM fragmented Class II (%)', 'COMPLEASM missing (%)', 'Genes with Frameshift (%)', 'Genes with Frameshift (N)']
            else:
                greater_better = ['Merfin QV*', 'COMPLEASM complete (%)', 'COMPLEASM single (%)', 'auN', 'Merfin Completness']
                smaller_better = ['COMPLEASM duplicated (%)', 'COMPLEASM fragmented Class I (%)', 'COMPLEASM fragmented Class II (%)', 'COMPLEASM missing (%)', 'Genes with Frameshift (%)', 'Genes with Frameshift (N)']
        complement = data_read[smaller_better].max(axis=0)

        print("SMAAAAAAAAAAAAAAAAALLLLLLLLLERRRRRRRRRRRRRR", data_read[smaller_better])
        print(complement)
        scaler_sb = MinMaxScaler()
        subset_sb = complement.sub(data_read[smaller_better])
        scaler_sb.fit(subset_sb)
        scaler_sb_df = scaler_sb.transform(subset_sb)
        # scaler_sb_df

        scaler = MinMaxScaler()
        scaler.fit(data_read[greater_better])
        scaler_df = scaler.transform(data_read[greater_better])

        score = np.mean(np.concatenate((scaler_df, scaler_sb_df), axis=1), axis=1)*100

        # data_read['Score'] = score

        column_types = {
            'correctness': {'cols': [], 'weights': []},
            'contiguity': {'cols': [], 'weights': []},
            'completeness': {'cols': [], 'weights': []},
            'score': {'cols': [], 'weights': []}
        }
        
        for type_error, weight in data_weights.items():
            if type_error != 'score':
                for metric in weight:
                    metric_key, metric_value = list(metric.items())[0]
                    column_types[type_error]['cols'].append(metric_key+"_norm")
                    column_types[type_error]['weights'].append(float(metric_value))
                if np.sum(column_types[type_error]['weights']) != 1:
                    print(f"Error: The weights for metric {type_error} do not sum up to 1. Please check the YAML weights file.")
                    break
            else:
                for metric in weight:
                    print("METRIC:", list(metric.items()))
                    metric_key, metric_value = list(metric.items())[0]
                    print(metric_key, metric_value)
                    column_types['score']['cols'].append("score_"+metric_key)
                    column_types['score']['weights'].append(float(metric_value))
                if np.sum(column_types['score']['weights']) != 1:
                    print(f"Error: The weights for {type_error} do not sum up to 1. Please check the YAML weights file.")
                    break
        
        data_read['score_correctness'] = np.average(data_read[column_types['correctness']['cols']], weights=column_types['correctness']['weights'], axis=1)*100
        data_read['score_contiguity'] = np.average(data_read[column_types['contiguity']['cols']], weights=column_types['contiguity']['weights'], axis=1)*100
        data_read['score_completeness'] = np.average(data_read[column_types['completeness']['cols']], weights=column_types['completeness']['weights'], axis=1)*100
        # print(data_read[column_types['score']['cols']])
        data_read['Score'] = np.average(data_read[column_types['score']['cols']], weights=column_types['score']['weights'], axis=1)

        data_read = data_read.drop(columns=["Error-free bases (%)", "REAPR erros", "REAPR low"] + column_types['correctness']['cols'] + column_types['contiguity']['cols'] + column_types['completeness']['cols'])

        data_read.to_csv(file_out, index=False, sep="\t")

        with open(file_out, 'r') as original:
            data = original.read()
        with open(file_out, 'w') as modified:
            modified.write("# id: 'assemblies_table'\n")
            modified.write("# section_anchor: 'assemblies_table'\n")
            modified.write("# description: 'a custom text introduction (a few sentences) for this section'\n")
            comments_table_columns = """# pconfig:
#     namespace: 'Cust Data'
# headers:
#     'Merfin Completeness':
#         title: 'My Column'
#         description: 'This is a longer hover text for my column'
#     'REAPR erros':
#         title: 'Second Column'
#         description: 'Hover description text'
#         format: '{:,.0f}'
"""
            # modified.write(comments_table_columns)
            modified.write("# section_name: 'Assembly Stats'\n# plot_type: 'table'\n" + data)
        print("Success ! The results summary table has been written ! \n The results can be view in:\n \t- Excel format in file evaluate_assembly/results.xlsx \n \t- HTML format in file evaluate_assembly/results.html \n \t- HTML heatmap in file evaluate_assembly/results_head.html \n \t- CSV format in file evaluate_assembly/results.csv")


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Transform record fasta header if it contains things like trailing whitespace or characters |':- that could break the REAPR execution.",
        epilog="Example: python check_header_fasta.py sample.fasta sample_valid.fasta",
    )
    parser.add_argument(
        "-g",
        "--genomes_ids",
        # nargs='+',
        metavar="GENOMES_IDS",
        # type=Path,
        help="Fasta input.",
    )
    parser.add_argument(
        "-a",
        "--ale_res",
        # nargs='+',
        metavar="ALE_RES",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "-r",
        "--reapr_res",
        # nargs='+',
        metavar="REAPR_RES",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "-b",
        "--busco_re_summary",
        # nargs='+',
        metavar="BUSCO_RE_SUMMARY",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "--craq_aqi_bedgraph",
        # nargs='+',
        metavar="craq_aqi_bedgraph",
        # type=Path,
        help="Craq output report.",
    )

    parser.add_argument(
        "-q",
        "--quast_res",
        # nargs='+',
        metavar="QUAST_RES",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "--merfin_qv_res",
        # nargs='+',
        metavar="MERFIN_QV_RES",
        # type=Path,
        help="Merfin hist log outputs.",
    )

    parser.add_argument(
        "--merfin_comp_res",
        # nargs='+',
        metavar="MERFIN_COMP_RES",
        # type=Path,
        help="Merfin completness output.",
    )

    parser.add_argument(
        "--compleasm_table",
        # nargs='+',
        metavar="COMPLEASM_TABLE",
        # type=Path,
        help="Compleasm full table input.",
    )

    parser.add_argument(
        "-f",
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "--busco",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        metavar="BUSCO",
        help="Active script to parse BUSCO results instead of COMPLEASM.",
    )

    parser.add_argument(
        "--yaml-weights-file",
        type=Path,
        help="Path to the YAML file containing score weights.",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    # if not args.file_out.is_file():
    #     logger.error(f"The given input file {args.file_in} was not found!")
    #     sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    # check_fasta(args.file_in, args.file_out)
    print(args)
    parse_results_to_table(args.genomes_ids, args.ale_res, args.reapr_res, args.busco_re_summary, args.quast_res, args.file_out, args.busco, args.merfin_qv_res, args.merfin_comp_res, args.compleasm_table, args.craq_aqi_bedgraph, args.yaml_weights_file)


if __name__ == "__main__":
    sys.exit(main())
        