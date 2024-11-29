#!/usr/bin/env python
import sys
from pathlib import Path
import json
from jinja2 import Environment, FileSystemLoader


def main():

	path_template = Path(sys.argv[1])
	path_base = path_template.parent

	reapr_bam = sys.argv[2]

	reapr_bai = sys.argv[3]

	ale_wig_base = sys.argv[4]

	ale_wig_depth = sys.argv[5]

	ale_wig_insert = sys.argv[6]

	ale_wig_kmer = sys.argv[7]

	ale_wig_place = sys.argv[8]

	reapr_score_per_base = sys.argv[9]

	file_out = sys.argv[10]

	env = Environment(loader=FileSystemLoader(searchpath=path_base))
	# env = Environment(loader=PackageLoader('app', 'templates'))
	env.filters['jsonify'] = json.dumps

	# Template file at ./app/templates/example.json
	template = env.get_template('track-config.json')

	page = {
		'reapr_bam': reapr_bam,
		'reapr_bai': reapr_bai,
		'ale_wig_base': ale_wig_base,
		'ale_wig_depth': ale_wig_depth,
		'ale_wig_insert': ale_wig_insert,
		'ale_wig_kmer': ale_wig_kmer,
		'ale_wig_place': ale_wig_place,
		'reapr_score_per_base': reapr_score_per_base
	}

	# print(template.render(page=page))
	# filename = path_base/"{}-track-config.json".format(meta_id)
	with open(file_out, mode="w", encoding="utf-8") as output:
		output.write(template.render(page=page))

if __name__ == "__main__":
    sys.exit(main())