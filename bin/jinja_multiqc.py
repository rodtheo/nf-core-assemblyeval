#!/usr/bin/env python
from jinja2 import Environment, FileSystemLoader
import sys
from pathlib import Path
import base64

def main():
    # python jinja_multiqc.py templa.html $meta.id $img_left $img_right
    meta_id = sys.argv[2]

    with open(sys.argv[3], "rb") as image_file_left:
        encoded_string = base64.b64encode(image_file_left.read())
        data_base64_left = encoded_string.decode()

    with open(sys.argv[4], "rb") as image_file_right:
        encoded_string = base64.b64encode(image_file_right.read())
        data_base64_right = encoded_string.decode()

    # asms = [
    #     {"id": meta_id, "img_left": "data:image/png;base64, "+data_base64_left, "img_right": "data:image/png;base64, "+data_base64_right}
    # ]

    asm = {"id": meta_id, "img_left": "data:image/png;base64, "+data_base64_left, "img_right": "data:image/png;base64, "+data_base64_right}

    path_template = Path(sys.argv[1])
    path_base = path_template.parent
    name_template = path_template.name
    env = Environment(loader=FileSystemLoader(path_base))
    # template = env.get_template("template_mqc.html")
    template = env.get_template(name_template)

    filename = path_base/"{}_jinja_out_mqc.html".format(meta_id)

    with open(filename, mode="w", encoding="utf-8") as output:
        output.write(template.render(asm=asm))


if __name__ == "__main__":
    sys.exit(main())