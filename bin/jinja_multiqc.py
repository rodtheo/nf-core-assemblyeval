#!/usr/bin/env python
from jinja2 import Environment, FileSystemLoader
import sys
from pathlib import Path
import base64

def main():

    with open(sys.argv[2], "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
        data_base64 = encoded_string.decode()

    asms = [
        {"id": "Bsubtilis", "img": "data:image/png;base64, "+data_base64},
        {"id": "Aterreus", "img": "data:image/png;base64, "+data_base64}
    ]

    path_template = Path(sys.argv[1])
    path_base = path_template.parent
    name_template = path_template.name
    env = Environment(loader=FileSystemLoader(path_base))
    print("AAAAA", path_base)
    # template = env.get_template("template_mqc.html")
    template = env.get_template(name_template)

    filename = path_base/"jinja_out_mqc.html"

    with open(filename, mode="w", encoding="utf-8") as output:
        output.write(template.render(asms=asms))


if __name__ == "__main__":
    sys.exit(main())