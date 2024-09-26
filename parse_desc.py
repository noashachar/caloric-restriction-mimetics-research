import re
from glob import glob

"""
In ".desc.txt" files each row is in a `key: value` format.
For example, a ".desc.txt" may look like:

Caloric Restriction Duration: 8 months
Deficit Size: 60% of daily intake
Organism: mouse
"""


def name_if_non_given(filename):
    _, category, number, *_ = re.split(r"[/.]", filename)
    return f"{category} Group {number}"


def parse_desc_file(filename: str):
    with open(filename) as file:
        lines = file.readlines()

    attributes = {}

    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue

        colon_idx = line.index(":")
        key = line[:colon_idx].strip()
        value = line[colon_idx + 1 :].strip()
        attributes[key] = value

    if "Name" not in attributes:
        attributes["Name"] = name_if_non_given(filename)

    return attributes


treatment_names = [
    parse_desc_file(filename)["Name"]
    for filename in sorted(glob("Research-Data/**/*.desc.txt"))
]
