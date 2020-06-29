# Import Python standard libraries
import csv
import os
from pathlib import Path

# Cache `resource` path for easy reusage
# TODO: move to main namespace?
RESOURCE_PATH = Path(os.path.realpath(__file__)).parent.parent / "resources"

# TODO: return dictionary with language names?
# TODO: remove spaces that are in data, as "m a i   d i r e  m a i", either
#       here in code (better) or in data
def read_wiktionary_data(languages):
    # Make sure at least two languages were requested
    if len(languages) < 2:
        raise ValueError("At least two languages must be requested.")

    # Collect data
    filename = RESOURCE_PATH / "wiktionary.tsv"
    with open(filename) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        data = [[row[lang] for lang in languages] for row in reader]

    # Remove entries without enough data
    data = [entry for entry in data if len([seq for seq in entry if seq]) >= 2]

    return data


if __name__ == "__main__":
    wikt = read_wiktionary_data(["Italian", "Russian"])
    for e in wikt:
        print(e)
