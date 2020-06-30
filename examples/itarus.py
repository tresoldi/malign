# Import Python standard libraries
import csv
import os
from pathlib import Path

import catcoocc

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

        data = []
        for row in reader:
            sentinel = True
            for lang in languages:
                if len(row[lang].strip().split(" ")) <= 2:
                    sentinel = False

            if sentinel:
                data.append([row[lang].split(" ") for lang in languages])

#    # Remove entries without enough data
#    data = [entry for entry in data if len([seq for seq in entry if seq]) >= 2]

    return data

def main():
    wikt = read_wiktionary_data(["Italian", "Russian"])

#    for x in wikt:
#        print(x, "".join(x[0]), "".join(x[1]))
#    return

    cooccs = catcoocc.collect_cooccs(wikt)
    scorer = catcoocc.scorer.CatScorer(cooccs)
    matrix = scorer.tresoldi() # theil_u()
    matrix = catcoocc.scorer.scale_scorer(matrix, nrange=(0, 10))
    matrix = catcoocc.scorer.invert_scorer(matrix)

    for c in sorted(set(cooccs)):
        print(c, matrix[c])


if __name__ == "__main__":
    main()
