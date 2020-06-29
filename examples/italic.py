import os
import csv
from pathlib import Path
from collections import defaultdict

import catcoocc

RESOURCE_PATH = Path(os.path.realpath(__file__)).parent.parent / "resources"

def load_data(languages):
    # Collect data
    filename = RESOURCE_PATH / "northeuralex_italic.tsv"

    # Collect cogids
    cogid = defaultdict(dict)
    with open(filename) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for row in reader:
            cogid[row['COGID']][row['LANGUAGE']]= row['ALIGNMENT'].split(' ')

    # Only keep cogids with the languages we want
    filter = {}
    for identifier, values in cogid.items():
        found =[lang in values for lang in languages]
        if all(found):
            filter[identifier] = {lang:values[lang] for lang in languages}

    data = []
    for identifier, values in filter.items():
        data.append( [values[lang] for lang in languages])

    return data

def main():
    # Load data
    data = load_data(['Italian', 'Spanish'])

    # Compute cooccs
    cooccs = catcoocc.collect_cooccs(data)
    scorer = catcoocc.scorer.CatScorer(cooccs)

    s = scorer.theil_u()
    s = catcoocc.scorer.scale_scorer(s, nrange=(0, 10))
    s = catcoocc.scorer.invert_scorer(s)

    for c in sorted(set(cooccs)):
        print(c, s[c])

if __name__ == "__main__":
    main()
