# MAlign

[![PyPI](https://img.shields.io/pypi/v/malign.svg)](https://pypi.org/project/malign)
[![Build Status](https://travis-ci.org/tresoldi/malign.svg?branch=master)](https://travis-ci.org/tresoldi/malign)
[![codecov](https://codecov.io/gh/tresoldi/malign/branch/master/graph/badge.svg)](https://codecov.io/gh/tresoldi/malign)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f6428290a03742e69a6a5cb512a99650)](https://www.codacy.com/manual/tresoldi/malign?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tresoldi/malign&amp;utm_campaign=Badge_Grade)

MALIGN is a library for multiple asymmetric alignments on different alphabets.
It is currently under initial research and development, but can already be
used to obtain multiple alignments for DNA sequences.


## Installation and usage

The library can be installed as any standard Python library with
`pip`, and used as demonstrated in the following snippet:

In any standard Python environment, `malign` can be installed with:

```bash
$ pip install malign
```

For most purposes, it is enough to pass the two sequences to be aligned,
along with a scorer, to the `get_aligns()` function:


```python
>> import malign
>> seq1 = ["A", "T", "T", "C", "G", "G", "A", "T"]
>> seq2 = ["T", "A", "C", "G", "G", "A", "T", "T", "T"]
>> graph = malign.compute_graph(seq1, seq2, malign.DNA_SCORER)
>> dest = "%i:%i" % (len(seq1), len(seq2))
>> aligns = malign.get_aligns(graph, ("0:0", dest), seq1, seq2, 4)
>> for idx, align in enumerate(aligns):
>>   print(" ".join(align[0][0]), align[1])
>>   print(" ".join(align[0][1]), "\n")
A T T C G G A - - T 84.0
T A - C G G A T T T

A T T C G G A - - T 84.0
T - A C G G A T T T

A T T C G G A T - - 84.0
T A - C G G A T T T

A T T C G G A T - - 84.0
T - A C G G A T T T
```

The library can also be used by means of the command-line `malign` tool:

```bash
$ malign --dna ATTCGGAT TACGGATTT
Alignment #0 (score: 84.00)
A T T C G G A - - T
T A - C G G A T T T

Alignment #1 (score: 84.00)
A T T C G G A - - T
T - A C G G A T T T

Alignment #2 (score: 84.00)
A T T C G G A T - -
T A - C G G A T T T

Alignment #3 (score: 84.00)
A T T C G G A T - -
T - A C G G A T T T
```

## Changelog

Version 0.1:
  - First release for internal announcement, testing, and community outreach

## Roadmap

Version 0.2:
  - Setup readthedocs
  - Sort in consistent and reproducible way all alignments, even when the
    score is the same
  - Deal with conflicting package versions due to `lingpy`, or write new
    NW implementation
  - Implement single-pass function with defaults, with scorer, graph,
    destnation, etc.

## Community guidelines

While the author can be contacted directly for support, it is recommended
that third parties use GitHub standard features, such as issues and
pull requests, to contribute, report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in
the `CONTRIBUTING.md` file.

## Author and citation

The library is developed by Tiago Tresoldi (tresoldi@shh.mpg.de).

The author has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation
programme (grant agreement
No. [ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en),
[Computer-Assisted Language Comparison](https://digling.org/calc/).

If you use `malign`, please cite it as:

  > Tresoldi, Tiago (2020). MALIGN, a library for multiple asymmetric alignments on different alphabets. Version 1.0. Jena. DOI: 10.5281/zenodo.3668870

  In BibTeX:

```bibtex
@misc{Tresoldi2020malign,
  author = {Tresoldi, Tiago},
  title = {MALIGN, a library for multiple asymmetric alignments on different alphabets. Version 0.1.},
  howpublished = {\url{https://github.com/tresoldi/malign}},
  address = {Jena},
  year = {2020},
}
```
