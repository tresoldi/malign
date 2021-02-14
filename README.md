# MAlign

[![PyPI](https://img.shields.io/pypi/v/malign.svg)](https://pypi.org/project/malign)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f6428290a03742e69a6a5cb512a99650)](https://www.codacy.com/manual/tresoldi/malign?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tresoldi/malign&amp;utm_campaign=Badge_Grade)

MALIGN is a library for performing multiple alignments on sequences of different
domains, allowing the usage of asymmetric scoring matrices. Multiple alignments
are actual multiple alignments, scoring according to the overall probability
of each alignment site, and not a succession of pairwise alignments gradually
combined.

While intended for linguistic usage mostly, it can be used for aligning any type of
sequential representation as long as the elements of each domain are hashable. It is
particularly suitable as a general-purpose tool for cases where there are no prior
hypotheses on the scoring matrices, which can be inferred or imputed (including
from incomplete data), or optimized from observable examples to find local and
global minima that can be used to explain the relationships between the sequences.

## Installation and usage

The library can be installed as any standard Python library with
`pip`, preferably within a virtual environment:

```bash
$ pip install malign
```

For most purposes, it is enough to pass the sequences to be aligned and
specify one of the available methods (currently `anw`, the default, 
and `yenksp`) to the `.multi_align()` function, along with the maximum
number of alignments to be returned (`k`):

```python
>> import malign                                                                                                      
>> alms = malign.multi_align(["ATTCGGAT", "TACGGATTT"], k=2)                                                   
>> print(malign.tabulate_alms(alms))                                                                                  
| Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |  #4  |  #5  |  #6  |  #7  |  #8  |  #9  |
|-------|-------|---------|------|------|------|------|------|------|------|------|------|------|
| 0     | A     |   -0.29 |  A   |  T   |  T   |  C   |  G   |  G   |  A   |  -   |  T   |  -   |
| 0     | B     |   -0.29 |  -   |  T   |  A   |  C   |  G   |  G   |  A   |  T   |  T   |  T   |
|       |       |         |      |      |      |      |      |      |      |      |      |      |
| 1     | A     |   -0.29 |  A   |  T   |  T   |  C   |  G   |  G   |  A   |  -   |  -   |  T   |
| 1     | B     |   -0.29 |  -   |  T   |  A   |  C   |  G   |  G   |  A   |  T   |  T   |  T   |
```

Scoring matrices can be either computed with the auxiliary methods, including various
optimizations, or read from JSON files:

```python
>> ita_rus = malign.ScoringMatrix()
>> ita_rus.load("docs/ita_rus.matrix")
>> alms = malign.multi_align(["Giacomo", "Яков"], k=4, method="anw", matrix=ita_rus)
>> print(malign.tabulate_alms(alms))
| Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |  #4  |  #5  |  #6  |  #7  |
|-------|-------|---------|------|------|------|------|------|------|------|------|
| 0     | A     |    2.86 |  G   |  i   |  a   |  c   |  o   |  m   |  o   |      |
| 0     | B     |    2.86 |  -   |  Я   |  -   |  к   |  о   |  в   |  -   |      |
|       |       |         |      |      |      |      |      |      |      |      |
| 1     | A     |    2.29 |  G   |  i   |  a   |  c   |  o   |  m   |  o   |      |
| 1     | B     |    2.29 |  -   |  Я   |  -   |  к   |  о   |  -   |  в   |      |
|       |       |         |      |      |      |      |      |      |      |      |
| 2     | A     |    2.12 |  G   |  i   |  a   |  c   |  o   |  m   |  o   |  -   |
| 2     | B     |    2.12 |  -   |  Я   |  -   |  к   |  о   |  -   |  -   |  в   |
|       |       |         |      |      |      |      |      |      |      |      |
| 3     | A     |    2.12 |  G   |  i   |  a   |  c   |  o   |  m   |  o   |  -   |
| 3     | B     |    2.12 |  -   |  Я   |  -   |  к   |  -   |  -   |  о   |  в   |
```

More complex examples, including for matrix imputation and optimization, can
be found in the documentation.

## Changelog

Version 0.1:
  - First release for an internal announcement, testing, and community outreach

Version 0.2:
  - Major revision with asymmetric Needleman-Wunsch and Yen's `k`-shortest path
    implementation
  - Added scoring matrix object
  - Sort alignments in consistent and reproducible ways, even when the alignment
    score is the same

Version 0.3
  - Code improvements, including type annotation, and some refactoring
  - Allowing usage with any hashable Python object (not only strings)
  - Add methods for matrix imputation  
  - Update of documentation
  - General preparations for public announcement

## TODO

Version 0.3:
  - Complete documentation and setup `readthedocs`
  - Consider implementation of UPGMA and NJ multiple alignment
  - Add function/method to visualize the graphs used for the `yenksp` methods
  - Implement blocks and local search in `anw` and `yenksp`, with different
    starting/ending positions
  - Implement memoization where possible
  - Consider expanding dumb_malign by adding random gaps (`pad_align`), as an additional
    baseline method
  - Allow `anw` to work within a threshold percentage of the best score
  - Implement a method combining the results of the different algorithms
  - Add methods and demonstration for matrix optimization

## Community guidelines

While the author can be contacted directly for support, it is recommended
that third parties use GitHub standard features, such as issues and
pull requests, to contribute, report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in
the `CONTRIBUTING.md` file.

## Author and citation

The library is developed by Tiago Tresoldi (tiago.tresoldi@lingfil.uu.se).

During the first stages of development, the author received funding from the
European Research Council (ERC) under the European Union’s Horizon 2020
research and innovation programme (grant agreement
No. [ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en),
[Computer-Assisted Language Comparison](https://digling.org/calc/).

If you use `malign`, please cite it as:

  > Tresoldi, Tiago (2021). MALIGN, a library for multiple asymmetric alignments on
  > different domains. Version 0.3. Uppsala: Uppsala Universitet.

In BibTeX:

```bibtex
@misc{Tresoldi2021malign,
  author = {Tresoldi, Tiago},
  title = {MALIGN, a library for multiple asymmetric alignments on different domains. Version 0.3},
  howpublished = {\url{https://github.com/tresoldi/malign}},
  address = {Uppsala},
  publisher = {Uppsala Universitet}
  year = {2021},
}
```
