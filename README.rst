MAlign
======

|PyPI| |Build Status| |codecov| |Codacy Badge|

MALIGN is a library for performing multiple alignments on sequences of
different alphabets. It allows each sequence to have its own domain,
which in turns allows to use asymmetric and sparse scoring matrices,
including on gaps, and to perform real, single-pass multiple alignment,
allowing to compute ``k``-best alignments. While intended for linguistic
usage mostly, it can be used for aligning any type of sequential
representation, and it is particularly suitable as a general-purpose
tool for cases where there are no prior hypotheses on the scoring
matrices.

Installation and usage
----------------------

The library can be installed as any standard Python library with
``pip``, and used as demonstrated in the following snippet:

In any standard Python environment, ``malign`` can be installed with:

.. code:: bash

   $ pip install malign

For most purposes, it is enough to pass the sequences to be aligned and
a method (such as ``anw`` or ``yenksp``) to the ``.multi_align()``
function:

.. code:: python

   >>> import malign                                                                                                      
   >>> alms = malign.multi_align(["ATTCGGAT", "TACGGATTT"], "anw", k=2)                                                   
   >>> print(malign.tabulate_alms(alms))                                                                                  
   | Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |  #4  |  #5  |  #6  |  #7  |  #8  |  #9  |
   |-------|-------|---------|------|------|------|------|------|------|------|------|------|------|
   | 0     | A     |   -0.29 |  A   |  T   |  T   |  C   |  G   |  G   |  A   |  -   |  T   |  -   |
   | 0     | B     |   -0.29 |  -   |  T   |  A   |  C   |  G   |  G   |  A   |  T   |  T   |  T   |
   |       |       |         |      |      |      |      |      |      |      |      |      |      |
   | 1     | A     |   -0.29 |  A   |  T   |  T   |  C   |  G   |  G   |  A   |  -   |  -   |  T   |
   | 1     | B     |   -0.29 |  -   |  T   |  A   |  C   |  G   |  G   |  A   |  T   |  T   |  T   |

Scoring matrices can be either computed with the auxiliary methods,
including various optimizations, or read from JSON files:

.. code:: python

   >>> ita_rus = malign.ScoringMatrix(filename="docs/ita_rus.matrix")
   >>> alms = malign.multi_align(["Giacomo", "Яков"], k=4, method="anw", matrix=ita_rus)
   >>> print(malign.tabulate_alms(alms))
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

The library can also be used by means of the command-line ``malign``
tool. If no matrix is provided, an identity one is used by default.

.. code:: bash

   $ ▶ malign baba,maa
   | Idx   | Seq   |   Score |  #0  |  #1  |  #2  |  #3  |
   |-------|-------|---------|------|------|------|------|
   | 0     | A     |   -0.47 |  b   |  a   |  b   |  a   |
   | 0     | B     |   -0.47 |  m   |  a   |  -   |  a   |

   $ ▶ malign --matrix docs/ita_rus.matrix -k 6 Giacomo,Яков
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
   |       |       |         |      |      |      |      |      |      |      |      |
   | 4     | A     |    2.12 |  G   |  i   |  a   |  c   |  o   |  m   |  -   |  o   |
   | 4     | B     |    2.12 |  -   |  Я   |  -   |  к   |  о   |  -   |  в   |  -   |
   |       |       |         |      |      |      |      |      |      |      |      |
   | 5     | A     |    2.12 |  G   |  i   |  a   |  c   |  o   |  -   |  m   |  o   |
   | 5     | B     |    2.12 |  -   |  Я   |  -   |  к   |  о   |  в   |  -   |  -   |

Changelog
---------

Version 0.1: - First release for internal announcement, testing, and
community outreach

Version 0.2: - Major revision with asymmetric Needleman-Wunsch and Yen’s
``k``-shortest path implementation. - Added scoring matrix object - Sort
alignments in consistent and reproducible ways, even when the alignment
score is the same

Roadmap
-------

Version 0.3: - Complete documentation and setup ``readthedocs`` - Add
new inference method to sparse matrices using impurity/entropy -
Describe matrix filling methods in more detail - Consider implementation
of UPGMA and NJ multiple alignment - Add function/method to visualize
the graphs used for the ``yenksp`` methods - Implement blocks and local
search in ``anw`` and ``yenksp``, with different starting/ending
positions - Implement memoization where possible - Consider expanding
dumb_malign by adding random gaps (``pad_align``), as an additional
baseline method - Allow ``anw`` to work within a threshold percentage of
the best score - Implement a method combining the results of the
different algorithms - Add methods and demonstration for matrix
optimization

Community guidelines
--------------------

While the author can be contacted directly for support, it is
recommended that third parties use GitHub standard features, such as
issues and pull requests, to contribute, report problems, or seek
support.

Contributing guidelines, including a code of conduct, can be found in
the ``CONTRIBUTING.md`` file.

Author and citation
-------------------

The library is developed by Tiago Tresoldi (tresoldi@shh.mpg.de).

The author has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation
programme (grant agreement No. \ `ERC Grant
#715618 <https://cordis.europa.eu/project/rcn/206320/factsheet/en>`__,
`Computer-Assisted Language Comparison <https://digling.org/calc/>`__.

If you use ``malign``, please cite it as:

   Tresoldi, Tiago (2020). MALIGN, a library for multiple asymmetric
   alignments on different alphabets. Version 0.2. Jena.

In BibTeX:

.. code:: bibtex

   @misc{Tresoldi2020malign,
     author = {Tresoldi, Tiago},
     title = {MALIGN, a library for multiple asymmetric alignments on different alphabets. Version 0.2},
     howpublished = {\url{https://github.com/tresoldi/malign}},
     address = {Jena},
     publisher = {Max Planck Institute for the Science of Human History}
     year = {2020},
   }

.. |PyPI| image:: https://img.shields.io/pypi/v/malign.svg
   :target: https://pypi.org/project/malign
.. |Build Status| image:: https://travis-ci.org/tresoldi/malign.svg?branch=master
   :target: https://travis-ci.org/tresoldi/malign
.. |codecov| image:: https://codecov.io/gh/tresoldi/malign/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/tresoldi/malign
.. |Codacy Badge| image:: https://api.codacy.com/project/badge/Grade/f6428290a03742e69a6a5cb512a99650
   :target: https://www.codacy.com/manual/tresoldi/malign?utm_source=github.com&utm_medium=referral&utm_content=tresoldi/malign&utm_campaign=Badge_Grade
