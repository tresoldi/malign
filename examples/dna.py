import malign


def main():
    # First test, pwdumb, with no scorer
    dna_seq1 = [c for c in "GATTACA"]
    dna_seq2 = [c for c in "ATAC"]
    alms = malign.multi_align([dna_seq1, dna_seq2], method="dumb")
    print("Experiment #1")
    print(malign.tabulate_alms(alms))

    # Second test, pwnw, default scorer
    alms = malign.multi_align([dna_seq1, dna_seq2], method="nw")
    print("Experiment #2")
    print(malign.tabulate_alms(alms))

    # Third test, pwkbest, default scorer
    alms = malign.multi_align([dna_seq1, dna_seq2], k=4, method="yenksp")
    print("Experiment #3")
    print(malign.tabulate_alms(alms))

    # Fourth test, pwkbest, DNA scorer
    scorer = malign.utils.DNA_MATRIX
    alms = malign.multi_align([dna_seq1, dna_seq2], k=4, method="yenksp", scorer=scorer)
    print("Experiment #4")
    print(malign.tabulate_alms(alms))

    # Fifth test, kbest, DNA scorer


#    scorer = malign.utils.DNA_MATRIX.copy()
#    graph = malign.yenksp.compute_graph(dna_seq1, dna_seq2, scorer)
#    dest = (len(dna_seq1), len(dna_seq2))
#    alms = malign.yenksp.align(graph, ((0, 0), dest), dna_seq1, dna_seq2, 4, gap_open=5)
#    print("Experiment #5")
#    print(malign.tabulate_alms(alms))

# Sixth test, kbest, bad DNA scorer
#    scorer[("T", "T")] = -50
#    graph = malign.yenksp.compute_graph(dna_seq1, dna_seq2, scorer)
#    dest = (len(dna_seq1), len(dna_seq2))
#    alms = malign.yenksp.align(graph, ((0, 0), dest), dna_seq1, dna_seq2, 4, gap_open=5)
#    print("Experiment #6")
#    print(malign.tabulate_alms(alms))


if __name__ == "__main__":
    main()
