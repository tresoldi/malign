import malign
from malign import tabulate_alms


# TODO: return an object and not a dictionary
# TODO: Scoring matrix should be callable/printable directly


def pairwise_dumb():
    # Pairwise dumb alignment
    seq_a = "tra"
    seq_b = "fatata"

    print("--------------------- PAIRWISE DUMB")
    alms = malign.multi_align([seq_a, seq_b], method="dumb")
    print(tabulate_alms(alms))


def multiwise_dumb():
    seqs = ["tra", "fra", "batata", "virp", "x"]

    print("--------------------- MULTIWISE DUMB")
    alms = malign.multi_align(seqs, method="dumb")
    print(tabulate_alms(alms))


def dna1():
    seq_a = "GATTACA"
    seq_b = "A"

    print("--------------------- DNA1")
    alms = malign.multi_align(
        [seq_a, seq_b], k=2, method="anw", matrix=malign.utils.DNA_MATRIX
    )
    print(tabulate_alms(alms))


def dna2():
    seq_a = "GATTACA"
    seq_b = "ATTT"

    print("--------------------- DNA2")
    alms = malign.multi_align(
        [seq_a, seq_b], k=2, method="anw", matrix=malign.utils.DNA_MATRIX
    )
    print(tabulate_alms(alms))


def modified_dna():
    modified_DNA = malign.utils.DNA_MATRIX.copy()
    modified_DNA["C", "T"] = -99.0

    seq_a = "GATTACA"
    seq_b = "ATTT"
    alms = malign.multi_align([seq_a, seq_b], k=4, method="anw", matrix=modified_DNA)

    print("--------------------- MODIFIED DNA #1")
    print(tabulate_alms(alms))

    print("--------------------- MODIFIED DNA #2")
    alms = malign.multi_align(
        ["GATTACA", "GATTATA"], k=2, method="anw", matrix=modified_DNA
    )
    print(tabulate_alms(alms))

    print("--------------------- MODIFIED DNA #3")
    alms = malign.multi_align(
        ["GATTATA", "GATTACA"], k=2, method="anw", matrix=modified_DNA
    )
    print(tabulate_alms(alms))


def with_ita_rus():
    ita_rus = malign.ScoringMatrix("docs/ita_rus.matrix")

    print("--------------------- ITA_RUS #1")
    alms = malign.multi_align(["atomo", "атом"], k=2, method="anw", matrix=ita_rus)
    print(tabulate_alms(alms))

    print("--------------------- ITA_RUS #2")
    alms = malign.multi_align(["Giacomo", "Яков"], k=4, method="anw", matrix=ita_rus)
    print(tabulate_alms(alms))


def with_ita_grk():
    ita_grk = malign.ScoringMatrix("docs/ita_grk.matrix")

    print("--------------------- ITA_GRK #1")
    alms = malign.multi_align(["atomo", "ατομο"], k=2, method="anw", matrix=ita_grk)
    print(tabulate_alms(alms))

    print("--------------------- ITA_GRK #2")
    alms = malign.multi_align(["Giacomo", "Ιακωβος"], k=4, method="anw", matrix=ita_grk)
    print(tabulate_alms(alms))


def with_full_matrix():
    ita_rus = malign.ScoringMatrix("docs/ita_rus.matrix")
    ita_grk = malign.ScoringMatrix("docs/ita_grk.matrix")

    # Combine the two matrices into a single one, add some points, show a couple of examples
    # TODO: move to function
    scores_ita_rus = {
        (key[0], key[1], None): value for key, value in ita_rus.scores.items()
    }
    scores_ita_grk = {
        (key[0], None, key[1]): value for key, value in ita_grk.scores.items()
    }
    scores = {**scores_ita_rus, **scores_ita_grk}

    full_matrix = malign.ScoringMatrix(scores)
    full_matrix["o", "в", "ο"] = -4.0
    full_matrix["o", "в", "ς"] = -10.0
    full_matrix["-", "в", "ς"] = -7.5
    full_matrix["-", "в", "ο"] = -6.5
    full_matrix["o", "-", "ο"] = 5.0
    full_matrix["o", "-", "ς"] = -3.0
    full_matrix["i", "-", "Ι"] = -4.0
    full_matrix["c", "к", "κ"] = 10.0
    full_matrix["-", "в", "ς"] = -10.0
    full_matrix["m", "в", "β"] = 10.0

    print("--------------------- FULL MATRIX #1")
    for key in [
        ("-", "к", "ο"),
        ("i", "а", "Ι"),
        ("m", "в", "β"),
        ("m", "-", "β"),
        ("-", "в", "ς"),
        ("-", "-", "ς"),
        ("o", "-", "ς"),
        ("-", "-", "ο"),
        ("-", "в", "ο"),
    ]:
        print(key, full_matrix[key])

    print("--------------------- FULL MATRIX #2")
    alms = malign.multi_align(
        ["atomo", "атом", "ατομο"], k=4, method="anw", matrix=full_matrix
    )
    print(tabulate_alms(alms))

    print("--------------------- FULL MATRIX #3")
    alms = malign.multi_align(
        ["Giacomo", "Яков", "Ιακωβος"], k=4, method="anw", matrix=full_matrix
    )
    print(tabulate_alms(alms))

    return full_matrix


def voldemort():
    seqs = ["VOLDEMORT", "WALDEMAR", "VLADIMIR", "VOLODYMIR"]
    voldemort_matrix = malign.utils.identity_matrix(seqs)

    print("--------------------- VOLDEMORT #1")
    for key in [
        ("T", "-", "-", "-"),
        ("T", "R", "R", "R"),
        ("-", "R", "R", "R"),
        ("R", "R", "R", "R"),
        ("O", "A", "I", "I"),
    ]:
        print(key, voldemort_matrix[key])

    print("--------------------- VOLDEMORT #2")
    alms = malign.multi_align(seqs, k=4, method="anw", matrix=voldemort_matrix)
    print(tabulate_alms(alms))

    return voldemort_matrix


def yenksp(full_matrix, voldemort_matrix):
    print("--------------------- YENKSP #1")
    alms = malign.multi_align(
        ["atomo", "атом", "ατομο"], k=2, method="yenksp", matrix=full_matrix
    )
    print(tabulate_alms(alms))

    print("--------------------- YENKSP #2")
    alms = malign.multi_align(
        ["Giacomo", "Яков", "Ιακωβος"], k=2, method="yenksp", matrix=full_matrix
    )
    print(tabulate_alms(alms))

    print("--------------------- YENKSP #3")
    alms = malign.multi_align(
        ["VOLDEMORT", "WALDEMAR", "VLADIMIR", "VOLODYMIR"],
        k=4,
        method="yenksp",
        matrix=voldemort_matrix,
    )
    print(tabulate_alms(alms))


def impute():
    s = {
        ("a", "A", "1"): -1,
        ("a", "A", "2"): 4,
        ("a", "A", "3"): 3,
        ("a", "B", "1"): 1,
        ("b", "A", "1"): -10,
        ("b", "A", "2"): 10,
        ("b", "A", "3"): 10,
        ("b", "A", "4"): 10,
        ("c", "B", "1"): 2,
        ("c", "B", "2"): 2,
        ("a", "-", "-"): -2,
        ("b", "-", "-"): -2,
        ("-", "A", "-"): -20,
        ("-", "B", "-"): -20,
        ("-", "-", "1"): -3,
        ("a", "B", "-"): -10,
        ("-", "A", "1"): -100,
        ("-", "A", "2"): -10,
        ("-", "B", "3"): -5,
    }
    m = malign.ScoringMatrix(s, impute_method="bayesian_ridge")

    print("--------------------- IMPUTE")
    for key in [
        ["a", "A", "2"],
        ["-", "-", "3"],
        ["c", "-", "4"],
        ["a", "A", None],
        ["a", None, "3"],
        [None, "B", "1"],
    ]:
        print(key, m[key])


def impute_ita_rus_grk():
    ita_rus = {
        ("a", "а"): 10,
        ("a", "т"): -4,
        ("a", "о"): 3,
        ("a", "м"): -3,
        ("a", "Я"): 5,
        ("a", "к"): -4,
        ("a", "в"): -4,
        ("a", "-"): -1,
        ("t", "а"): -4,
        ("t", "т"): 10,
        ("t", "о"): -4,
        ("t", "м"): -3,
        ("t", "Я"): -4,
        ("t", "к"): -1,
        ("t", "в"): -3,
        ("t", "-"): -3,
        ("o", "а"): 4,
        ("o", "т"): -4,
        ("o", "о"): 10,
        ("o", "м"): -5,
        ("o", "Я"): 2,
        ("o", "к"): -4,
        ("o", "в"): -2,
        ("o", "-"): 0,
        ("m", "а"): -2,
        ("m", "т"): -1,
        ("m", "о"): -4,
        ("m", "м"): 9,
        ("m", "Я"): -3,
        ("m", "к"): -2,
        ("m", "в"): 2,
        ("m", "-"): 1,
        ("G", "а"): -3,
        ("G", "т"): 1,
        ("G", "о"): -3,
        ("G", "м"): -1,
        ("G", "Я"): 4,
        ("G", "к"): 3,
        ("G", "в"): -2,
        ("G", "-"): 2,
        ("i", "а"): 5,
        ("i", "т"): -3,
        ("i", "о"): 4,
        ("i", "м"): 0,
        ("i", "Я"): 6,
        ("i", "к"): -3,
        ("i", "в"): -3,
        ("i", "-"): 2,
        ("c", "а"): -4,
        ("c", "т"): -1,
        ("c", "о"): -5,
        ("c", "м"): -4,
        ("c", "Я"): -5,
        ("c", "к"): 4,
        ("c", "в"): 1,
        ("c", "-"): -3,
        ("-", "а"): 1,
        ("-", "т"): -2,
        ("-", "о"): -1,
        ("-", "м"): 2,
        ("-", "Я"): 2,
        ("-", "к"): -4,
        ("-", "в"): 2,
        ("-", "-"): 0,
    }

    ita_grk = {
        ("a", "α"): 10,
        ("a", "-"): -5,
        ("-", "α"): -5,
        ("o", "α"): 4,
        ("i", "α"): 5,
        ("t", "τ"): 10,
        ("c", "τ"): 2,
        ("a", "ο"): 6,
        ("o", "ο"): 10,
        ("o", "-"): -10,
        ("-", "ο"): 10,
        ("m", "μ"): 10,
        ("a", "Ι"): 2,
        ("i", "Ι"): 7,
        ("t", "κ"): 2,
        ("c", "κ"): 8,
        ("a", "ω"): 4,
        ("o", "ω"): 10,
        ("m", "β"): 4,
        ("v", "β"): 5,
        ("s", "ς"): 6,
    }

    ita_rus_m = malign.ScoringMatrix(ita_rus, impute_method="bayesian_ridge")
    ita_grk_m = malign.ScoringMatrix(ita_grk, impute_method="bayesian_ridge")

    print("--------------------- IMPUTE ITA RUS GRK")
    print("""ita_rus_m["c", "Я"]""", ita_rus_m["c", "Я"])  # provided
    print("""ita_grk_m["t", "μ"]""", ita_grk_m["t", "μ"])  # inferred

    # we provide a single point ita/rus/grk
    # TODO: use builder
    scores_ita_rus = {
        (key[0], key[1], None): value for key, value in ita_rus_m.scores.items()
    }
    scores_ita_grk = {
        (key[0], None, key[1]): value for key, value in ita_grk_m.scores.items()
    }
    scores = {**scores_ita_rus, **scores_ita_grk}
    scores["t", "т", "τ"] = 10

    irg_m = malign.ScoringMatrix(scores, impute_method="bayesian_ridge")

    for key in [
        ["c", "Я", None],
        ["t", None, "μ"],
        [None, "Я", "μ"],
        ["a", "а", "α"],
        ["a", "а", "-"],
        ["a", "Я", "μ"],
    ]:
        print(f"irg_m[{key}]", irg_m[key])

    irg_m.save("ita_rus_grk.matrix")

    m2 = malign.ScoringMatrix("ita_rus_grk.matrix")

    for key in [
        ["c", "Я", None],
        ["t", None, "μ"],
        [None, "Я", "μ"],
        ["a", "а", "α"],
        ["a", "а", "-"],
        ["a", "Я", "μ"],
        ["a", "Я", None],
    ]:
        print(f"m2[{key}]", m2[key])


def scoring_matrix():
    # Define a scoring matrix for a simple, two domain system
    # alphabet a is {'-', 'a', 'b', 'c'}, alphabet b is {'-', 'X', 'Y'}
    vectors = {
        ("-", "-"): 0.0,
        ("-", "X"): -3.0,
        ("-", "Y"): -9.0,
        ("a", "-"): -3.0,
        ("a", "X"): 0.0,
        ("a", "Y"): 8.0,
        ("b", "-"): -5.0,
        ("b", "X"): 4.0,
        ("b", "Y"): 4.0,
        ("c", "-"): 2.0,
        ("c", "X"): -1.0,
        ("c", "Y"): 7.0,
    }

    print("--------------------- SCORING MATRIX #1")
    matrix = malign.ScoringMatrix(vectors, impute_method="bayesian_ridge")
    print(matrix.tabulate())

    v = vectors.copy()
    v.pop(("a", "-"))
    v.pop(("c", "Y"))
    matrix = malign.ScoringMatrix(v, impute_method="bayesian_ridge")
    print("--------------------- SCORING MATRIX #2")
    print(matrix.tabulate())

    # Define a scoring matrix for a simple, three domain system
    # alphabet a is {'-', 'a', 'b', 'c'}, alphabet b is {'-', 'X', 'Y'}, alphabet c is {'-', 'i', 'j'}
    # note that we already comment out some entries to make gaps
    vectors = {
        ("-", "-", "-"): 0.0,
        #    ('-', '-', 'i') : -4.0,
        ("-", "-", "j"): -8.0,
        ("-", "X", "-"): -5.0,
        ("-", "Y", "-"): -5.0,
        ("a", "-", "-"): 0.0,
        #    ('a', 'X', '-') :  0.0,
        ("a", "Y", "-"): 8.0,
        ("b", "-", "-"): 0.0,
        ("b", "X", "-"): 4.0,
        ("b", "Y", "-"): 4.0,
        ("c", "-", "-"): 0.0,
        ("c", "X", "-"): -1.0,
        ("c", "Y", "-"): -7.0,
        ("-", "X", "i"): -3.0,
        ("-", "Y", "i"): -9.0,
        ("a", "-", "i"): -3.0,
        #    ('a', 'X', 'i') :  0.0,
        #    ('a', 'Y', 'i') :  8.0,
        ("b", "-", "i"): -5.0,
        ("b", "X", "i"): 4.0,
        ("b", "Y", "i"): 4.0,
        ("c", "-", "i"): 2.0,
        ("c", "X", "i"): -1.0,
        ("c", "Y", "i"): 7.0,
        ("-", "X", "j"): -5.0,
        ("-", "Y", "j"): -6.0,
        ("a", "-", "j"): 3.0,
        ("a", "X", "j"): 8.0,
        ("a", "Y", "j"): 8.0,
        ("b", "-", "j"): -6.0,
        ("b", "X", "j"): 5.0,
        ("b", "Y", "j"): 4.0,
        ("c", "-", "j"): -5.0,
        #    ('c', 'X', 'j') :  6.0,
        ("c", "Y", "j"): 7.0,
    }

    # TODO: have no fill method

    print("Standard fill method")
    matrix = malign.ScoringMatrix(vectors)
    print("--------------------- SCORING MATRIX #3")
    print(matrix.tabulate())

    print("Distance fill method")
    matrix = malign.ScoringMatrix(vectors, impute_method="bayesian_ridge")
    print("--------------------- SCORING MATRIX #4")
    print(matrix.tabulate())

    vector_01 = {
        ("-", "X"): -3.0,
        ("a", "-"): -3.0,
        ("a", "X"): 0.0,
        ("a", "Y"): 8.0,
        ("b", "-"): -5.0,
        ("b", "Y"): 4.0,
        ("c", "X"): -1.0,
        ("c", "Y"): 7.0,
    }
    vector_02 = {
        ("a", "-"): -4.0,
        ("a", "i"): 2.0,
        ("a", "j"): 2.0,
        ("b", "i"): -5.0,
        ("b", "j"): 9.0,
        ("c", "-"): -7.0,
        ("c", "j"): 4.0,
    }

    print("--------------------- SCORING MATRIX #5 (0, 1)")
    matrix01 = malign.ScoringMatrix(vector_01)
    print(matrix01.tabulate())

    print("--------------------- SCORING MATRIX #5 (0, 2)")
    matrix02 = malign.ScoringMatrix(vector_02)
    print(matrix02.tabulate())

    print("--------------------- SCORING MATRIX #5 (0, 1, 22)")
    # TODO: use compounder
    scores_01 = {
        (key[0], key[1], None): value for key, value in matrix01.scores.items()
    }
    scores_02 = {
        (key[0], None, key[1]): value for key, value in matrix02.scores.items()
    }
    scores = {**scores_01, **scores_02}

    matrix = malign.ScoringMatrix(scores, impute_method="bayesian_ridge")
    print(matrix.tabulate())

    print("--------------------- SCORING MATRIX #6 (IDENTITY)")
    seqs = ["ab", "aab", "bbb"]
    id_matrix = malign.utils.identity_matrix(seqs, match=2, gap_score=-3)
    print(id_matrix.tabulate())


############
pairwise_dumb()
multiwise_dumb()
dna1()
dna2()
modified_dna()
with_ita_rus()
with_ita_grk()
full_matrix = with_full_matrix()
voldemort_matrix = voldemort()
yenksp(full_matrix, voldemort_matrix)
impute()
impute_ita_rus_grk()
scoring_matrix()
