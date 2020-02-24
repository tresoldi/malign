# TODO: rename for pairwise and multiple
def align(seq_a, seq_b):
    # Temporary code for scaffolding
    if len(seq_a) < len(seq_b):
        seq_a, seq_b = seq_b, seq_a

    # Number of padding characters
    num_pad = len(seq_a) - len(seq_b)

    # One character for each, plus padding
    alm_a = [token for token in seq_a]
    alm_b = [token for token in seq_b]
    alm_pad = ["-"] * num_pad
    alm_b += alm_pad

    return alm_a, alm_b
