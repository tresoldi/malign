import numpy as np

# TODO: expand dumb_malign by adding random gaps, call this pad_align? maybe with swaps?
# TODO: receive scorer
def dumb_malign(seqs, gap="-", **kwargs):
    # Obtain the longest sequence length
    max_length = max([len(seq) for seq in seqs])

    # Pad all sequences in `alm`
    alm = {"seqs": []}
    scores = []
    for seq in seqs:
        # Computer lengths and bads
        num_pad = max_length - len(seq)
        left_pad_len = int(num_pad / 2)
        right_pad_len = num_pad - left_pad_len
        left_pad = [gap] * left_pad_len
        right_pad = [gap] * right_pad_len

        # Append the padded sequence and the score, here computed from the
        # number of gaps
        alm["seqs"].append([*left_pad, *list(seq), *right_pad])
        scores.append(1.0 - (num_pad / max_length))

    # Add overall score
    alm["score"] = np.mean(scores)

    # The `dumb` method will always return a single aligment, but we
    # still return a list for compatibility with other methods
    return [alm]
