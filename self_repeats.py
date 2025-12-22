import pandas as pd
import numpy as np
import sys


def load_test_seq(fn: str, n: int):
    df = pd.read_csv(fn)
    return df['text'][n-1]


def duplicate(motif: str, target: str):
    n = len(target)
    m = len(motif)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    dp[0][0] = 0
    for i in range(1, n + 1):
        dp[i][0] = -i
    for j in range(1, m + 1):
        dp[0][j] = -j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            v = dp[i-1][j] - 1
            h = dp[i][j-1] - 1
            d = dp[i-1][j-1] + (1 if target[i-1] == motif[j-1] else -1)
            dp[i][j] = max(v, h, d)
        dp[i][0] = dp[i][m]
        for j in range(1, m + 1):
            dp[i][j] = max(dp[i][j], dp[i][j-1] - 1)
        # FIXME: Is it correct?
        # if dp[i].max() < 0:
        #     break
    max_score, ext = -n * m, -1
    for i in range(1, n + 1):
        # FIXME: I ignore the incomplete copy that might be better for easy implementation
        if dp[i][m] > max_score:
            max_score = dp[i][m]
            ext = i
    if max_score < 0:
        max_score = 0
        ext = 0
    return max_score, ext


def solve(text: str):
    # Let F(i) denote the DP score of S[0,i] and is calculated as:
    # F(i) = F(i-1) + 1
    #   diagonal transfer without repeats.
    # F(i) = F(k) + dup(S[i-l+1, i], S[k+1, i] + G
    #   a repeat motif length of l copies a certain number of times to be aligned with S[k+1, i],
    #   where 1 < l <= L (maximum repeat length allowed) and 0 <= k <= i - 2 * l + 1
    #   G is the initial penalty of a repeat event.

    # To drive the optimal path away from diagonal, a repeat is given a bonus score B each time of copy;
    # hence, repeats with shorter motif length and greater copy times are preferred in all possible representations.
    # Currently, I use a simple scoring matrix that doesn't support affine-gap penalty for dup function.
    # match = 1, mismatch = -1, gap (insertion / deletion) = -1
    # B = 1 and G = -5

    min_repeat_len = 10
    open_penalty = -5
    max_motif_len = 100
    n = len(text)
    dp = np.zeros(n, dtype=int)
    for i in range(0, min_repeat_len):
        dp[i] = i + 1
    for i in range(min_repeat_len, n):
        dp[i] = dp[i-1] + 1
        for l in range(2, min(max_motif_len, i // 2)):
            motif = text[i - l + 1: i + 1]
            target = text[0: i - l + 1]
            motif = motif[::-1]
            target = target[::-1]
            # print('motif=%s [%d,%d), target=%s [%d,%d)' % (motif, i-l+1, i+1, target, 0, i-l+1))
            w_score, ext = duplicate(motif, target)
            # print('w_score=%d, ext=%d, start=%d' % (w_score, ext, i - l - ext))
            # No reward for itself (it is not a copy)
            w = w_score + l + open_penalty
            if i - l - ext == -1:  # Extend to the start
                dp[i] = max(dp[i], w)
            else:
                dp[i] = max(dp[i], dp[i - l - ext] + w)
    n_repeats = 0
    for i in range(2, n):
        if dp[i] != dp[i-1] + 1:
            n_repeats += 1
    print('%d repeats identified' % (n_repeats))


if __name__ == '__main__':
    seq = load_test_seq('./motif.csv', int(sys.argv[1]))
    print('sequence length: %d' % (len(seq)))
    solve(seq)
