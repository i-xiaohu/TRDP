import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt


def load_test_seq(fn: str, n: int):
    df = pd.read_csv(fn)
    flank_l = int(df['flank_l'][n-1])
    motif = df['motif'][n-1]
    m = len(motif)
    seq = df['text'][n-1]
    text_len = len(seq)
    periods = int(df['periods'][n-1])
    # FIXME: it is inaccurate; use points given by WDP
    k = flank_l
    points = [-1] * text_len
    for i in range(1, periods):
        k += m
        points[k] = k - m
    return seq, points


def gary_benson(que: str, points_q: list, ref: str, points_r: list):
    n, m = len(que), len(ref)
    mat_s, mis_p, gap_o, gap_e = 1, -4, -6, -1
    E = np.zeros((n + 1, m + 1), dtype=int)
    H = np.zeros((n + 1, m + 1), dtype=int)
    bt = np.zeros((n + 1, m + 1), dtype=int)
    cp_pool = [0, 0, 0]
    DIAGNOAL, VERTICAL, HORIZONTAL, COPY = 0, 1, 2, 3
    H[0][0] = 0
    for j in range(1, m + 1):
        E[0][j] = H[0][j] = gap_o + gap_e * j
    for i in range(1, n + 1):
        H[i][0] = F = gap_o + gap_e * i
        bt[i][0] = VERTICAL
        for j in range(1, m + 1):
            F = max(H[i][j-1] + gap_o, F) + gap_e
            E[i][j] = max(H[i-1][j] + gap_o, E[i-1][j]) + gap_e
            M = H[i-1][j-1] + (mat_s if que[i-1] == ref[j-1] else mis_p)
            # Copy que motif multiple times to align with ref suffix
            C1, motif_len1, ext1 = min(F, E[i][j], M) - 100, 0, -1
            smodel_score, copy_score = -1, -1
            if points_q[i-1] != -1:
                motif = que[points_q[i-1]: i-1]
                motif_len1 = len(motif)
                suffix = ref[0: j]
                # Backward DP extension
                motif = motif[::-1]
                suffix = suffix[::-1]
                n2, m2 = len(suffix), len(motif)
                dp = np.zeros((n2 + 1, m2 + 1), dtype=int)
                # I use linear gap penalty for inner DP
                mat_s2, mis_p2, gap_p2, open_p2, copy_p2 = 1, -4, -1, -10, -3
                dp[0][0] = 0
                for j2 in range(1, m2 + 1):
                    dp[0][j2] = gap_p2 * j2
                for i2 in range(1, n2 + 1):
                    dp[i2][0] = gap_p2 * i2
                for i2 in range(1, n2 + 1):
                    for j2 in range(1, m2 + 1):
                        h = dp[i2][j2-1] + gap_p2
                        v = dp[i2-1][j2] + gap_p2
                        d = mat_s2 if suffix[i2-1] == motif[j2-1] else mis_p2
                        if j2 == 1:
                            h += copy_p2
                            d += copy_p2
                        dp[i2][j2] = max(h, v, d)
                    dp[i2][0] = dp[i2][m2]
                    for j2 in range(1, m2 + 1):
                        h = dp[i2][j2-1] + gap_p2
                        if j2 == 1:
                            h += copy_p2
                        dp[i2][j2] = max(dp[i2][j2], h)
                    if dp[i2].max() < 0:
                        n2 = i2
                        break
                for e in range(1, n2 + 1):
                    c = H[i][j-e] + open_p2 + dp[e][m2] - copy_p2
                    if c > C1:
                        smodel_score = H[i][j-e]
                        copy_score = dp[e][m2]
                        C1 = c
                        ext1 = e
            H[i][j] = max(M, E[i][j], F, C1)
            if H[i][j] == F:
                bt[i][j] = HORIZONTAL
            elif H[i][j] == M:
                bt[i][j] = DIAGNOAL
            elif H[i][j] == E[i][j]:
                bt[i][j] = VERTICAL
            elif H[i][j] == C1:
                # print('i=%d, j=%d, motif=%s, copy_score=%d, ext=%d' % (i, j, motif, C1, ext1))
                bt[i][j] = len(cp_pool)
                cp_pool.append((C1, smodel_score, copy_score, ext1))

    print('Score: %d' % H[n][m])
    n_del, n_ins, n_copy = 0, 0, 0
    i, j = n, m
    while i >= 0 and j >= 0:
        if bt[i][j] == HORIZONTAL:
            n_del += 1
            j -= 1
        elif bt[i][j] == VERTICAL:
            n_ins += 1
            i -= 1
        elif bt[i][j] == DIAGNOAL:
            i, j = i - 1, j - 1
        else:
            n_copy += 1
            k = bt[i][j]
            c1, ss, cs, ext = cp_pool[k]
            print('C=%d, s_model_score=%d, copy_score=%d, copy length=%d' % (c1, ss, cs, ext))
            j -= ext
    print('%d deletions, %d insertions, %d copies' % (n_del, n_ins, n_copy))


if __name__ == '__main__':
    seq1, point1 = load_test_seq('./motif.csv', int(sys.argv[1]))
    seq2, point2 = load_test_seq('./motif.csv', int(sys.argv[2]))
    print('Seq1 length: %d; Seq2 length: %d' % (len(seq1), len(seq2)))
    print('Points1: %d; Points2: %d' % (sum([1 for i in point1 if i != -1]), sum([1 for i in point2 if i != -1])))
    gary_benson(seq1, point1, seq2, point2)
