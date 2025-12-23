import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt


def load_test_seq(fn: str, n: int):
    df = pd.read_csv(fn)
    return df['text'][n-1]


def visualize(que: str, ref: str, bt: np):
    n, m = len(que), len(ref)
    height = 8
    width = round(height * m / n)
    plt.figure(figsize=(width, height), dpi=350)
    ax = plt.gca()
    # axes scale equally in the gridded figure
    ax.set_aspect('equal', adjustable='box')
    # put origin point on left-top corner
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.xlim(0, m + 1)
    plt.ylim(n + 1, 0)
    i, j = n, m
    DIAGNOAL, VERTICAL, HORIZONTAL = 0, 1, 2
    COLOR_GAP, COLOR_DIA = 'orange', 'green'
    n_del, n_ins = 0, 0
    while i > 0 or j > 0:
        x0, y0 = j, i
        if bt[i][j] == DIAGNOAL:
            x1, y1 = j-1, i-1
            plt.plot([x0, x1], [y0, y1], color=COLOR_DIA)
            i, j = i-1, j-1
        elif bt[i][j] == VERTICAL:
            x1, y1 = j, i-1
            plt.plot([x0, x1], [y0, y1], color=COLOR_GAP)
            i, j = i-1, j
            n_del += 1
        else:
            x1, y1 = j-1, i
            plt.plot([x0, x1], [y0, y1], color=COLOR_GAP)
            i, j = i, j-1
            n_ins += 1

    print('%d deletions and %d insertions' % (n_del, n_ins))
    ax.set_title('Global Smith-Waterman Alignment')
    ax.set_xlabel('Sequence A (len=%d)' % m)
    ax.set_ylabel('Sequence B (len=%d)' % n)
    plt.grid()
    plt.plot([0, 0], [0, 0], color=COLOR_GAP, label='indel')
    plt.plot([0, 0], [0, 0], color=COLOR_DIA, label='mis/match')
    plt.legend()
    plt.savefig("./gsw_matrix.png")


def smith_waterman(que: str, ref: str):
    n, m = len(que), len(ref)
    mat_s, mis_p, gap_o, gap_e = 1, -4, -6, -1
    E = np.zeros((n + 1, m + 1), dtype=int)
    H = np.zeros((n + 1, m + 1), dtype=int)
    bt = np.zeros((n + 1, m + 1), dtype=int)
    DIAGNOAL, VERTICAL, HORIZONTAL = 0, 1, 2
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
            H[i][j] = max(M, E[i][j], F)
            if H[i][j] == F:
                bt[i][j] = HORIZONTAL
            elif H[i][j] == M:
                bt[i][j] = DIAGNOAL
            else:
                bt[i][j] = VERTICAL
    print('Score: %d' % H[n][m])
    visualize(que, ref, bt)


if __name__ == '__main__':
    seq1 = load_test_seq('./motif.csv', int(sys.argv[1]))
    seq2 = load_test_seq('./motif.csv', int(sys.argv[2]))
    print('Seq1 length: %d; Seq2 length: %d' % (len(seq1), len(seq2)))
    smith_waterman(seq1, seq2)
