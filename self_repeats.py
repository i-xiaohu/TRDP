import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt


def load_test_seq(fn: str, n: int):
    df = pd.read_csv(fn)
    return df['text'][n-1]


def duplicate(motif: str, target: str):
    n = len(target)
    m = len(motif)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    dp[0][0] = 0
    mat_s, mis_p, gap_p = 1, -3, -2
    for i in range(1, n + 1):
        dp[i][0] = i * gap_p
    for j in range(1, m + 1):
        dp[0][j] = j * gap_p
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            v = dp[i-1][j] + gap_p
            h = dp[i][j-1] + gap_p
            d = dp[i-1][j-1] + (mat_s if target[i-1] == motif[j-1] else gap_p)
            dp[i][j] = max(v, h, d)
        dp[i][0] = dp[i][m]
        for j in range(1, m + 1):
            dp[i][j] = max(dp[i][j], dp[i][j-1] + gap_p)
        # FIXME: Is it correct?
        if dp[i].max() < 0:
            break
    max_score, ext = n * m * min(gap_p, mis_p), -1
    for i in range(1, n + 1):
        # FIXME: I ignore incomplete - potentially optimal - copy for easy implementation
        if dp[i][m] > max_score:
            max_score = dp[i][m]
            ext = i
    if max_score <= 0:
        max_score = 0
        ext = 0
    else:
        x, y, n_copy = ext, m, 0
        bt = [(x, y)]
        while x > 0 and y > 0:
            if y == 1:
                if dp[x][1] == dp[x][m] + gap_p:
                    n_copy += 1
                    bt.append((x, m))
                elif dp[x][1] == dp[x-1][1] + gap_p:
                    bt.append((x-1, 1))
                else:
                    n_copy += 1
                    if x == 1:
                        bt.append((x-1, 0))
                    else:
                        bt.append((x-1, m))
            elif dp[x][y] == dp[x][y-1] + gap_p:
                bt.append((x, y-1))
            elif dp[x][y] == dp[x-1][y] + gap_p:
                bt.append((x-1, y))
            else:
                bt.append((x-1, y-1))
            x, y = bt[-1]
        max_score += n_copy
    return max_score, ext


def visualize(text: str, bt):
    n = len(text)
    width = height = 8
    plt.figure(figsize=(width, height), dpi=350)
    ax = plt.gca()
    # axes scale equally in the gridded figure
    ax.set_aspect('equal', adjustable='box')
    # put origin point on left-top corner
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.xlim(0, n + 1)
    plt.ylim(n + 1, 0)

    copy_label, non_rep_label, tr_label = True, True, True
    i = n - 1
    x0, y0 = i, i
    while i >= 0:
        p = bt[i]
        if p != i-1:
            motif = p[2][::-1]
            mlen, tlen = len(motif), i - p[0]
            x1, y1 = p[0], p[0]
            if tr_label:
                plt.plot([x0, x1], [y0, y1], color='red', label='tandem repeat')
                tr_label = False
            else:
                plt.plot([x0, x1], [y0, y1], color='red')
            plt.text(x1, y0+2, 'motif=%s, repeat_range=[%d,%d)' % (motif, p[0]+1, i+1))
            plt.plot([x0, x0-tlen+mlen], [y0, y0], color='red', linestyle='--')
            repeat_time = tlen // mlen
            px, py = x0-tlen+mlen, y0
            for k in range(0, repeat_time - 1):
                nx, ny = px - mlen, py
                plt.plot([px, nx], [py, ny], color='red', linestyle='--')
                px, py = nx, ny
                nx, ny = px + mlen, py - mlen
                if copy_label:
                    plt.plot([px, nx], [py, ny], color='green', label='copy')
                    copy_label = False
                else:
                    plt.plot([px, nx], [py, ny], color='green')
                px, py = nx, ny
            x0, y0 = x1, y1
            i = p[0]
        else:
            x1, y1 = i - 1, i - 1
            if non_rep_label:
                plt.plot([x0, x1], [y0, y1], color='blue', label='non-repetitive')
                non_rep_label = False
            else:
                plt.plot([x0, x1], [y0, y1], color='blue')
            i -= 1
            x0, y0 = x1, y1
    ax.set_title('Tandem Repeat Self Alignment')
    x_ticks = [i * 10 for i in range(0, n // 10 + 1)]
    y_ticks = [i * 10 for i in range(0, n // 10 + 1)]
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xlabel('Sequence S (len=%d)' % n)
    ax.set_ylabel('Sequence S (len=%d)' % n)
    plt.grid()
    plt.legend()
    plt.savefig("./self_matrix.png")


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
    min_motif_len = 3
    max_motif_len = 100
    n = len(text)
    dp = np.zeros(n, dtype=int)
    bt = dict()
    for i in range(0, min_repeat_len):
        dp[i] = i + 1
        bt[i] = i - 1
    # TODO: the loop is O(N^2 * L^2)
    for i in range(min_repeat_len, n):
        print('i = %d' % i)
        dp[i] = dp[i-1] + 1
        bt[i] = i - 1
        activate_early_stop = False
        fail_streak = 0
        for l in range(min_motif_len, min(max_motif_len, i // 2)):
            motif = text[i - l + 1: i + 1]
            target = text[0: i - l + 1]
            motif = motif[::-1]
            target = target[::-1]
            # print('motif=%s [%d,%d), target=%s [%d,%d)' % (motif, i-l+1, i+1, target, 0, i-l+1))
            w_score, ext = duplicate(motif, target)
            if w_score == 0 or ext + l < min_repeat_len:
                fail_streak += 1
                continue
            if activate_early_stop and fail_streak > 10:
                break
            fail_streak = 0
            # print('w_score=%d, ext=%d, start=%d' % (w_score, ext, i - l - ext))
            # No reward for itself (it is not a copy)
            w = w_score + l + open_penalty
            if i - l - ext == -1:  # Extend to the start
                if w >= dp[i]:
                    activate_early_stop = True
                    dp[i] = w
                    bt[i] = (-1, w_score, motif, target)
            else:
                if dp[i - l - ext] + w > dp[i]:
                    activate_early_stop = True
                    dp[i] = dp[i - l - ext] + w
                    bt[i] = (i - l - ext, w_score, motif, target)

    n_repeats = 0
    i = n - 1
    while i >= 0:
        p = bt[i]
        if p != i-1:
            n_repeats += 1
            print('Repeat %d: motif=%s, target=%s [%d,%d), score=%d' % (n_repeats, p[2][::-1], p[3][::-1], p[0]+1, i+1, p[1]))
            i = p[0]
        else:
            i -= 1

    visualize(text, bt)


if __name__ == '__main__':
    seq = load_test_seq('./motif.csv', int(sys.argv[1]))
    print('sequence length: %d' % (len(seq)))
    solve(seq)
