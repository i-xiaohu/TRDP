import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys


def load_test_seq(fn: str, n: int):
    df = pd.read_csv(fn)
    return df['text'][n-1]


def dotplot_dp(seq: str, dp: np, fn: str):
    n = m = len(seq)
    height = width = 8
    plt.figure(figsize=(width, height), dpi=350)
    # axes scale equally in the gridded figure
    plt.gca().set_aspect('equal', adjustable='box')
    # put origin point on left-top corner
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlim(0, m + 1)
    plt.ylim(n + 1, 0)
    # show pattern and text on axes
    x_ticks = [i for i in range(1, m + 1)]
    x_tick_labels = ['%s [%d]' % (seq[i], i + 1) for i in range(m)]
    y_ticks = [i for i in range(1, n + 1)]
    y_tick_labels = ['[%d] %s' % (i + 1, seq[i]) for i in range(n)]
    plt.xticks(x_ticks, x_tick_labels)
    plt.yticks(y_ticks, y_tick_labels)
    # FIXME: how to prevent the main diagonal from dominating the color rendering
    black_color = dp.max()
    for y in range(0, n + 1):
        for x in range(0, m + 1):
            cell_color = round(255 * (1 - dp[y][x] / black_color))
            hex_color = f"#{cell_color:02X}{cell_color:02X}{cell_color:02X}"
            rect = patches.Rectangle((x, y), 1, 1, color=hex_color)
            plt.gca().add_patch(rect)
    # show dp values
    for y in range(0, n + 1):
        for x in range(0, m + 1):
            plt.text(x + 0.5, y + 0.5, '%d' % dp[y][x])
    plt.grid()
    plt.savefig("./%s.png" % fn)


# The result is like mosaic without clear pattern indicating tandem repeats
def dotplot(seq: str, fn: str):
    n = m = len(seq)
    height = width = 8
    plt.figure(figsize=(width, height), dpi=350)
    # axes scale equally in the gridded figure
    plt.gca().set_aspect('equal', adjustable='box')
    # put origin point on left-top corner
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xlim(0, m)
    plt.ylim(n, 0)
    for y in range(0, n):
        for x in range(0, m):
            if seq[y] == seq[x]:
                rect = patches.Rectangle((x, y), 1, 1, color='green')
                plt.gca().add_patch(rect)
    plt.grid()
    plt.savefig("./%s.png" % fn)


def self_alignment(seq: str):
    n = len(seq)
    dp = np.zeros((n + 1, n + 1), dtype=int)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            v = dp[i-1][j] - 1
            h = dp[i][j-1] - 1
            d = dp[i-1][j-1] + (1 if seq[i-1] == seq[j-1] else 0)
            dp[i][j] = max(v, h, d, 0)

    # dotplot(seq, 'dotplot')
    dotplot_dp(seq, dp, 'dotplot_dp')


def simple_search(text: str):
    # It is easy to identify all perfect tandem repeats.
    # However, it is much harder to define them if they contain mismatches.
    # One aiming to identify all approximate tandem repeats must precisely define the property of them,
    # with notations like the minimum repeat unit length and the maximum divergence.
    # Without a clear definition, an exact algorithm for the objective would be impossible.

    # Even we do have a standard for testing whether two substrings are considered tandemly repeat,
    # it doesn't apply to a series of repeat units, in which each two adjacent units are similar but the
    # head and tail units share no common (for example, AATT | AATC | AAGC | GAGC | GCGC).
    # A reasonable method is to pairwise compare repeat units under the same criteria, but it would be time-consuming.

    # Another phenomenon that further complicates the problem is the nested repeats, meaning the repeat unit itself
    # is a repeat. At first I thought the problem was trivial, because we can simply prohibit such repeat units.
    # Unfortunately I discovered that it is only workable for perfect repeats, not for biologically meaningful ones.
    # Repeat units of different size may yield different divergence levels. A classical example is the monomers and HORs
    # in centromeres. Which type of repeats will be discovered is up to the preference for unit length or similarity.
    # Still, we can not guarantee the completeness of final results even if we prioritize smaller unit length.
    # Consider this repeat ACT | AGT | ACT | AGT with 1 mismatch for each pair of units, another possible combination
    # is ACTAGT | ACTAGT, leading to a perfect repeat.

    # Under such complex background, I decide to first implement a somehow brute-force algorithm with global wraparound
    # DP as a subroutine. I hope this approach will exactly identifies all tandem repeats of our concern.
    # For a string S length of n, I traverse the string backward to find potential repeats. In the iteration of i, I
    # enumerate the unit length k and determine if a suffix of S[1:i] can be formed by multiple copies of S[i-k+1, i].
    # This can be done with global WDP in O(N) time. If the alignment score between S[q:i] and S[i-k+1, i] is above
    # a threshold, i.e., their divergence is below a certain value, I report the result.To avoid redundant computation
    # and reports, I skip the substring S[q:i] and start the process at q - 1 if S[q:i] is the shortest suffix among
    # identified repeats (check the correctness). In the worst case, the time complexity is O(N^3). But banded DP will
    # easily reduce it to O(k * N^2) where k is the maximum allowed mismatches within repeats. If the input string is
    # abundant with repeats, the program will run faster by skipping repeats.

    # By the way, the self-alignment techniques like Z and anti-Z transfer are likely the same approach described above
    # but concealed in DP matrix. Consider all repeats as sub-diagonals in the matrix, all essential computations also
    # occur in my program, except for the pointless competition with the main diagonal. The anti-Z is remarkably similar
    # to my jump operation.

    min_unit_len = 3
    min_repeat_len = 10
    max_motif_len = 100
    n = len(text)
    i = 0
    while i + min_repeat_len < n:
        has_repeat = False
        jump = 0
        for k in range(min_unit_len, i // 2):
            # Alignment between S[1:i] and S[i-k+1, i]
            unit = '@' + text[i: i-k: -1]
            substring = '@' + text[i-k:0:-1]
            # print('Alignment between %s and %s' % (unit, substring))
            m = len(substring) - 1
            dp = np.zeros((m + 1, k + 1), dtype=int)
            for x in range(1, m + 1):
                dp[x][0] = x
            for y in range(1, k + 1):
                dp[0][y] = y
            for x in range(1, m + 1):
                # Pass 1
                for y in range(1, k + 1):
                    h = dp[x][y-1] + 1
                    v = dp[x-1][y] + 1
                    if y == 1:
                        d = min(dp[x-1][k], dp[x-1][0]) + (0 if substring[x] == unit[y] else 1)
                    else:
                        d = dp[x-1][y-1] + (0 if substring[x] == unit[y] else 1)
                    dp[x][y] = min(h, v, d)
                # Pass 2
                dp[x][0] = dp[x][k]
                for y in range(1, k + 1):
                    dp[x][y] = min(dp[x][y], dp[x][y-1] + 1)
            # print('Global score: %d' % dp[m][k])
            for x in range(m, 0, -1):
                if dp[x][k] <= x * max_divergence:
                    print('unit: %s' % unit[1:])
                    print('substring: %s' % substring[1:x+1])
                    print('dp(%d,%d): %d' % (x, k, dp[x][k]))
                    jump = max(jump, i - k - x)
                    has_repeat = True
                    break
        if has_repeat:
            i = jump
        else:
            i -= 1


if __name__ == '__main__':
    target_seq = load_test_seq('motif.csv', int(sys.argv[1]))
    simple_search(target_seq)
