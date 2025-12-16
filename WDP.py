import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def cost(a: str, b: str):
    return 1 if a == b else -1


def load_test_seq(fn: str, n: int):
    df = pd.read_csv(fn)
    return df['motif'][n-1], df['text'][n-1]


def view_matrix(pattern: str, text: str, dp: np, bt: list, fn: str):
    if dp.ndim != 2:
        print('Only support 2D matrix')
        return False
    n = len(text)
    m = len(pattern)
    height = 8
    width = int(height * (m / n)) + 1
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
    x_tick_labels = ['%s' % (pattern[i]) for i in range(m)]
    y_ticks = [i for i in range(1, n + 1)]
    y_tick_labels = ['%s' % (text[i]) for i in range(n)]
    plt.xticks(x_ticks, x_tick_labels)
    plt.yticks(y_ticks, y_tick_labels)
    # show dp values
    for y in range(0, n + 1):
        for x in range(0, m + 1):
            plt.text(x + 0.5, y + 0.5, '%d' % dp[y][x])
    plt.grid()
    # show backtrace
    for i in range(1, len(bt)):
        cy, cx = bt[i-1]
        py, px = bt[i]
        plt.arrow(cx + 0.5, cy + 0.5, px - cx, py - cy, head_width=0.2)
    plt.savefig("./%s.png" % fn)


def local_brute_force(pattern: str, text: str):
    n = len(text)
    m = len(pattern)
    dp = np.zeros((n + 1, n * m + 1), dtype=int)
    con = ''.join([pattern for i in range(n)])  # n copies
    for i in range(1, n + 1):
        for j in range(1, n * m + 1):
            v = dp[i-1][j] - 1
            h = dp[i][j-1] - 1
            d = dp[i-1][j-1] + cost(text[i-1], con[j-1])
            dp[i][j] = max(v, h, d, 0)
    print('Local brute-force DP: best score=%d' % dp.max())

    x, y = -1, -1
    for i in range(1, n + 1):
        # Display the backtrace path of the first period (optimal segment)
        for j in range(1, n * m + 1):
            if dp.max() == dp[i][j]:
                x, y = i, j
                break
        if x != -1 or y != -1:
            break
    print('Best score achieved at (%d, %d) with motif copied %d times' % (x, y, (y + m - 1) // m))
    # Keep the rendered DP matrix as tiny as possible
    tiny_dp = np.zeros((n + 1, y + 5), dtype=int)
    for i in range(0, n + 1):
        for j in range(0, y + 5):
            tiny_dp[i][j] = dp[i][j]
    con = con[0: y + 4]
    bt = [(x, y)]
    while x > 0 and y > 0 and dp[x][y] != 0:
        if dp[x][y] == dp[x][y-1] - 1:
            bt.append((x, y-1))
        elif dp[x][y] == dp[x-1][y] - 1:
            bt.append((x-1, y))
        elif dp[x][y] == dp[x-1][y-1] + cost(text[x-1], con[y-1]):
            bt.append((x-1, y-1))
        x, y = bt[-1]
    while x > 0 and dp[x][y] != 0:
        bt.append((x-1, y))
        x, y = bt[-1]
    while y > 0 and dp[x][y] != 0:
        bt.append((x, y-1))
        x, y = bt[-1]
    view_matrix(con, text, tiny_dp, bt, 'local_force')


def local_wraparound(pattern: str, text: str):
    n = len(text)
    m = len(pattern)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        dp[i-1][0] = dp[i-1][m]  # It is correct for local alignment
        # First pass: the optimal alignments ending with a diagonal or vertical transfer are correctly calculated
        for j in range(1, m + 1):
            v = dp[i-1][j] - 1
            h = dp[i][j-1] - 1
            d = dp[i-1][j-1] + cost(text[i-1], pattern[j-1])
            dp[i][j] = max(v, h, d, 0)
        # Second pass: the optimal alignments ending with a horizontal transfer are correctly calculated
        dp[i][0] = dp[i][m]
        for j in range(1, m + 1):
            dp[i][j] = max(dp[i][j], dp[i][j-1] - 1)
    print('Local wraparound DP: best score=%d' % dp.max())

    x, y = -1, -1
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if dp.max() == dp[i][j]:
                x, y = i, j
                break
        if x != -1 or y != -1:
            break
    bt = [(x, y)]
    while x > 0 and y > 0 and dp[x][y] != 0:
        # If it contains two periods, some backtrace arrows might overlap on the same row.
        # For example, two gaps between GACTG and ACTACTACTACT cross the segment boundary.
        if y == 1:
            if dp[x][1] == dp[x][m] - 1:
                bt.append((x, m))
            elif dp[x][1] == dp[x-1][m] - 1:
                bt.append((x-1, 1))
            else:
                bt.append((x-1, m))
        elif dp[x][y] == dp[x][y-1] - 1:
            bt.append((x, y-1))
        elif dp[x][y] == dp[x-1][y] - 1:
            bt.append((x-1, y))
        elif dp[x][y] == dp[x-1][y-1] + cost(text[x-1], pattern[y-1]):
            bt.append((x-1, y-1))
        x, y = bt[-1]
    while x > 0 and dp[x][y] != 0:
        bt.append((x-1, y))
        x, y = bt[-1]
    while y > 0 and dp[x][y] != 0:
        bt.append((x, y-1))
        x, y = bt[-1]
    view_matrix(pattern, text, dp, bt, 'local_wraparound')


def global_brute_force(pattern: str, text: str):
    n = len(text)
    m = len(pattern)
    dp = np.zeros((n + 1, n * m + 1), dtype=int)
    con = ''.join([pattern for i in range(n)])  # n copies
    dp[0][0] = 0
    for i in range(1, n + 1):
        dp[i][0] = -i
    for j in range(1, n * m + 1):
        dp[0][j] = -j
    for i in range(1, n + 1):
        for j in range(1, n * m + 1):
            v = dp[i-1][j] - 1
            h = dp[i][j-1] - 1
            d = dp[i-1][j-1] + cost(text[i-1], con[j-1])
            dp[i][j] = max(v, h, d)

    print('Global brute-force DP: best score=%d' % dp.max())


def global_wraparound(pattern: str, text: str):
    n = len(text)
    m = len(pattern)
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
            if j > 1:
                d = dp[i-1][j-1] + cost(text[i-1], pattern[j-1])
            else:
                d = max(dp[i-1][m], dp[i-1][0]) + cost(text[i-1], pattern[j-1])
            dp[i][j] = max(v, h, d)
        dp[i][0] = dp[i][m]
        for j in range(1, m + 1):
            dp[i][j] = max(dp[i][j], dp[i][j-1] - 1)
    print('Global wraparound DP: best score=%d' % dp.max())


if __name__ == '__main__':
    motif, target = load_test_seq('./motif.csv', 2)
    local_brute_force(motif, target)
    local_wraparound(motif, target)
