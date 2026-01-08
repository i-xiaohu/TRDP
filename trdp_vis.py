import sys
import matplotlib.pyplot as plt


def visualize(test_id: int):
    points = []
    with open('self_%d.txt' % test_id, 'r') as f:
        for line in f:
            x, y = line.strip().split()
            x, y = int(x), int(y)
            points.append((x, y))

    width = height = 8
    plt.figure(figsize=(width, height), dpi=350)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    plt.plot([j for (i, j) in points], [i for (i, j) in points], color='black')
    plt.savefig('self_%d.png' % test_id)


if __name__ == '__main__':
    visualize(int(sys.argv[1]))
