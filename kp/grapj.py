import matplotlib.pyplot as plt


def read_p(name='input.txt'):
    with open(name) as f:
        lines = f.readlines()
        n = int(lines[0])
        pts = map(
            lambda x: (int(x[0]), int(x[1])),
            map(
                lambda x: x.strip().split(),
                lines[1:n + 1]
            )
        )
        return list(pts)


def read_l(name='output.txt'):
    with open(name) as f:
        lines = f.readlines()
        n = int(lines[0])
        ls = map(
            lambda x: ((int(x[0]), int(x[1])), (int(x[2]), int(x[3]))),
            map(
                lambda x: x.strip().split(),
                lines[1:n + 1]
            )
        )
        return list(ls)


def plot_p(p, show=True):
    p = p[:]
    p.append(p[0])
    xs, ys = zip(*p)
    plt.plot(xs, ys)
    if show:
        plt.show()


def plot_l(l, show=True):
    l = l[:]
    for (a, b), (c, d) in l:
        plt.plot([a, c], [b, d], color='r')
    if show:
        plt.show()


if __name__ == '__main__':
    plot_p(read_p('input.txt'), show=False)
    plot_l(read_l('output.txt'))
    # plot_p(read_p('input_pub.txt'), show=False)
    # plot_l(read_l())
