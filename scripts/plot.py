from scipy import stats
import numpy as np

alpha = 0.5


def plot_depth_norm(data, axs):
    p = {
        "xlabel": "Mean-cell count",
        "ylabel": "Sum-cell count",
        # "xscale": "symlog",
        # "yscale": "symlog",
    }

    for idx, (ax, (title, matrix)) in enumerate(zip(axs, data.items())):
        if idx > 0:
            p.update({"ylabel": ""})
        x = data["raw"].mean(1)
        y = matrix.sum(1)
        ax.scatter(x, y, edgecolor="k", facecolor="#D43F3A", alpha=alpha)

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        minx = min(x)
        maxx = max(x)
        x = np.linspace(minx, maxx, 10)
        y = slope * x + intercept
        ax.plot(
            x,
            y,
            label=f"$r^2$ : {r_value**2:,.2f}\n$r$ : {r_value:,.2f}",
            color="k",
            linewidth=2,
        )
        ax.legend(fontsize=15)

        p.update({"title": title})
        ax.set(**p)

    return axs


def plot_depth_dist(data, axs):
    p = {
        "xlabel": "Cell counts",
        "ylabel": "Frequency",
    }
    # mn = min(data["raw"].sum())
    # mx = max

    for idx, (title, matrix) in enumerate(data.items()):
        ax = axs[idx]
        if idx > 0:
            p.update({"ylabel": ""})

        x = np.sort(matrix.sum(1))[::-1]

        mean = x.mean()
        stdev = np.sqrt(np.var(data["raw"].sum(1)))

        close = np.all(np.allclose(x, x[0]))

        if close:
            ax.hist([x[0]] * len(x), edgecolor="k", facecolor="#D43F3A")
        else:
            ax.hist(x, edgecolor="k", facecolor="#D43F3A")
        p.update({"title": title, "xlim": (mean - 5 * stdev, mean + 5 * stdev)})
        ax.set(**p)
    return axs


def plot_knee(data, axs):
    p = {
        "xscale": "linear",
        "yscale": "symlog",
        "xlabel": "Cell counts",
        "ylabel": "Cell rank",
    }

    for idx, (title, matrix) in enumerate(data.items()):
        if idx > 0:
            p.update({"ylabel": ""})
        ax = axs[idx]

        x = np.sort(matrix.sum(1))[::-1]
        y = np.arange(x.shape[0])

        ax.scatter(x, y, edgecolor="k", facecolor="#D43F3A", alpha=alpha)
        p.update({"title": title})
        ax.set(**p)
    return axs


def plot_pc_depth(data, pcs, axs):
    p = {
        "xlabel": "Summed-cell UMI counts",
        "ylabel": "PC1",
    }

    for idx, (title, pc) in enumerate(pcs.items()):
        if idx > 0:
            p.update({"ylabel": ""})
        ax = axs[idx]

        # raw sumed counts for all
        matrix = data["raw"]

        x = matrix.sum(1)
        y = pc[:, 0]
        ax.scatter(x, y, edgecolor="k", facecolor="#D43F3A", alpha=alpha)

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        minx = min(x)
        maxx = max(x)
        x = np.linspace(minx, maxx, 10)
        y = slope * x + intercept
        ax.plot(
            x,
            y,
            label=f"$r^2$ : {r_value**2:,.2f}\n$r$ : {r_value:,.2f}",
            color="k",
            linewidth=2,
        )
        ax.legend(fontsize=15)

        p.update({"title": title})
        ax.set(**p)
    return axs


def plot_mean_var(data, axs):
    p = {
        "xlabel": "Gene mean",
        "ylabel": "Gene variance",
        "xscale": "symlog",
        "yscale": "symlog",
        "ylim": (-1, 1e8),
    }

    for idx, (title, matrix) in enumerate(data.items()):
        if idx > 0:
            p.update({"ylabel": ""})
        ax = axs[idx]

        x = np.mean(matrix, axis=0)
        y = np.var(matrix, axis=0)
        ax.scatter(x, y, edgecolor="k", facecolor="#D43F3A", alpha=alpha)
        p.update({"title": title})
        ax.set(**p)
    return axs


def plot_monotone(data, rv, axs):

    p = {"xlabel": "Spearman r", "ylabel": "Frequency"}

    for idx, (label, m) in enumerate(data.items()):
        if idx > 0:
            p.update({"ylabel": ""})
        ax = axs[idx]
        x = rv[:, idx]

        ax.hist(x, edgecolor="k", facecolor="#D43F3A")
        p.update({"title": label})
        ax.set(**p)
    return axs


def plot_example_gene(data, axs):
    argmax = np.argmax(data["raw"].sum(0))
    p = {
        "xlabel": "Gene UMI counts",
        "ylabel": "Frequency",
    }

    for idx, (title, matrix) in enumerate(data.items()):
        if idx > 0:
            p.update({"ylabel": ""})
        ax = axs[idx]

        g = matrix[:, argmax]
        ax.hist(g, edgecolor="k", facecolor="#D43F3A")
        p.update({"title": title})
        ax.set(**p)
    return axs
