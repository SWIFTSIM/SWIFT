"""The reception slots on the sphere of directions, for 12, 8 and 16 slots."""

import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def frame(angles_deg=(24.0, 36.0, 66.0)) -> np.ndarray:
    """Return R = Rz Ry Rx, the rotation from the box to the slot frame."""
    ax, ay, az = np.radians(angles_deg)
    rx = np.array(
        [[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]]
    )
    ry = np.array(
        [[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]]
    )
    rz = np.array(
        [[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]]
    )
    return rz @ ry @ rx


def faces(K: int) -> list:
    """Return (slot, colour index, corners) of each cell, corners in frame coordinates."""
    e = np.eye(3)
    out = []
    if K == 12:
        for m in range(3):
            i, j = [a for a in range(3) if a != m]
            for sa in (0, 1):
                for sb in (0, 1):
                    si, sj = 2 * sa - 1, 2 * sb - 1
                    corner = lambda s: (s * e[m] + si * e[i] + sj * e[j]) / np.sqrt(3)
                    quad = [si * e[i], corner(1), sj * e[j], corner(-1)]
                    out.append((4 * m + 2 * sa + sb, m, quad))
    else:
        for o in range(8):
            s = np.array([1 if o >> k & 1 else -1 for k in range(3)], float)
            if K == 8:
                tri = [s[0] * e[0], s[1] * e[1], s[2] * e[2], s[2] * e[2]]
                out.append((o, bin(o).count("1") % 2, tri))
            else:
                mid = (s[0] * e[0] + s[1] * e[1]) / np.sqrt(2)
                out.append((o, 0, [s[0] * e[0], mid, s[2] * e[2], s[2] * e[2]]))
                out.append((o + 8, 1, [mid, s[1] * e[1], s[2] * e[2], s[2] * e[2]]))
    return out


def patch(quad: list, n: int = 14) -> list:
    """Split a spherical quadrilateral into small spherical patches."""
    a, b, c, d = [np.asarray(v, float) for v in quad]
    u = np.linspace(0, 1, n + 1)
    grid = np.array(
        [
            [
                (1 - s) * (1 - t) * a + s * (1 - t) * b + s * t * c + (1 - s) * t * d
                for t in u
            ]
            for s in u
        ]
    )
    grid /= np.linalg.norm(grid, axis=2)[..., None]
    return [
        [grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]]
        for i in range(n)
        for j in range(n)
    ]


def arc(p: np.ndarray, q: np.ndarray, n: int = 30) -> np.ndarray:
    """Return points on the great arc from p to q."""
    t = np.linspace(0, 1, n)[:, None]
    x = (1 - t) * p + t * q
    return x / np.linalg.norm(x, axis=1)[:, None]


def view_vector(elev: float, azim: float) -> np.ndarray:
    """Return the unit vector towards the viewer."""
    e, a = np.radians(elev), np.radians(azim)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def draw(ax, K: int, elev: float, azim: float, labels: bool, R: np.ndarray) -> None:
    """Draw the visible cells of one partition on a 3D axis."""
    colours = [st.SKY, st.ORANGE, st.GREEN]
    v = view_vector(elev, azim)
    polys, cols = [], []
    for k, ci, quad in faces(K):
        for p in patch(quad):
            if np.mean(p, axis=0) @ v > 0.02:
                polys.append(p)
                cols.append(colours[ci])
    ax.add_collection3d(
        Poly3DCollection(polys, facecolors=cols, edgecolors=cols, linewidths=0.5)
    )
    for k, ci, quad in faces(K):
        q = [np.asarray(x, float) / np.linalg.norm(x) for x in quad]
        for a, b in zip(q, q[1:] + q[:1]):
            if np.allclose(a, b):
                continue
            pts = arc(a, b) * 1.003
            vis = pts @ v > 0.0
            if vis.sum() > 1:
                ax.plot(*pts[vis].T, color=st.BLACK, lw=0.7)
        centre = np.mean(q, axis=0)
        centre /= np.linalg.norm(centre)
        if labels and centre @ v > 0.15:
            ax.text(
                *(1.04 * centre),
                str(k),
                ha="center",
                va="center",
                fontsize=8,
                zorder=10,
            )
    if labels:
        for i, name in enumerate(("y_0", "y_1", "y_2")):
            e = np.eye(3)[i]
            if e @ v > -0.2:
                ax.plot(*np.array([e, 1.22 * e]).T, color=st.BLACK, lw=0.9)
                ax.text(
                    *(1.42 * e + (0.0, 0.0, 0.16 if i == 0 else 0.0)),
                    f"${name}$",
                    ha="center",
                    va="center",
                )
        for i in range(3):
            for s in (1, -1):
                d = s * R[:, i]
                if d @ v > 0.1:
                    ax.scatter(
                        *d * 1.005, color=st.BLACK, s=7, depthshade=False, zorder=11
                    )
    u = np.cross(v, [0.0, 0.0, 1.0])
    u /= np.linalg.norm(u)
    w = np.cross(v, u)
    t = np.linspace(0, 2 * np.pi, 200)[:, None]
    ax.plot(*(np.cos(t) * u + np.sin(t) * w).T, color=st.BLACK, lw=0.7)
    ax.set_proj_type("ortho")
    ax.view_init(elev=elev, azim=azim)
    lim = 1.15 if labels else 1.0
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim, lim))
    ax.set_box_aspect((1, 1, 1), zoom=1.25 if labels else 1.3)
    ax.set_axis_off()


def main() -> None:
    """Draw the figure."""
    st.setup()
    R = frame()
    fig = plt.figure(figsize=(st.TEXT_WIDTH, 2.35))
    boxes = [
        (0.0, 0.0, 0.36, 1.0),
        (0.33, 0.0, 0.36, 1.0),
        (0.70, 0.5, 0.26, 0.48),
        (0.70, 0.0, 0.26, 0.48),
    ]
    views = [
        (22, 35, 12, True),
        (-22, 215, 12, True),
        (22, 35, 8, False),
        (22, 35, 16, False),
    ]
    titles = [
        "12 slots, front",
        "12 slots, back",
        "8 slots (octants)",
        "16 slots (half-octants)",
    ]
    for box, (el, az, K, lab), title in zip(boxes, views, titles):
        ax = fig.add_axes(box, projection="3d", computed_zorder=False)
        draw(ax, K, el, az, lab, R)
        ax.set_title(title, y=1.0 if lab else 0.95, fontsize=8)
    st.save(fig, "isrf_slots.pdf")


if __name__ == "__main__":
    main()
