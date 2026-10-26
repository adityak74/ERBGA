"""Regenerate the README figures from a real ERBGA run.

    python docs/make_figures.py

Writes docs/karate-light.svg and docs/karate-dark.svg: Zachary's karate club
coloured by the communities ERBGA finds (modularity, seed 0). The layout is a
seeded Fruchterman-Reingold, so the output is deterministic.
"""

from pathlib import Path

import numpy as np

from erbga import detect_communities
from erbga.api import load_graph

OUT = Path(__file__).resolve().parent
W, H, PAD = 760, 470, 54

# Validated as a 4-colour all-pairs set in both modes (dataviz validate_palette).
# Colour is never the only cue: each community also has its own shape and label.
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", edge="#a3a29b",
                  series=["#2a78d6", "#eda100", "#e87ba4", "#008300"]),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", edge="#5f5e58",
                 series=["#3987e5", "#c98500", "#d55181", "#008300"]),
}
SHAPES = ["circle", "square", "triangle", "diamond"]
FONT = "-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"


def layout(n, src, dst, seed=3, iters=600):
    rng = np.random.default_rng(seed)
    pos = rng.uniform(-1, 1, (n, 2))
    k = np.sqrt(4.0 / n)
    temp = 0.2
    for _ in range(iters):
        delta = pos[:, None, :] - pos[None, :, :]
        dist = np.maximum(np.linalg.norm(delta, axis=2), 1e-3)
        force = (delta / dist[..., None] * (k * k / dist)[..., None]).sum(axis=1)
        d = pos[src] - pos[dst]
        dl = np.maximum(np.linalg.norm(d, axis=1), 1e-3)
        pull = d / dl[:, None] * (dl * dl / k)[:, None]
        np.add.at(force, src, -pull)
        np.add.at(force, dst, pull)
        length = np.maximum(np.linalg.norm(force, axis=1), 1e-9)
        pos += force / length[:, None] * np.minimum(length, temp)[:, None]
        temp *= 0.992
    pos -= pos.min(axis=0)
    pos /= pos.max(axis=0)
    return pos


def marker(shape, x, y, r, fill, ring):
    common = f'fill="{fill}" stroke="{ring}" stroke-width="2"'
    if shape == "circle":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" {common}/>'
    if shape == "square":
        s = r * 1.7
        return (f'<rect x="{x - s / 2:.1f}" y="{y - s / 2:.1f}" width="{s:.1f}" '
                f'height="{s:.1f}" rx="2" {common}/>')
    if shape == "triangle":
        h = r * 1.25
        pts = f"{x:.1f},{y - h:.1f} {x + h:.1f},{y + h * 0.8:.1f} {x - h:.1f},{y + h * 0.8:.1f}"
        return f'<polygon points="{pts}" {common}/>'
    h = r * 1.35
    pts = f"{x:.1f},{y - h:.1f} {x + h:.1f},{y:.1f} {x:.1f},{y + h:.1f} {x - h:.1f},{y:.1f}"
    return f'<polygon points="{pts}" {common}/>'


def render(theme, g, labels, pos, fitness):
    t = THEMES[theme]
    top, bottom = 78, H - 62
    xy = np.column_stack([PAD + pos[:, 0] * (W - 2 * PAD), top + pos[:, 1] * (bottom - top)])
    # communities ordered by size so slot 1 is the largest
    order = np.argsort(-np.bincount(labels), kind="stable")
    slot = np.empty_like(order)
    slot[order] = np.arange(len(order))
    sizes = np.bincount(labels)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'viewBox="0 0 {W} {H}" role="img" aria-labelledby="t d" '
             f'font-family="{FONT}">',
             '<title id="t">Karate club communities found by ERBGA</title>',
             f'<desc id="d">Zachary\'s karate club, 34 members and 78 friendships, split by ERBGA '
             f'into {len(sizes)} connected communities with modularity {fitness:.4f}, the proven '
             f'optimum. Solid edges are inside a community, dashed edges cross between '
             f'communities.</desc>',
             f'<rect width="{W}" height="{H}" rx="12" fill="{t["surface"]}"/>']
    for e in range(g.m):
        a, b = g.src[e], g.dst[e]
        inside = labels[a] == labels[b]
        dash = "" if inside else ' stroke-dasharray="4 4"'
        opacity = "0.9" if inside else "0.55"
        parts.append(f'<line x1="{xy[a, 0]:.1f}" y1="{xy[a, 1]:.1f}" x2="{xy[b, 0]:.1f}" '
                     f'y2="{xy[b, 1]:.1f}" stroke="{t["edge"]}" stroke-width="1.5" '
                     f'stroke-opacity="{opacity}"{dash}/>')
    for v in range(g.n):
        s = slot[labels[v]]
        r = 7 + 0.45 * g.degree[v]
        parts.append(marker(SHAPES[s], xy[v, 0], xy[v, 1], r, t["series"][s], t["surface"]))

    # legend row: shape + colour + text label for each community
    y = H - 26
    column = (W - 2 * PAD) / len(order)
    for s, c in enumerate(order):
        x = PAD + s * column
        parts.append(marker(SHAPES[s], x + 7, y - 4, 6, t["series"][s], t["surface"]))
        text = f"Community {s + 1} ({sizes[c]})"
        parts.append(f'<text x="{x + 20}" y="{y}" font-size="13" fill="{t["ink2"]}">{text}</text>')
    parts.append(f'<text x="{PAD}" y="30" font-size="15" font-weight="600" fill="{t["ink"]}">'
                 f'Zachary\'s karate club · modularity {fitness:.4f} (proven optimum)</text>')
    parts.append(f'<text x="{PAD}" y="50" font-size="12" fill="{t["ink2"]}">'
                 f'Solid edges stay inside a community; dashed edges cross between them. '
                 f'Every community is connected.</text>')
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main():
    g = load_graph("karate")
    result = detect_communities(g, "modularity", seed=0)
    pos = layout(g.n, g.src, g.dst)
    for theme in THEMES:
        path = OUT / f"karate-{theme}.svg"
        path.write_text(render(theme, g, result.labels, pos, result.fitness))
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
