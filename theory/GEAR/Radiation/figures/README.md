# Figures for the GEAR radiation theory document

Schematics of the HII neighbour-search mechanics in
`src/runner_radiation_feedback.c`. Each `figN_*.py` script regenerates its
`figN_*.pdf` (and a `.png` preview) standalone; `style.py` holds the shared
palette/layout helpers.

Regenerate all (needs matplotlib):

```bash
for f in fig*.py; do python3 "$f"; done
```

or run the document build with figure regeneration: `../run.sh --figures`.

- `fig1_traversal_structure` — star-centric traversal: self link descends the
  star's own `radiation_level` region, pair links descend the up-to-26
  neighbour regions; no neighbour-of-neighbour hops.
- `fig2_case1_independent_merge` — a radius covering several wired leaf cells:
  each leaf is tested independently and candidates merge into one global
  distance-sorted buffer (nearest-`max_ngbs` across the union).
- `fig3_sorted_scan_nonadjacent` — sorted scan against a non-adjacent
  receiver: sort entries and `di_star` live on the same global 1D axis;
  adjacency never enters. Includes the periodic-wrap panel.
- `fig4_any_axis_conservative` — a 1D projection never exceeds the 3D
  distance, so the slab window is conservative along any sorted axis.
- `fig5_case2_stencil_gap` — an unwired cell two regions away is structurally
  invisible; the three guards (split gate, rebuild criterion, in-pass clamp)
  bound the gap and the next rebuild coarsens `radiation_level`.
- `fig6_contrast_dopair_vs_radiation` — classic two-sided lockstep DOPAIR
  (adjacency-bound) vs this one-sided star-centric scan (receiver's sorts
  only).

All six figures are included in `02_task_graph.tex`: `fig1` and `fig2` in
"Run-time recursion and neighbour search", `fig3`, `fig4`, and `fig6` in
"Sorted gather mechanics", and `fig5` in "Sorted gather mechanics" alongside
the two-regions-away gap it illustrates.
