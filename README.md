# Electrical topology detection from telemetry

Recover the tree structure of an electrical installation (which feeder is
fed by which) using only power measurements from sensors installed on some of
its feeders. No wiring diagram is needed.

The result is the **induced tree** over the monitored feeders: each monitored
feeder is linked to its nearest *monitored* ancestor, skipping unmonitored
feeders in between.

```
Complete installation          Induced tree (what you get)

        D1*                            D1
   ┌────┼────┬────┐                 ┌───┴───┐
  D2*  D3   D4   D5                 D2      D9
 ┌─┴─┐   │    │    │                 │
D6  D7*  D8  D9*  D10               D7          (* = monitored)
```

## How it works

1. **Candidate parents.** Feeder `i` can be the parent of feeder `j` only if
   its power is at least `j`'s at every time step (a parent carries the power
   of its children).
2. **Scores.** Among the candidates, score each pair by the Pearson
   correlation of the *variations* of power between consecutive time steps.
   Working on variations removes most of the load profile shared by all
   feeders.
3. **Assignment.** Each feeder picks the best-scoring parent. The feeder with
   no candidate parent is the root (main feeder).

Two variants are implemented:

| Variant | Description |
|---|---|
| `base` | Choose the highest-correlation candidate for each feeder. This is the optimum of the assignment integer program, whose LP relaxation is already integral. |
| `capacity` | Same objective, plus the constraint that at every retained instant the powers of a feeder's children sum to at most the power of the feeder. This couples the choices, so it is solved as a mixed-integer linear program (HiGHS through SciPy). |

The method, with the formulation and its properties, is described in
`2025_topology_detection.tex`.

## Requirements

- Python 3.9+
- `numpy`, `pandas`, `scipy` (version 1.9 or later, for `scipy.optimize.milp`)

```
pip install numpy pandas scipy
```

## Input

A CSV file with one row per measurement and three columns, in this order:

| timestamp | feeder id | active power (W) |
|---|---|---|
| 2026-11-02T00:00:00 | F0 | 18450.2 |
| 2026-11-02T00:00:00 | F1 | 6210.7 |
| 2026-11-02T00:10:00 | F0 | 18102.9 |
| ... | ... | ... |

Assumptions:

- a single main feeder, and it is monitored;
- all values are non-negative (no injection, e.g. no solar panels);
- the sensors are synchronized.

Only timestamps where **every** feeder has a reading are used. If sensor
clocks are slightly offset, align them with `--round` (for example
`--round 10min`).

The file is expected to have a header row. Use `--no-header` if it does not.

## Usage

```
python detect_topology.py telemetry.csv
```

Options:

| Option | Default | Meaning |
|---|---|---|
| `--no-header` | off | the CSV has no header row |
| `--round FREQ` | none | floor timestamps to this frequency, e.g. `10min` |
| `--eps W` | 0 | tolerance in W for the dominance test and the capacity constraints |
| `--m-prime N` | 200 | number of evenly spaced instants used in the capacity constraints |
| `--time-limit S` | 60 | time limit of the capacity MILP, in seconds |
| `--out FILE` | `topology.csv` | output file |

## Output

The script prints the edges (`parent -> child`) with their correlation, and the
objective value, for each variant:

```
9 feeders, 1008 complete time steps (0 incomplete steps dropped)
root: F0

=== base ===
F0 -> F1   (rho = 0.615)
F1 -> F4   (rho = 0.800)
F4 -> F9   (rho = 0.511)
...
objective = 4.8805

=== capacity ===
...
feeders whose parent differs between variants: none
```

It also writes `topology.csv` with the columns `feeder`, `parent_base` and
`parent_capacity` (empty for the root).

Feeders on which the two variants disagree are the ones to look at first.

## Notes and limitations

- **One root required.** If the script finds zero or several feeders with no
  candidate parent, it stops: the main feeder may be missing, or the data
  inconsistent.
- **Identical series.** Two feeders with identical series dominate each other;
  the script keeps only the lower-indexed one as a possible parent, which
  avoids cycles.
- **Amount of data.** The method is designed for about a week of data at
  10-minute intervals or finer. With less data the dominance test rules out
  fewer candidates and the correlations are noisier, so treat the result with
  more caution.
- **Shared load profiles.** If many loads follow the same smooth profile (for
  example a common daily cycle), its variations remain in the data and can bias
  the correlations toward higher-level feeders.
- **Capacity instants.** The `m'` instants of the capacity variant are
  currently evenly spaced. This is a placeholder for a better selection rule.
- **Reactive power.** The paper describes how reactive power can strengthen the
  scores and provide a consistency check. The script does not use it yet.

## Files

| File | Content |
|---|---|
| `detect_topology.py` | the implementation (both variants) |
| `2025_topology_detection.tex` | the paper describing the method |
