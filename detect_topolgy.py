"""Detect the induced tree topology of monitored feeders from power telemetry.

Input CSV, one row per measurement: timestamp, feeder id, active power (W).
Assumptions: a single main feeder (the root, monitored), all values >= 0.

Both variants of the method are run:
  base      : each feeder takes, among the feeders that dominate it at every
              time step, the parent with the highest correlation of power
              variations. The integer program decouples by feeder, so its LP
              relaxation is integral and the optimum is a per-feeder argmax.
  capacity  : same objective, plus the constraint that at every retained
              instant the powers of the children of a feeder sum to at most
              the power of that feeder. Solved as a MILP (HiGHS via SciPy).

Usage: python detect_topology.py telemetry.csv [--eps 5] [--round 10min]
"""
import argparse

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


def load_matrix(path, has_header, round_freq):
    """Read the CSV and return (feeder ids, matrix C of shape n x m)."""
    df = pd.read_csv(path, header=0 if has_header else None,
                     usecols=[0, 1, 2], names=["timestamp", "feeder", "power"])
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    if round_freq:
        df["timestamp"] = df["timestamp"].dt.floor(round_freq)
    df["feeder"] = df["feeder"].astype(str)
    table = df.pivot_table(index="feeder", columns="timestamp",
                           values="power", aggfunc="mean")
    n_all = table.shape[1]
    table = table.dropna(axis=1, how="any")  # keep instants seen by every feeder
    print(f"{table.shape[0]} feeders, {table.shape[1]} complete time steps "
          f"({n_all - table.shape[1]} incomplete steps dropped)")
    if table.shape[1] < 3:
        raise SystemExit("Not enough synchronized data (try --round 10min).")
    return table.index.to_numpy(), table.to_numpy(dtype=float)


def candidate_parents(C, eps):
    """dom[i, j] is True if feeder i can be the parent of feeder j, i.e.
    C[i, t] >= C[j, t] - eps for all t."""
    n = C.shape[0]
    dom = np.zeros((n, n), dtype=bool)
    for i in range(n):
        dom[i] = np.all(C[i][None, :] >= C - eps, axis=1)
    np.fill_diagonal(dom, False)
    # Identical series dominate each other (cycles): keep only i < j as parent.
    mutual = dom & dom.T
    rows, cols = np.indices((n, n))
    dom[mutual & (rows > cols)] = False
    return dom


def correlations(C):
    """Pearson correlation between the power variations of every feeder pair."""
    dC = np.diff(C, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.corrcoef(dC)
    return np.nan_to_num(rho, nan=0.0)  # constant series -> correlation 0


def solve_base(dom, rho):
    """Base variant: per-feeder argmax of the correlation among candidates.
    parent[j] = -1 for the root."""
    n = dom.shape[0]
    score = np.where(dom, rho, -np.inf)
    parent = np.full(n, -1)
    for j in range(n):
        if dom[:, j].any():
            parent[j] = np.argmax(score[:, j])
    return parent


def solve_capacity(C, dom, rho, m_prime, eps, time_limit):
    """Capacity variant: MILP with sum of children <= parent at m' instants.
    Returns the parent vector, or None if the solver finds no solution."""
    n, m = C.shape
    t_idx = np.unique(np.linspace(0, m - 1, min(m_prime, m)).astype(int))
    Cp = C[:, t_idx]
    mp = Cp.shape[1]

    pairs = np.argwhere(dom)  # one binary variable x(i, j) per row (i, j)
    nv = len(pairs)
    cost = -rho[pairs[:, 0], pairs[:, 1]]  # milp minimizes

    # Each non-root feeder has exactly one parent.
    children = np.unique(pairs[:, 1])
    child_row = np.searchsorted(children, pairs[:, 1])
    a_one = coo_matrix((np.ones(nv), (child_row, np.arange(nv))),
                       shape=(len(children), nv)).tocsr()

    # Capacity: for each parent i and instant t, sum_j x(i,j) C'(j,t) <= C'(i,t) + eps.
    rows = (pairs[:, 0][:, None] * mp + np.arange(mp)[None, :]).ravel()
    cols = np.repeat(np.arange(nv), mp)
    data = Cp[pairs[:, 1]].ravel()
    a_cap = coo_matrix((data, (rows, cols)), shape=(n * mp, nv)).tocsr()
    ub = (Cp + eps).ravel()
    keep = np.diff(a_cap.indptr) > 0  # drop rows of feeders with no candidate child
    a_cap, ub = a_cap[keep], ub[keep]

    res = milp(c=cost,
               constraints=[LinearConstraint(a_one, 1, 1),
                            LinearConstraint(a_cap, -np.inf, ub)],
               integrality=np.ones(nv), bounds=Bounds(0, 1),
               options={"time_limit": time_limit})
    if res.x is None:
        print("Capacity variant: no solution found -", res.message,
              "(try a larger --eps)")
        return None
    if res.status != 0:
        print("Capacity variant: solution not proven optimal -", res.message)
    chosen = pairs[res.x > 0.5]
    parent = np.full(n, -1)
    parent[chosen[:, 1]] = chosen[:, 0]
    return parent


def show(name, parent, feeders, rho):
    """Print the edges (parent -> child) and the objective value."""
    print(f"\n=== {name} ===")
    non_root = np.where(parent >= 0)[0]
    for j in non_root:
        print(f"{feeders[parent[j]]} -> {feeders[j]}   (rho = {rho[parent[j], j]:.3f})")
    print(f"objective = {rho[parent[non_root], non_root].sum():.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="telemetry file: timestamp, feeder id, power (W)")
    ap.add_argument("--no-header", action="store_true", help="CSV has no header row")
    ap.add_argument("--round", default=None,
                    help="floor timestamps to this frequency, e.g. 10min")
    ap.add_argument("--eps", type=float, default=0.0,
                    help="tolerance in W for dominance and capacity (sensor accuracy)")
    ap.add_argument("--m-prime", type=int, default=200,
                    help="number of evenly spaced instants kept for capacity constraints")
    ap.add_argument("--time-limit", type=float, default=60.0, help="MILP time limit (s)")
    ap.add_argument("--out", default="topology.csv", help="output CSV")
    args = ap.parse_args()

    feeders, C = load_matrix(args.csv, not args.no_header, args.round)
    dom = candidate_parents(C, args.eps)
    rho = correlations(C)

    roots = np.where(~dom.any(axis=0))[0]
    if len(roots) != 1:
        raise SystemExit(f"Expected exactly one root, found {len(roots)}: "
                         f"{list(feeders[roots])}. The main feeder may be missing, "
                         f"or the data inconsistent (try a larger --eps).")
    print(f"root: {feeders[roots[0]]}")

    parent_base = solve_base(dom, rho)
    show("base", parent_base, feeders, rho)

    parent_cap = solve_capacity(C, dom, rho, args.m_prime, args.eps, args.time_limit)
    if parent_cap is not None:
        show("capacity", parent_cap, feeders, rho)
        diff = np.where(parent_base != parent_cap)[0]
        print(f"\nfeeders whose parent differs between variants: "
              f"{list(feeders[diff]) if len(diff) else 'none'}")

    out = pd.DataFrame({
        "feeder": feeders,
        "parent_base": [feeders[p] if p >= 0 else "" for p in parent_base],
        "parent_capacity": ([feeders[p] if p >= 0 else "" for p in parent_cap]
                            if parent_cap is not None else ""),
    })
    out.to_csv(args.out, index=False)
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()