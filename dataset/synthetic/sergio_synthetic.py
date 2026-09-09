#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate two SERGIO spliced/unspliced synthetic scRNA-seq datasets and save as Loom.

Outputs
-------
1) sergio_5class_600genes.loom
   - 5 classes
   - 600 genes
   - cells/class = [200, 400, 600, 800, 1000]
   - total cells = 3000
   - exactly 100 designed marker genes/class

2) sergio_10class_1100genes.loom
   - 10 classes
   - 1100 genes
   - cells/class = [200, 200, 300, 300, 500, 500, 800, 1000, 1000, 1200]
   - total cells = 6000
   - exactly 100 designed marker genes/class

Loom orientation: gene x cell
Downstream access:
    S = ds.layers['spliced'][:,:]
    U = ds.layers['unspliced'][:,:]
    true_labels = ds.ca['subclass_label']
    # barcode = ds.col_attrs['barcode']
    g_names = ds.ra['gene_name']
"""

from pathlib import Path
import csv
import sys
import tempfile
import numpy as np
import loompy

# SERGIO source uses deprecated NumPy aliases.
if not hasattr(np, "int"):
    np.int = int
if not hasattr(np, "float"):
    np.float = float


def load_sergio():
    """Find a cloned SERGIO repository and import SERGIO.sergio.sergio."""
    here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    candidates = [Path.cwd(), here, here / "SERGIO", Path.cwd() / "SERGIO", here.parent / "SERGIO"]
    for root in candidates:
        if (root / "SERGIO" / "sergio.py").exists():
            sys.path.insert(0, str(root))
            break
    try:
        from SERGIO.sergio import sergio
    except ImportError as exc:
        raise ImportError(
            "Cannot import SERGIO. First run:\n"
            "git clone https://github.com/PayamDiba/SERGIO.git\n"
            "Then place this script inside the cloned repo or next to the cloned SERGIO folder."
        ) from exc
    return sergio


sergio = load_sergio()

# ============================================================
# Settings
# ============================================================
SEED = 2026
OUTPUT_DIR = Path("sergio_loom_output")

MARKERS_PER_CLASS = 100
N_SECONDARY = 15                    # 1 MR + 15 secondary + 84 downstream = 100
BACKGROUND_GENES = 100

# SERGIO biological/stochastic settings
NOISE_U = 0.25
NOISE_S = 0.20
NOISE_TYPE = "dpd"
DECAY_RANGE = (0.80, 1.00)
SPLICE_RATIO_RANGE = (3.0, 4.2)
DT = 0.01
SAMPLING_STATE = 1
HILL_COOP_STATE = 2

# Class-specific GRN strength
MR_HIGH = (6.0, 7.2)
MR_LOW = (0.30, 0.55)
K_SECONDARY = (5.5, 7.0)
K_DOWNSTREAM = (2.5, 3.8)
K_BACKGROUND = (0.25, 0.45)

# Technical noise / medium-high count regime
LIBRARY_SIGMA = 0.18
DROPOUT_SHAPE = 1.0
DROPOUT_PERCENTILE = 20
MAX_EXPECTED_TOTAL_PER_ENTRY = 500.0
MAX_OBSERVED_TOTAL_PER_ENTRY = 800

DATASETS = [
    dict(
        name="sergio_5class_600genes",
        n_classes=5,
        n_genes=600,
        cells_per_class=[200, 400, 600, 800, 1000],
        target_library_size=20000,
        seed=SEED,
    ),
    dict(
        name="sergio_10class_1100genes",
        n_classes=10,
        n_genes=1100,
        cells_per_class=[200, 200, 300, 300, 500, 500, 800, 1000, 1000, 1200],
        target_library_size=35000,
        seed=SEED + 1,
    ),
]


# ============================================================
# GRN design
# ============================================================
def make_design(n_classes, n_genes):
    expected = n_classes * MARKERS_PER_CLASS + BACKGROUND_GENES
    if n_genes != expected:
        raise ValueError(f"Expected {expected} genes, got {n_genes}.")

    # One class-specific master regulator per class.
    master = list(range(n_classes))
    secondary, downstream, markers = {}, {}, {}

    gene_name = np.empty(n_genes, dtype=object)
    marker_class = np.full(n_genes, -1, dtype=np.int32)
    gene_role = np.empty(n_genes, dtype=object)

    for c, g in enumerate(master):
        gene_name[g] = f"C{c}_MR"
        marker_class[g] = c
        gene_role[g] = "master_regulator"

    cursor = n_classes
    n_downstream = MARKERS_PER_CLASS - 1 - N_SECONDARY  # 84

    for c in range(n_classes):
        secondary[c] = list(range(cursor, cursor + N_SECONDARY))
        cursor += N_SECONDARY
        downstream[c] = list(range(cursor, cursor + n_downstream))
        cursor += n_downstream
        markers[c] = [master[c]] + secondary[c] + downstream[c]

        for j, g in enumerate(secondary[c], 1):
            gene_name[g] = f"C{c}_SEC_{j:02d}"
            marker_class[g] = c
            gene_role[g] = "secondary_regulator"
        for j, g in enumerate(downstream[c], 1):
            gene_name[g] = f"C{c}_MARKER_{j:03d}"
            marker_class[g] = c
            gene_role[g] = "downstream_marker"

    background = list(range(cursor, n_genes))
    if len(background) != BACKGROUND_GENES:
        raise RuntimeError("Background gene number is incorrect.")
    for j, g in enumerate(background, 1):
        gene_name[g] = f"BG_{j:03d}"
        gene_role[g] = "background"

    for c in range(n_classes):
        assert len(markers[c]) == 100
        assert np.sum(marker_class == c) == 100

    return dict(
        master=master,
        secondary=secondary,
        downstream=downstream,
        markers=markers,
        background=background,
        gene_name=gene_name,
        marker_class=marker_class,
        gene_role=gene_role,
    )


def write_grn_files(folder, design, n_classes, rng):
    """
    With shared_coop_state > 0, SERGIO target rows are:
      target, n_reg, reg1,...,regN, K1,...,KN
    Master-regulator rows are:
      MR, rate_class0,...,rate_classN
    """
    folder = Path(folder)
    target_file = folder / "targets.csv"
    mr_file = folder / "master_regulators.csv"

    # Class-specific MR basal production.
    with mr_file.open("w", newline="") as f:
        w = csv.writer(f)
        for c, mr in enumerate(design["master"]):
            rates = rng.uniform(*MR_LOW, size=n_classes)
            rates[c] = rng.uniform(*MR_HIGH)
            w.writerow([mr] + [f"{x:.8f}" for x in rates])

    with target_file.open("w", newline="") as f:
        w = csv.writer(f)
        for c in range(n_classes):
            mr = design["master"][c]

            # MR -> 15 secondary marker regulators
            for target in design["secondary"][c]:
                w.writerow([target, 1, mr, f"{rng.uniform(*K_SECONDARY):.8f}"])

            # Two secondary regulators -> 84 downstream markers
            sec = np.asarray(design["secondary"][c], dtype=int)
            for target in design["downstream"][c]:
                r1, r2 = rng.choice(sec, 2, replace=False)
                k1, k2 = rng.uniform(*K_DOWNSTREAM, size=2)
                w.writerow([target, 2, int(r1), int(r2), f"{k1:.8f}", f"{k2:.8f}"])

        # 100 shared/background genes: weakly regulated by all class MRs.
        for target in design["background"]:
            ks = rng.uniform(*K_BACKGROUND, size=n_classes)
            w.writerow(
                [target, n_classes]
                + design["master"]
                + [f"{x:.8f}" for x in ks]
            )

    return target_file, mr_file


# ============================================================
# Technical noise and UMI count generation
# ============================================================
def to_umi_counts(sim, U_expr, S_expr, target_library_size):
    U_expr = np.maximum(np.asarray(U_expr, dtype=np.float64), 0.0)
    S_expr = np.maximum(np.asarray(S_expr, dtype=np.float64), 0.0)

    # Set mean U+S library size using SERGIO's dynamics-specific library-size effect.
    sigma = LIBRARY_SIGMA
    mu = np.log(float(target_library_size)) - 0.5 * sigma**2
    _, U_lam, S_lam = sim.lib_size_effect_dynamics(
        U_expr, S_expr, mean=mu, scale=sigma
    )

    # Moderate dropout, separately on unspliced and spliced.
    keep_U, keep_S = sim.dropout_indicator_dynamics(
        U_lam, S_lam,
        shape=DROPOUT_SHAPE,
        percentile=DROPOUT_PERCENTILE,
    )
    U_lam *= keep_U
    S_lam *= keep_S

    # One global scaling preserves relative U/S and class/gene structure.
    peak_lambda = float(np.max(U_lam + S_lam))
    if peak_lambda > MAX_EXPECTED_TOTAL_PER_ENTRY:
        factor = MAX_EXPECTED_TOTAL_PER_ENTRY / peak_lambda
        U_lam *= factor
        S_lam *= factor

    # Poisson UMI sampling. No element-wise clipping is used.
    # If an extremely rare draw exceeds 800, globally shrink and resample.
    for _ in range(20):
        U, S = sim.convert_to_UMIcounts_dynamics(U_lam, S_lam)
        peak = int(np.max(U + S))
        if peak <= MAX_OBSERVED_TOTAL_PER_ENTRY:
            return U.astype(np.int32), S.astype(np.int32)
        factor = min(0.90, 0.90 * MAX_OBSERVED_TOTAL_PER_ENTRY / peak)
        U_lam *= factor
        S_lam *= factor

    raise RuntimeError("Failed to enforce max(S+U) <= 800.")


# ============================================================
# Loom writer
# ============================================================
def save_as_loom(path, U, S, labels, barcodes, design):
    """Save all matrices in gene x cell orientation."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    U = np.asarray(U, dtype=np.int32)
    S = np.asarray(S, dtype=np.int32)
    X = U + S

    row_attrs = {
        "gene_name": np.asarray(design["gene_name"], dtype=str),
        "marker_class": np.asarray(design["marker_class"], dtype=np.int32),
        "gene_role": np.asarray(design["gene_role"], dtype=str),
    }
    col_attrs = {
        "subclass_label": np.asarray(labels, dtype=np.int32),
        "barcode": np.asarray(barcodes, dtype=str),
    }

    # Main matrix is total counts; additional layers are S and U.
    loompy.create(str(path), X, row_attrs=row_attrs, col_attrs=col_attrs)
    with loompy.connect(str(path), mode="r+") as ds:
        ds.layers["spliced"] = S
        ds.layers["unspliced"] = U

    # Exact downstream-format validation.
    with loompy.connect(str(path), mode="r") as ds:
        S2 = ds.layers["spliced"][:, :]
        U2 = ds.layers["unspliced"][:, :]
        y2 = ds.ca["subclass_label"]
        g2 = ds.ra["gene_name"]
        assert S2.shape == S.shape
        assert U2.shape == U.shape
        assert len(y2) == U.shape[1]
        assert len(g2) == U.shape[0]
        assert int(np.max(S2 + U2)) <= 800


# ============================================================
# One complete dataset
# ============================================================
def simulate_one(cfg):
    name = cfg["name"]
    n_classes = cfg["n_classes"]
    n_genes = cfg["n_genes"]
    cells_per_class = np.asarray(cfg["cells_per_class"], dtype=np.int32)
    seed = cfg["seed"]

    if len(cells_per_class) != n_classes:
        raise ValueError(
            f"{name}: cells_per_class has {len(cells_per_class)} entries, "
            f"but n_classes={n_classes}."
        )
    if np.any(cells_per_class <= 0):
        raise ValueError(f"{name}: every class must contain at least one cell.")

    # SERGIO's number_sc is a single value shared across bins/classes.
    # Simulate each class at the largest requested class size and then
    # independently subsample each class to the requested size.
    n_sc_sergio = int(cells_per_class.max())
    total_cells = int(cells_per_class.sum())

    print("\n" + "=" * 72)
    print(name)
    print(f"classes={n_classes}, genes={n_genes}")
    print(f"requested cells/class={cells_per_class.tolist()}")
    print(f"SERGIO cells/bin={n_sc_sergio}; final total cells={total_cells}")
    print("markers/class=100")
    print("=" * 72)

    # SERGIO itself uses global np.random.
    np.random.seed(seed)
    rng = np.random.default_rng(seed)

    design = make_design(n_classes, n_genes)
    decays = rng.uniform(*DECAY_RANGE, size=n_genes)
    splice_ratio = rng.uniform(*SPLICE_RATIO_RANGE, size=n_genes)

    # Classes are independent SERGIO bins for clustering benchmarking.
    # A zero bifurcation matrix still uses SERGIO's splicing dynamics.
    bifurcation = np.zeros((n_classes, n_classes), dtype=float)

    with tempfile.TemporaryDirectory(prefix=f"{name}_") as tmp:
        target_file, mr_file = write_grn_files(tmp, design, n_classes, rng)

        sim = sergio(
            number_genes=n_genes,
            number_bins=n_classes,
            number_sc=n_sc_sergio,
            noise_params=NOISE_U,
            noise_type=NOISE_TYPE,
            decays=decays,
            dynamics=True,
            sampling_state=SAMPLING_STATE,
            dt=DT,
            bifurcation_matrix=bifurcation,
            noise_params_splice=NOISE_S,
            noise_type_splice=NOISE_TYPE,
            splice_ratio=splice_ratio,
            dt_splice=DT,
        )

        sim.build_graph(
            input_file_taregts=str(target_file),
            input_file_regs=str(mr_file),
            shared_coop_state=HILL_COOP_STATE,
        )

        sim.simulate_dynamics()

        # SERGIO output: class/bin x gene x cell
        U_expr, S_expr = sim.getExpressions_dynamics()
        U3, S3 = to_umi_counts(sim, U_expr, S_expr, cfg["target_library_size"])

    # SERGIO produces the same number of cells in every bin.
    # Keep a random subset from each bin so final class sizes match the
    # requested unequal cell counts exactly.
    U_blocks = []
    S_blocks = []
    label_blocks = []
    barcode_blocks = []

    for c, n_keep in enumerate(cells_per_class):
        n_keep = int(n_keep)
        if U3[c].shape[1] < n_keep or S3[c].shape[1] < n_keep:
            raise RuntimeError(
                f"Class {c}: SERGIO returned only {U3[c].shape[1]} cells, "
                f"cannot keep requested {n_keep}."
            )

        keep = rng.choice(U3[c].shape[1], size=n_keep, replace=False)
        U_blocks.append(U3[c][:, keep])
        S_blocks.append(S3[c][:, keep])
        label_blocks.append(np.full(n_keep, c, dtype=np.int32))
        barcode_blocks.extend(
            [f"class{c}_cell{int(i):04d}" for i in keep]
        )

    # gene x all_cells
    U = np.concatenate(U_blocks, axis=1)
    S = np.concatenate(S_blocks, axis=1)
    labels = np.concatenate(label_blocks)
    barcodes = np.asarray(barcode_blocks, dtype=object)

    # Shuffle cells so classes are not stored as consecutive blocks.
    order = rng.permutation(U.shape[1])
    U, S = U[:, order], S[:, order]
    labels, barcodes = labels[order], barcodes[order]
    X = U + S

    # Strict checks.
    expected_shape = (n_genes, total_cells)
    assert U.shape == expected_shape
    assert S.shape == expected_shape
    assert int(X.max()) <= 800
    for c, n_expected in enumerate(cells_per_class):
        assert np.sum(design["marker_class"] == c) == 100
        assert np.sum(labels == c) == int(n_expected)

    out = OUTPUT_DIR / f"{name}.loom"
    save_as_loom(out, U, S, labels, barcodes, design)

    # QC summary.
    lib = X.sum(axis=0)
    print(f"Saved              : {out.resolve()}")
    print(f"shape              : {X.shape}  (gene x cell)")
    print(f"max spliced        : {int(S.max())}")
    print(f"max unspliced      : {int(U.max())}")
    print(f"max spliced+unsp.  : {int(X.max())}")
    print(f"mean library       : {lib.mean():.2f}")
    print(f"median library     : {np.median(lib):.2f}")
    print(f"zero fraction      : {np.mean(X == 0):.4f}")

    # Check the designed markers are enriched in their own class.
    print("Marker QC:")
    for c in range(n_classes):
        genes = design["markers"][c]
        own = X[np.ix_(genes, labels == c)].mean(axis=1)
        other = X[np.ix_(genes, labels != c)].mean(axis=1)
        fc = (own + 1.0) / (other + 1.0)
        print(f"  class {c}: n_marker={len(genes)}, median_FC={np.median(fc):.3f}")

    return out


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    files = [simulate_one(cfg) for cfg in DATASETS]

    print("\nGenerated Loom files:")
    for x in files:
        print(x.resolve())

    print("\nRead with exactly:")
    print("import loompy")
    print("ds = loompy.connect('sergio_loom_output/sergio_5class_600genes.loom')")
    print("S = ds.layers['spliced'][:,:]")
    print("U = ds.layers['unspliced'][:,:]")
    print("true_labels = ds.ca['subclass_label']")
    print("# barcode = ds.col_attrs['barcode']")
    print("g_names = ds.ra['gene_name']")


if __name__ == "__main__":
    main()
