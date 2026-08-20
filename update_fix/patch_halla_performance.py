#!/usr/bin/env python3

"""
Apply the MTD Explorer performance patch to HAllA 0.8.20.

The upstream HAllA 0.8.20 implementation evaluates Spearman correlations
through Python callbacks for both the X-by-Y association matrix and the
within-dataset hierarchical-clustering distances.  Large host-expression
matrices therefore spend most of their time repeatedly calling
scipy.stats.spearmanr one pair at a time.

For finite, non-constant continuous matrices using the Spearman metric, this
patch uses the mathematical identity

    Spearman(x, y) == Pearson(rank(x), rank(y))

and computes ranks once per feature.  The X-by-Y correlation matrix is then
computed with vectorized matrix multiplication, the corresponding asymptotic
Spearman p-values are computed from the same coefficients, and hierarchical
clustering uses SciPy's compiled correlation-distance implementation on the
ranked rows.

The original HAllA code is retained as an automatic fallback for other
metrics, forced permutation tests, missing values, constant features, or any
input that is not suitable for the fast path.

The patch also avoids materializing HAllA's unused square-form hierarchical
distance matrix.  The condensed distance vector used by scipy.cluster.linkage
is preserved.

The operation is idempotent, creates backups before changing upstream HAllA
modules, validates Python syntax, and runs numerical equivalence microtests.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import py_compile
import shutil
import sys
from importlib import metadata
from pathlib import Path

EXPECTED_HALLA_VERSION = "0.8.20"
PATCH_VERSION = "1"

MAIN_IMPORT_MARKER = "# MTD Explorer performance patch: fast Spearman helpers"
MAIN_FAST_MARKER = "# MTD Explorer performance patch: vectorized Spearman X-by-Y path"
HIERARCHY_IMPORT_MARKER = "# MTD Explorer performance patch: fast Spearman hierarchy helper"
HIERARCHY_FAST_MARKER = "# MTD Explorer performance patch: compiled Spearman hierarchy path"
HIERARCHY_MEMORY_MARKER = "# MTD Explorer performance patch: avoid unused square distance matrix"
HELPER_MARKER = f"# MTD Explorer HAllA fast Spearman helper v{PATCH_VERSION}"

HELPER_SOURCE = r'''# MTD Explorer HAllA fast Spearman helper v1
"""Fast, numerically equivalent Spearman helpers for HAllA 0.8.20.

This module is installed by update_fix/patch_halla_performance.py.
It intentionally activates only for finite, non-constant 2-D numeric matrices.
All unsupported cases fall back to upstream HAllA behavior in the patched
callers.
"""

import numpy as np
import scipy.spatial.distance as spd
from scipy.stats import distributions, rankdata


def _as_float_matrix(matrix):
    """Return a 2-D float64 NumPy view/copy suitable for numerical work."""
    data = np.asarray(matrix, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError("Expected a 2-D matrix.")
    return data


def can_use_fast_spearman(*matrices):
    """Return True only when the vectorized path is equivalent to HAllA."""
    if not matrices:
        return False

    sample_count = None
    for matrix in matrices:
        try:
            data = _as_float_matrix(matrix)
        except (TypeError, ValueError):
            return False

        if data.shape[0] == 0 or data.shape[1] < 3:
            return False

        if sample_count is None:
            sample_count = data.shape[1]
        elif data.shape[1] != sample_count:
            return False

        # Upstream HAllA removes missing values independently for each pair.
        # A single matrix-wide rank transform is not equivalent in that case,
        # so missing/non-finite data deliberately use the original code.
        if not np.isfinite(data).all():
            return False

        # Upstream HAllA returns (rho=0, p=1) for constant features.  Keeping
        # them on the original path avoids changing that special-case behavior.
        if np.any(np.ptp(data, axis=1) == 0):
            return False

    return True


def _rank_rows(matrix):
    """Rank each feature row with SciPy's average-tie convention."""
    data = _as_float_matrix(matrix)

    # SciPy 1.8 (the HAllA 0.8.20 environment used by MTD Explorer) supports
    # axis=1.  The fallback keeps the helper usable with older compatible
    # SciPy builds without changing the rank convention.
    try:
        return np.asarray(
            rankdata(data, method="average", axis=1),
            dtype=np.float64,
        )
    except TypeError:
        ranked = np.empty(data.shape, dtype=np.float64)
        for index in range(data.shape[0]):
            ranked[index, :] = rankdata(data[index, :], method="average")
        return ranked


def _center_and_normalize_rows(ranks):
    centered = ranks - ranks.mean(axis=1, keepdims=True)
    norms = np.sqrt(np.sum(centered * centered, axis=1))
    return centered, norms


def fast_spearman_similarity(X, Y):
    """Return the HAllA X-by-Y Spearman similarity matrix."""
    if not can_use_fast_spearman(X, Y):
        raise ValueError("Input is not eligible for the fast Spearman path.")

    x_ranks = _rank_rows(X)
    y_ranks = _rank_rows(Y)

    x_centered, x_norms = _center_and_normalize_rows(x_ranks)
    y_centered, y_norms = _center_and_normalize_rows(y_ranks)

    similarities = x_centered.dot(y_centered.T)
    similarities /= x_norms[:, None]
    similarities /= y_norms[None, :]

    # Mathematical correlations are bounded by [-1, 1].  Clipping only removes
    # possible last-bit floating-point excursions from vectorized arithmetic.
    np.clip(similarities, -1.0, 1.0, out=similarities)
    return similarities


def fast_spearman_pvalues(similarities, sample_count):
    """Return SciPy-compatible two-sided asymptotic Spearman p-values."""
    rho = np.asarray(similarities, dtype=np.float64)
    if sample_count < 3:
        raise ValueError("Spearman p-values require at least three samples.")

    dof = sample_count - 2
    clipped = np.clip(rho, -1.0, 1.0)

    # This is the same Student-t transformation used by scipy.stats.spearmanr
    # for the two-sided asymptotic p-value in the pinned SciPy generation.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = dof / ((clipped + 1.0) * (1.0 - clipped))
        np.clip(ratio, 0.0, None, out=ratio)
        t_stat = clipped * np.sqrt(ratio)

    pvalues = 2.0 * distributions.t.sf(np.abs(t_stat), dof)

    # Eligible fast-path inputs are non-constant and finite.  This guard is only
    # defensive against an unexpected numerical NaN.
    if np.isnan(pvalues).any():
        pvalues = np.where(np.isnan(pvalues), 1.0, pvalues)

    return pvalues


def fast_spearman_condensed_distance(
    matrix,
    sim2dist_set_abs=True,
    sim2dist_func=None,
):
    """Return HAllA's condensed 1-|Spearman| hierarchy distance vector.

    Memory is kept close to the unavoidable condensed O(n^2) vector.  The
    vector returned by scipy.spatial.distance.pdist is transformed in place
    whenever HAllA is using its default similarity-to-distance conversion.
    """
    if not can_use_fast_spearman(matrix):
        raise ValueError("Input is not eligible for the fast Spearman path.")

    ranks = _rank_rows(matrix)

    # SciPy's compiled correlation metric returns 1 - Pearson.  Applied to
    # ranks, Pearson is Spearman, so convert the condensed vector back to rho.
    condensed = spd.pdist(ranks, metric="correlation")
    np.subtract(1.0, condensed, out=condensed)
    np.clip(condensed, -1.0, 1.0, out=condensed)

    if sim2dist_func is not None:
        converted = np.asarray(sim2dist_func(condensed), dtype=np.float64)
        if converted.shape != condensed.shape:
            raise ValueError(
                "Custom similarity-to-distance function changed array shape."
            )
        return converted

    if sim2dist_set_abs:
        np.abs(condensed, out=condensed)

    np.subtract(1.0, condensed, out=condensed)
    return condensed
'''


def fail(message: str) -> None:
    print(f"[ERROR] {message}", file=sys.stderr)
    raise SystemExit(1)


def locate_halla_package() -> tuple[Path, str]:
    try:
        installed_version = metadata.version("halla")
    except metadata.PackageNotFoundError:
        fail("The halla Python package is not installed.")

    if installed_version != EXPECTED_HALLA_VERSION:
        fail(
            "This performance patch is intended for HAllA "
            f"{EXPECTED_HALLA_VERSION}, but version {installed_version} was found."
        )

    spec = importlib.util.find_spec("halla")
    if spec is None or not spec.submodule_search_locations:
        fail("Could not locate the installed HAllA package directory.")

    package_dir = Path(next(iter(spec.submodule_search_locations))).resolve()
    if not package_dir.is_dir():
        fail(f"HAllA package directory was not found: {package_dir}")

    return package_dir, installed_version


def compile_module(path: Path) -> None:
    try:
        py_compile.compile(str(path), doraise=True)
    except py_compile.PyCompileError as error:
        fail(f"Python syntax validation failed for {path}: {error}")


def backup_once(path: Path) -> None:
    backup = path.with_name(path.name + ".mtd_performance_original")
    if not backup.exists():
        shutil.copy2(path, backup)
        print(f"[INFO] Original HAllA module backed up to: {backup}")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        fail(
            f"Expected exactly one unpatched {label} block, but found {count}. "
            "No changes were written."
        )
    return text.replace(old, new, 1)


def patch_main(main_file: Path) -> None:
    text = main_file.read_text(encoding="utf-8")

    if MAIN_IMPORT_MARKER in text and MAIN_FAST_MARKER in text:
        print("[OK] HAllA main.py performance patch is already installed.")
        return

    if MAIN_IMPORT_MARKER in text or MAIN_FAST_MARKER in text:
        fail("A partial HAllA main.py performance patch was detected.")

    old_import = (
        "from .utils.stats import get_pvalue_table, pvalues2qvalues, test_pvalue_run_time\n"
    )
    new_import = old_import + (
        f"{MAIN_IMPORT_MARKER}\n"
        "from .utils.mtd_fast_spearman import (\n"
        "    can_use_fast_spearman,\n"
        "    fast_spearman_similarity,\n"
        "    fast_spearman_pvalues,\n"
        ")\n"
    )
    text = replace_once(text, old_import, new_import, "main.py import")

    start_anchor = "        # obtain similarity matrix\n"
    end_anchor = "        # obtain q-values\n"
    start = text.find(start_anchor)
    end = text.find(end_anchor, start + len(start_anchor))
    if start < 0 or end < 0:
        fail("Could not locate the HAllA pairwise similarity/p-value block.")

    original_block = text[start:end]
    required_fragments = (
        "self.similarity_table = spd.cdist(X, Y, metric=get_similarity_function(dist_metric))",
        "test_pvalue_run_time(X, Y, pdist_metric=dist_metric",
        "self.pvalue_table = get_pvalue_table(X, Y, pdist_metric=dist_metric",
    )
    for fragment in required_fragments:
        if fragment not in original_block:
            fail(
                "The HAllA pairwise block does not match the expected 0.8.20 "
                f"implementation; missing: {fragment}"
            )

    new_block = f'''        {MAIN_FAST_MARKER}\n        fast_spearman = (\n            dist_metric == 'spearman'\n            and not self.force_permutations\n            and can_use_fast_spearman(X, Y)\n        )\n\n        # obtain similarity matrix\n        self.logger.log_message('Generating the similarity table...')\n        if fast_spearman:\n            self.logger.log_message(\n                '[MTD Explorer] Using optimized vectorized Spearman similarity.'\n            )\n            self.similarity_table = fast_spearman_similarity(X, Y)\n        else:\n            self.similarity_table = spd.cdist(\n                X, Y, metric=get_similarity_function(dist_metric)\n            )\n\n        # obtain p-values\n        self.logger.log_message('Generating the p-value table...')\n        confp = config.permute\n        if fast_spearman:\n            self.logger.log_message(\n                '[MTD Explorer] Reusing Spearman coefficients for vectorized p-values.'\n            )\n            self.pvalue_table = fast_spearman_pvalues(\n                self.similarity_table, X.shape[1]\n            )\n        else:\n            extrapolated_time, timing_message = test_pvalue_run_time(\n                X, Y, pdist_metric=dist_metric,\n                permute_func=confp['func'], permute_iters=confp['iters'],\n                permute_speedup=confp['speedup'],\n                alpha=config.stats['fdr_alpha'],\n                force_perms=self.force_permutations,\n                num_threads=self.num_threads,\n                seed=self.seed\n            )\n            if extrapolated_time > 10 and self.verbose:\n                self.logger.log_message(timing_message)\n\n            self.pvalue_table = get_pvalue_table(\n                X, Y, pdist_metric=dist_metric,\n                permute_func=confp['func'], permute_iters=confp['iters'],\n                permute_speedup=confp['speedup'],\n                alpha=config.stats['fdr_alpha'],\n                no_progress=self.no_progress,\n                force_permutations=self.force_permutations,\n                num_threads=self.num_threads,\n                seed=self.seed\n            )\n\n'''

    text = text[:start] + new_block + text[end:]

    backup_once(main_file)
    main_file.write_text(text, encoding="utf-8")
    print(f"[INFO] Patched HAllA main module: {main_file}")


def patch_hierarchy(hierarchy_file: Path) -> None:
    text = hierarchy_file.read_text(encoding="utf-8")

    if (
        HIERARCHY_IMPORT_MARKER in text
        and HIERARCHY_FAST_MARKER in text
        and HIERARCHY_MEMORY_MARKER in text
    ):
        print("[OK] HAllA hierarchy.py performance patch is already installed.")
        return

    if any(
        marker in text
        for marker in (
            HIERARCHY_IMPORT_MARKER,
            HIERARCHY_FAST_MARKER,
            HIERARCHY_MEMORY_MARKER,
        )
    ):
        fail("A partial HAllA hierarchy.py performance patch was detected.")

    old_import = "from .utils.similarity import get_similarity_function, similarity2distance\n"
    new_import = old_import + (
        f"{HIERARCHY_IMPORT_MARKER}\n"
        "from .utils.mtd_fast_spearman import (\n"
        "    can_use_fast_spearman,\n"
        "    fast_spearman_condensed_distance,\n"
        ")\n"
    )
    text = replace_once(text, old_import, new_import, "hierarchy.py import")

    old_block = '''        self.distance_matrix = similarity2distance(spd.pdist(matrix, metric=treemetric),\n                                                   sim2dist_set_abs,\n                                                   sim2dist_func)\n        self.distance_matrix = np.clip(self.distance_matrix, a_min=0, a_max=None)\n        self.distance_matrix_sqr = spd.squareform(self.distance_matrix)\n        self._generate_hierarchical_clusters(linkage_method)\n'''

    new_block = f'''        {HIERARCHY_FAST_MARKER}\n        if pdist_metric == 'spearman' and can_use_fast_spearman(matrix):\n            self.distance_matrix = fast_spearman_condensed_distance(\n                matrix,\n                sim2dist_set_abs=sim2dist_set_abs,\n                sim2dist_func=sim2dist_func,\n            )\n        else:\n            self.distance_matrix = similarity2distance(\n                spd.pdist(matrix, metric=treemetric),\n                sim2dist_set_abs,\n                sim2dist_func\n            )\n\n        # Match the upstream non-negative clipping without allocating another\n        # O(n^2) condensed array.\n        np.maximum(self.distance_matrix, 0, out=self.distance_matrix)\n\n        {HIERARCHY_MEMORY_MARKER}\n        # distance_matrix_sqr is never consumed by HAllA 0.8.20.  Preserve the\n        # attribute for compatibility without allocating the square O(n^2) copy.\n        self.distance_matrix_sqr = None\n        self._generate_hierarchical_clusters(linkage_method)\n'''

    text = replace_once(text, old_block, new_block, "hierarchical distance")

    backup_once(hierarchy_file)
    hierarchy_file.write_text(text, encoding="utf-8")
    print(f"[INFO] Patched HAllA hierarchy module: {hierarchy_file}")


def write_helper(helper_file: Path) -> None:
    if helper_file.exists():
        current = helper_file.read_text(encoding="utf-8")
        if current == HELPER_SOURCE:
            print("[OK] HAllA fast Spearman helper is already installed.")
            return
        fail(
            "An unexpected existing mtd_fast_spearman.py was found. "
            "Refusing to overwrite it."
        )

    helper_file.write_text(HELPER_SOURCE, encoding="utf-8")
    print(f"[INFO] Installed HAllA fast Spearman helper: {helper_file}")


def load_helper_for_test(helper_file: Path):
    spec = importlib.util.spec_from_file_location(
        "mtd_halla_fast_spearman_validation",
        str(helper_file),
    )
    if spec is None or spec.loader is None:
        fail(f"Could not load helper module for validation: {helper_file}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_numerical_microtests(helper_file: Path) -> None:
    try:
        import numpy as np
        import scipy.spatial.distance as spd
        from scipy.stats import spearmanr
    except Exception as error:
        fail(f"Could not import NumPy/SciPy for patch validation: {error}")

    helper = load_helper_for_test(helper_file)

    # Deterministic data deliberately contain ties, which are common in the
    # microbiome input and exercise SciPy's average-rank behavior.
    X = np.array(
        [
            [1, 1, 2, 3, 5, 8, 13, 13],
            [8, 7, 7, 5, 4, 3, 2, 1],
            [0, 2, 1, 4, 3, 7, 6, 5],
            [2, 4, 4, 4, 8, 10, 12, 14],
        ],
        dtype=float,
    )
    Y = np.array(
        [
            [3, 3, 4, 5, 9, 9, 10, 12],
            [9, 8, 6, 6, 5, 3, 2, 0],
            [1, 3, 2, 5, 4, 8, 7, 6],
        ],
        dtype=float,
    )

    if not helper.can_use_fast_spearman(X, Y):
        fail("Fast Spearman helper unexpectedly rejected valid microtest data.")

    fast_rho = helper.fast_spearman_similarity(X, Y)
    fast_p = helper.fast_spearman_pvalues(fast_rho, X.shape[1])

    ref_rho = np.empty(fast_rho.shape, dtype=float)
    ref_p = np.empty(fast_p.shape, dtype=float)
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            ref_rho[i, j], ref_p[i, j] = spearmanr(X[i, :], Y[j, :])

    rho_diff = float(np.max(np.abs(fast_rho - ref_rho)))
    p_diff = float(np.max(np.abs(fast_p - ref_p)))

    reference_condensed = 1.0 - np.abs(
        spd.pdist(X, metric=lambda a, b: spearmanr(a, b)[0])
    )
    fast_condensed = helper.fast_spearman_condensed_distance(X)
    distance_diff = float(
        np.max(np.abs(fast_condensed - reference_condensed))
    )

    tolerance = 1e-12
    if rho_diff > tolerance:
        fail(
            "Fast Spearman similarity failed equivalence microtest: "
            f"max abs diff={rho_diff:.3e}"
        )
    if p_diff > tolerance:
        fail(
            "Fast Spearman p-value failed equivalence microtest: "
            f"max abs diff={p_diff:.3e}"
        )
    if distance_diff > tolerance:
        fail(
            "Fast Spearman hierarchy distance failed equivalence microtest: "
            f"max abs diff={distance_diff:.3e}"
        )

    invalid_nan = X.copy()
    invalid_nan[0, 0] = np.nan
    invalid_constant = X.copy()
    invalid_constant[0, :] = 1.0

    if helper.can_use_fast_spearman(invalid_nan):
        fail("Fast Spearman helper did not reject missing values.")
    if helper.can_use_fast_spearman(invalid_constant):
        fail("Fast Spearman helper did not reject a constant feature.")

    print("[OK] Fast Spearman numerical equivalence microtests passed.")
    print(f"[INFO] Max similarity difference: {rho_diff:.3e}")
    print(f"[INFO] Max p-value difference:    {p_diff:.3e}")
    print(f"[INFO] Max distance difference:   {distance_diff:.3e}")


def validate_patch(package_dir: Path) -> None:
    main_file = package_dir / "main.py"
    hierarchy_file = package_dir / "hierarchy.py"
    helper_file = package_dir / "utils" / "mtd_fast_spearman.py"

    for path in (main_file, hierarchy_file, helper_file):
        if not path.is_file():
            fail(f"Required patched HAllA file is missing: {path}")

    main_text = main_file.read_text(encoding="utf-8")
    hierarchy_text = hierarchy_file.read_text(encoding="utf-8")
    helper_text = helper_file.read_text(encoding="utf-8")

    for marker in (MAIN_IMPORT_MARKER, MAIN_FAST_MARKER):
        if main_text.count(marker) != 1:
            fail(f"Expected exactly one marker in main.py: {marker}")

    for marker in (
        HIERARCHY_IMPORT_MARKER,
        HIERARCHY_FAST_MARKER,
        HIERARCHY_MEMORY_MARKER,
    ):
        if hierarchy_text.count(marker) != 1:
            fail(f"Expected exactly one marker in hierarchy.py: {marker}")

    if HELPER_MARKER not in helper_text:
        fail("HAllA fast Spearman helper marker is missing.")

    expected_sha = hashlib.sha256(HELPER_SOURCE.encode("utf-8")).hexdigest()
    observed_sha = hashlib.sha256(helper_text.encode("utf-8")).hexdigest()
    if observed_sha != expected_sha:
        fail("HAllA fast Spearman helper contents do not match this patch version.")

    for path in (main_file, hierarchy_file, helper_file):
        compile_module(path)

    run_numerical_microtests(helper_file)

    print("[OK] HAllA 0.8.20 performance patch validation passed.")
    print(f"[INFO] HAllA package: {package_dir}")


def apply_patch(package_dir: Path) -> None:
    main_file = package_dir / "main.py"
    hierarchy_file = package_dir / "hierarchy.py"
    helper_file = package_dir / "utils" / "mtd_fast_spearman.py"

    for path in (main_file, hierarchy_file):
        if not path.is_file():
            fail(f"Expected HAllA 0.8.20 module was not found: {path}")

    # Write the helper first; main.py/hierarchy.py will not import it until a
    # future HAllA process starts, and validation only occurs after both callers
    # have been patched successfully.
    write_helper(helper_file)
    patch_main(main_file)
    patch_hierarchy(hierarchy_file)
    validate_patch(package_dir)

    print("[OK] HAllA 0.8.20 performance patch applied.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Apply or validate the MTD Explorer HAllA 0.8.20 "
            "Spearman performance patch."
        )
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate the patch without modifying HAllA.",
    )
    args = parser.parse_args()

    package_dir, installed_version = locate_halla_package()
    print(f"[INFO] HAllA version: {installed_version}")
    print(f"[INFO] HAllA package directory: {package_dir}")

    if args.check:
        validate_patch(package_dir)
    else:
        apply_patch(package_dir)


if __name__ == "__main__":
    main()
