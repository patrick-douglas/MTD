#!/usr/bin/env python3

"""
Apply the MTD Explorer correlation-performance patch to HAllA 0.8.20.

HAllA 0.8.20 evaluates Pearson and Spearman correlations through Python
callbacks for the X-by-Y association matrix and for within-dataset
hierarchical-clustering distances.  On large host-expression matrices this
can require millions to hundreds of millions of Python-level scipy.stats
calls.

For finite, non-constant matrices this patch provides fast paths for:

* Spearman: rank each feature once, then use vectorized Pearson correlation
  on the ranks (Spearman == Pearson(rank(x), rank(y))).
* Pearson: center/normalize each feature once, then use vectorized matrix
  multiplication directly.
* P-values: reuse the already computed correlation coefficients instead of
  recomputing every pair with scipy.stats.spearmanr/pearsonr.
* Hierarchical clustering: use SciPy's compiled ``pdist(...,
  metric="correlation")`` implementation (on ranks for Spearman; raw values
  for Pearson) instead of a Python callback for every pair.

The original HAllA code remains the automatic fallback for other metrics,
forced permutation tests, missing/non-finite values, constant features, and
inputs that are not suitable for the fast paths.

The patch also avoids materializing HAllA's unused square-form hierarchical
distance matrix.  The condensed distance vector required by
scipy.cluster.linkage is preserved.

This script can safely upgrade the previous MTD Explorer Spearman-only patch
(v1) by restoring the original HAllA 0.8.20 main.py/hierarchy.py backups and
then applying this v2 patch.  The operation is idempotent, validates Python
syntax, and runs numerical equivalence microtests for both Spearman and
Pearson.
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
PATCH_VERSION = "2"

# v2 markers
MAIN_IMPORT_MARKER = "# MTD Explorer performance patch v2: fast correlation helpers"
MAIN_FAST_MARKER = "# MTD Explorer performance patch v2: vectorized Pearson/Spearman X-by-Y path"
HIERARCHY_IMPORT_MARKER = "# MTD Explorer performance patch v2: fast correlation hierarchy helpers"
HIERARCHY_FAST_MARKER = "# MTD Explorer performance patch v2: compiled Pearson/Spearman hierarchy path"
HIERARCHY_MEMORY_MARKER = "# MTD Explorer performance patch v2: avoid unused square distance matrix"
HELPER_MARKER = f"# MTD Explorer HAllA fast correlation helper v{PATCH_VERSION}"
HELPER_FILENAME = "mtd_fast_correlation.py"

# Previous Spearman-only patch markers.  These are used only for a controlled
# v1 -> v2 upgrade of an already patched halla0820 environment.
V1_MAIN_IMPORT_MARKER = "# MTD Explorer performance patch: fast Spearman helpers"
V1_MAIN_FAST_MARKER = "# MTD Explorer performance patch: vectorized Spearman X-by-Y path"
V1_HIERARCHY_IMPORT_MARKER = "# MTD Explorer performance patch: fast Spearman hierarchy helper"
V1_HIERARCHY_FAST_MARKER = "# MTD Explorer performance patch: compiled Spearman hierarchy path"
V1_HIERARCHY_MEMORY_MARKER = "# MTD Explorer performance patch: avoid unused square distance matrix"
V1_HELPER_FILENAME = "mtd_fast_spearman.py"
V1_HELPER_SHA256 = "99b1d67be4014af91353732f76b87ed213cba7c9acb8e3c3a760ab3043781c2b"

HELPER_SOURCE = r'''# MTD Explorer HAllA fast correlation helper v2
"""Fast Pearson/Spearman helpers for HAllA 0.8.20.

Installed by update_fix/patch_halla_performance.py.

The fast paths deliberately activate only for finite, non-constant 2-D
numeric matrices.  Unsupported cases remain on upstream HAllA behavior in the
patched callers.
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


def _basic_fast_path_eligibility(*matrices):
    """Check conditions shared by the Pearson and Spearman fast paths."""
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
        # Matrix-wide vectorization is not equivalent in that situation.
        if not np.isfinite(data).all():
            return False

        # Upstream HAllA returns (correlation=0, p=1) for constant features.
        if np.any(np.ptp(data, axis=1) == 0):
            return False

    return True


def can_use_fast_spearman(*matrices):
    """Return True when the vectorized Spearman path is safe."""
    return _basic_fast_path_eligibility(*matrices)


def can_use_fast_pearson(*matrices):
    """Return True when the vectorized Pearson path is safe.

    A conservative near-constant check keeps numerically delicate rows on
    scipy.stats.pearsonr, matching the upstream HAllA path rather than risking
    a change from cancellation during matrix-wide centering.
    """
    if not _basic_fast_path_eligibility(*matrices):
        return False

    for matrix in matrices:
        data = _as_float_matrix(matrix)
        means = data.mean(axis=1)
        centered_norms = np.linalg.norm(data - means[:, None], axis=1)

        # SciPy has historically warned/fallen into a numerically delicate
        # regime for nearly constant inputs.  The exact threshold has varied
        # across releases, so use a conservative guard and let upstream HAllA
        # handle such rows pair-by-pair.
        scale = np.maximum(np.abs(means), 1.0)
        if np.any(centered_norms <= 1e-12 * scale):
            return False

    return True


def _rank_rows(matrix):
    """Rank each feature row with SciPy's average-tie convention."""
    data = _as_float_matrix(matrix)
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


def _center_and_normalize_rows(data):
    centered = data - data.mean(axis=1, keepdims=True)
    norms = np.sqrt(np.sum(centered * centered, axis=1))
    return centered, norms


def _cross_correlation(X, Y):
    x_data = _as_float_matrix(X)
    y_data = _as_float_matrix(Y)

    x_centered, x_norms = _center_and_normalize_rows(x_data)
    y_centered, y_norms = _center_and_normalize_rows(y_data)

    similarities = x_centered.dot(y_centered.T)
    similarities /= x_norms[:, None]
    similarities /= y_norms[None, :]

    # Remove only last-bit excursions outside the mathematical bounds.
    np.clip(similarities, -1.0, 1.0, out=similarities)
    return similarities


def fast_pearson_similarity(X, Y):
    """Return the HAllA X-by-Y Pearson similarity matrix."""
    if not can_use_fast_pearson(X, Y):
        raise ValueError("Input is not eligible for the fast Pearson path.")
    return _cross_correlation(X, Y)


def fast_spearman_similarity(X, Y):
    """Return the HAllA X-by-Y Spearman similarity matrix."""
    if not can_use_fast_spearman(X, Y):
        raise ValueError("Input is not eligible for the fast Spearman path.")
    return _cross_correlation(_rank_rows(X), _rank_rows(Y))


def _fast_correlation_pvalues(similarities, sample_count):
    """Two-sided correlation p-values from already computed coefficients.

    For n >= 3, the usual Student-t transformation is mathematically
    equivalent to the null distribution used by scipy.stats.pearsonr and is
    also the transformation used for scipy.stats.spearmanr's asymptotic
    two-sided p-value in the HAllA 0.8.20 SciPy generation.
    """
    corr = np.asarray(similarities, dtype=np.float64)
    if sample_count < 3:
        raise ValueError("Correlation p-values require at least three samples.")

    dof = sample_count - 2
    clipped = np.clip(corr, -1.0, 1.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = dof / ((clipped + 1.0) * (1.0 - clipped))
        np.clip(ratio, 0.0, None, out=ratio)
        t_stat = clipped * np.sqrt(ratio)

    pvalues = 2.0 * distributions.t.sf(np.abs(t_stat), dof)

    if np.isnan(pvalues).any():
        pvalues = np.where(np.isnan(pvalues), 1.0, pvalues)

    return pvalues


def fast_pearson_pvalues(similarities, sample_count):
    return _fast_correlation_pvalues(similarities, sample_count)


def fast_spearman_pvalues(similarities, sample_count):
    return _fast_correlation_pvalues(similarities, sample_count)


def _condensed_correlation_distance(
    matrix,
    sim2dist_set_abs=True,
    sim2dist_func=None,
):
    """Convert compiled correlation distance to HAllA's condensed distance."""
    condensed = spd.pdist(matrix, metric="correlation")

    # scipy correlation distance is 1-r.  Recover r so HAllA's
    # similarity2distance semantics can be applied exactly.
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


def fast_pearson_condensed_distance(
    matrix,
    sim2dist_set_abs=True,
    sim2dist_func=None,
):
    """Return HAllA's condensed Pearson hierarchy distance vector."""
    if not can_use_fast_pearson(matrix):
        raise ValueError("Input is not eligible for the fast Pearson path.")
    return _condensed_correlation_distance(
        _as_float_matrix(matrix),
        sim2dist_set_abs=sim2dist_set_abs,
        sim2dist_func=sim2dist_func,
    )


def fast_spearman_condensed_distance(
    matrix,
    sim2dist_set_abs=True,
    sim2dist_func=None,
):
    """Return HAllA's condensed Spearman hierarchy distance vector."""
    if not can_use_fast_spearman(matrix):
        raise ValueError("Input is not eligible for the fast Spearman path.")
    return _condensed_correlation_distance(
        _rank_rows(matrix),
        sim2dist_set_abs=sim2dist_set_abs,
        sim2dist_func=sim2dist_func,
    )
'''


def fail(message: str) -> None:
    print(f"[ERROR] {message}", file=sys.stderr)
    raise SystemExit(1)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def _verify_upstream_backup(path: Path, kind: str) -> None:
    if not path.is_file():
        fail(f"Required original HAllA backup is missing: {path}")

    text = path.read_text(encoding="utf-8")
    if "MTD Explorer performance patch" in text:
        fail(f"Original HAllA backup unexpectedly contains an MTD patch: {path}")

    if kind == "main":
        required = (
            "self.similarity_table = spd.cdist(X, Y, metric=get_similarity_function(dist_metric))",
            "self.pvalue_table = get_pvalue_table(X, Y, pdist_metric=dist_metric",
        )
    else:
        required = (
            "self.distance_matrix = similarity2distance(spd.pdist(matrix, metric=treemetric)",
            "self.distance_matrix_sqr = spd.squareform(self.distance_matrix)",
        )

    for fragment in required:
        if fragment not in text:
            fail(
                f"Original HAllA backup does not match expected 0.8.20 {kind}.py: "
                f"missing {fragment!r}"
            )


def upgrade_v1_if_needed(package_dir: Path) -> None:
    """Restore pristine HAllA modules before upgrading the MTD v1 patch."""
    main_file = package_dir / "main.py"
    hierarchy_file = package_dir / "hierarchy.py"
    old_helper = package_dir / "utils" / V1_HELPER_FILENAME

    main_text = main_file.read_text(encoding="utf-8")
    hierarchy_text = hierarchy_file.read_text(encoding="utf-8")

    v1_main_markers = (V1_MAIN_IMPORT_MARKER, V1_MAIN_FAST_MARKER)
    v1_hierarchy_markers = (
        V1_HIERARCHY_IMPORT_MARKER,
        V1_HIERARCHY_FAST_MARKER,
        V1_HIERARCHY_MEMORY_MARKER,
    )

    present_main = [marker in main_text for marker in v1_main_markers]
    present_hierarchy = [marker in hierarchy_text for marker in v1_hierarchy_markers]
    any_v1 = any(present_main) or any(present_hierarchy) or old_helper.exists()

    if not any_v1:
        return

    if not all(present_main) or not all(present_hierarchy):
        fail("A partial previous MTD HAllA performance patch was detected.")

    if not old_helper.is_file():
        fail(f"Previous MTD HAllA helper is missing: {old_helper}")

    if sha256_file(old_helper) != V1_HELPER_SHA256:
        fail(
            "The existing Spearman-only helper does not match the known v1 "
            "MTD patch. Refusing an automatic upgrade."
        )

    main_backup = main_file.with_name(main_file.name + ".mtd_performance_original")
    hierarchy_backup = hierarchy_file.with_name(
        hierarchy_file.name + ".mtd_performance_original"
    )

    _verify_upstream_backup(main_backup, "main")
    _verify_upstream_backup(hierarchy_backup, "hierarchy")

    shutil.copy2(main_backup, main_file)
    shutil.copy2(hierarchy_backup, hierarchy_file)
    old_helper.unlink()

    print("[INFO] Previous Spearman-only HAllA performance patch v1 detected.")
    print("[INFO] Restored pristine HAllA 0.8.20 modules from validated backups.")
    print("[INFO] Removed the v1 helper before applying performance patch v2.")


def patch_main(main_file: Path) -> None:
    text = main_file.read_text(encoding="utf-8")

    if MAIN_IMPORT_MARKER in text and MAIN_FAST_MARKER in text:
        print("[OK] HAllA main.py performance patch v2 is already installed.")
        return

    if MAIN_IMPORT_MARKER in text or MAIN_FAST_MARKER in text:
        fail("A partial HAllA main.py performance patch v2 was detected.")

    old_import = (
        "from .utils.stats import get_pvalue_table, pvalues2qvalues, test_pvalue_run_time\n"
    )
    new_import = old_import + (
        f"{MAIN_IMPORT_MARKER}\n"
        "from .utils.mtd_fast_correlation import (\n"
        "    can_use_fast_pearson,\n"
        "    can_use_fast_spearman,\n"
        "    fast_pearson_similarity,\n"
        "    fast_pearson_pvalues,\n"
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

    new_block = f'''        {MAIN_FAST_MARKER}\n        fast_spearman = (\n            dist_metric == 'spearman'\n            and not self.force_permutations\n            and can_use_fast_spearman(X, Y)\n        )\n        fast_pearson = (\n            dist_metric == 'pearson'\n            and not self.force_permutations\n            and can_use_fast_pearson(X, Y)\n        )\n\n        # obtain similarity matrix\n        self.logger.log_message('Generating the similarity table...')\n        if fast_spearman:\n            self.logger.log_message(\n                '[MTD Explorer] Using optimized vectorized Spearman similarity.'\n            )\n            self.similarity_table = fast_spearman_similarity(X, Y)\n        elif fast_pearson:\n            self.logger.log_message(\n                '[MTD Explorer] Using optimized vectorized Pearson similarity.'\n            )\n            self.similarity_table = fast_pearson_similarity(X, Y)\n        else:\n            self.similarity_table = spd.cdist(\n                X, Y, metric=get_similarity_function(dist_metric)\n            )\n\n        # obtain p-values\n        self.logger.log_message('Generating the p-value table...')\n        confp = config.permute\n        if fast_spearman:\n            self.logger.log_message(\n                '[MTD Explorer] Reusing Spearman coefficients for vectorized p-values.'\n            )\n            self.pvalue_table = fast_spearman_pvalues(\n                self.similarity_table, X.shape[1]\n            )\n        elif fast_pearson:\n            self.logger.log_message(\n                '[MTD Explorer] Reusing Pearson coefficients for vectorized p-values.'\n            )\n            self.pvalue_table = fast_pearson_pvalues(\n                self.similarity_table, X.shape[1]\n            )\n        else:\n            extrapolated_time, timing_message = test_pvalue_run_time(\n                X, Y, pdist_metric=dist_metric,\n                permute_func=confp['func'], permute_iters=confp['iters'],\n                permute_speedup=confp['speedup'],\n                alpha=config.stats['fdr_alpha'],\n                force_perms=self.force_permutations,\n                num_threads=self.num_threads,\n                seed=self.seed\n            )\n            if extrapolated_time > 10 and self.verbose:\n                self.logger.log_message(timing_message)\n\n            self.pvalue_table = get_pvalue_table(\n                X, Y, pdist_metric=dist_metric,\n                permute_func=confp['func'], permute_iters=confp['iters'],\n                permute_speedup=confp['speedup'],\n                alpha=config.stats['fdr_alpha'],\n                no_progress=self.no_progress,\n                force_permutations=self.force_permutations,\n                num_threads=self.num_threads,\n                seed=self.seed\n            )\n\n'''

    text = text[:start] + new_block + text[end:]

    backup_once(main_file)
    main_file.write_text(text, encoding="utf-8")
    print(f"[INFO] Patched HAllA main module with Pearson/Spearman fast paths: {main_file}")


def patch_hierarchy(hierarchy_file: Path) -> None:
    text = hierarchy_file.read_text(encoding="utf-8")

    if (
        HIERARCHY_IMPORT_MARKER in text
        and HIERARCHY_FAST_MARKER in text
        and HIERARCHY_MEMORY_MARKER in text
    ):
        print("[OK] HAllA hierarchy.py performance patch v2 is already installed.")
        return

    if any(
        marker in text
        for marker in (
            HIERARCHY_IMPORT_MARKER,
            HIERARCHY_FAST_MARKER,
            HIERARCHY_MEMORY_MARKER,
        )
    ):
        fail("A partial HAllA hierarchy.py performance patch v2 was detected.")

    old_import = "from .utils.similarity import get_similarity_function, similarity2distance\n"
    new_import = old_import + (
        f"{HIERARCHY_IMPORT_MARKER}\n"
        "from .utils.mtd_fast_correlation import (\n"
        "    can_use_fast_pearson,\n"
        "    can_use_fast_spearman,\n"
        "    fast_pearson_condensed_distance,\n"
        "    fast_spearman_condensed_distance,\n"
        ")\n"
    )
    text = replace_once(text, old_import, new_import, "hierarchy.py import")

    old_block = '''        self.distance_matrix = similarity2distance(spd.pdist(matrix, metric=treemetric),\n                                                   sim2dist_set_abs,\n                                                   sim2dist_func)\n        self.distance_matrix = np.clip(self.distance_matrix, a_min=0, a_max=None)\n        self.distance_matrix_sqr = spd.squareform(self.distance_matrix)\n        self._generate_hierarchical_clusters(linkage_method)\n'''

    new_block = f'''        {HIERARCHY_FAST_MARKER}\n        if pdist_metric == 'spearman' and can_use_fast_spearman(matrix):\n            self.distance_matrix = fast_spearman_condensed_distance(\n                matrix,\n                sim2dist_set_abs=sim2dist_set_abs,\n                sim2dist_func=sim2dist_func,\n            )\n        elif pdist_metric == 'pearson' and can_use_fast_pearson(matrix):\n            self.distance_matrix = fast_pearson_condensed_distance(\n                matrix,\n                sim2dist_set_abs=sim2dist_set_abs,\n                sim2dist_func=sim2dist_func,\n            )\n        else:\n            self.distance_matrix = similarity2distance(\n                spd.pdist(matrix, metric=treemetric),\n                sim2dist_set_abs,\n                sim2dist_func\n            )\n\n        # Match the upstream non-negative clipping without allocating another\n        # O(n^2) condensed array.\n        np.maximum(self.distance_matrix, 0, out=self.distance_matrix)\n\n        {HIERARCHY_MEMORY_MARKER}\n        # distance_matrix_sqr is never consumed by HAllA 0.8.20.  Preserve the\n        # attribute for compatibility without allocating the square O(n^2) copy.\n        self.distance_matrix_sqr = None\n        self._generate_hierarchical_clusters(linkage_method)\n'''

    text = replace_once(text, old_block, new_block, "hierarchical distance")

    backup_once(hierarchy_file)
    hierarchy_file.write_text(text, encoding="utf-8")
    print(
        f"[INFO] Patched HAllA hierarchy module with Pearson/Spearman fast paths: "
        f"{hierarchy_file}"
    )


def write_helper(helper_file: Path) -> None:
    if helper_file.exists():
        current = helper_file.read_text(encoding="utf-8")
        if current == HELPER_SOURCE:
            print("[OK] HAllA fast correlation helper v2 is already installed.")
            return
        fail(
            f"An unexpected existing {HELPER_FILENAME} was found. "
            "Refusing to overwrite it."
        )

    helper_file.write_text(HELPER_SOURCE, encoding="utf-8")
    print(f"[INFO] Installed HAllA fast Pearson/Spearman helper: {helper_file}")


def load_helper_for_test(helper_file: Path):
    spec = importlib.util.spec_from_file_location(
        "mtd_halla_fast_correlation_validation",
        str(helper_file),
    )
    if spec is None or spec.loader is None:
        fail(f"Could not load helper module for validation: {helper_file}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _max_abs_difference(left, right) -> float:
    import numpy as np

    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def run_numerical_microtests(helper_file: Path) -> None:
    try:
        import numpy as np
        import scipy.spatial.distance as spd
        from scipy.stats import pearsonr, spearmanr
    except Exception as error:
        fail(f"Could not import NumPy/SciPy for patch validation: {error}")

    helper = load_helper_for_test(helper_file)

    # Deterministic data contain ties to exercise Spearman average-rank
    # behavior and nontrivial Pearson values.
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
    if not helper.can_use_fast_pearson(X, Y):
        fail("Fast Pearson helper unexpectedly rejected valid microtest data.")

    # Spearman X-by-Y similarity and p-values.
    fast_s_rho = helper.fast_spearman_similarity(X, Y)
    fast_s_p = helper.fast_spearman_pvalues(fast_s_rho, X.shape[1])
    ref_s_rho = np.empty(fast_s_rho.shape, dtype=float)
    ref_s_p = np.empty(fast_s_p.shape, dtype=float)

    # Pearson X-by-Y similarity and p-values.
    fast_p_rho = helper.fast_pearson_similarity(X, Y)
    fast_p_p = helper.fast_pearson_pvalues(fast_p_rho, X.shape[1])
    ref_p_rho = np.empty(fast_p_rho.shape, dtype=float)
    ref_p_p = np.empty(fast_p_p.shape, dtype=float)

    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            ref_s_rho[i, j], ref_s_p[i, j] = spearmanr(X[i, :], Y[j, :])
            ref_p_rho[i, j], ref_p_p[i, j] = pearsonr(X[i, :], Y[j, :])

    spearman_rho_diff = _max_abs_difference(fast_s_rho, ref_s_rho)
    spearman_p_diff = _max_abs_difference(fast_s_p, ref_s_p)
    pearson_rho_diff = _max_abs_difference(fast_p_rho, ref_p_rho)
    pearson_p_diff = _max_abs_difference(fast_p_p, ref_p_p)

    reference_s_condensed = 1.0 - np.abs(
        spd.pdist(X, metric=lambda a, b: spearmanr(a, b)[0])
    )
    fast_s_condensed = helper.fast_spearman_condensed_distance(X)
    spearman_distance_diff = _max_abs_difference(
        fast_s_condensed, reference_s_condensed
    )

    reference_p_condensed = 1.0 - np.abs(
        spd.pdist(X, metric=lambda a, b: pearsonr(a, b)[0])
    )
    fast_p_condensed = helper.fast_pearson_condensed_distance(X)
    pearson_distance_diff = _max_abs_difference(
        fast_p_condensed, reference_p_condensed
    )

    tolerance = 1e-12
    checks = (
        ("Spearman similarity", spearman_rho_diff),
        ("Spearman p-value", spearman_p_diff),
        ("Spearman hierarchy distance", spearman_distance_diff),
        ("Pearson similarity", pearson_rho_diff),
        ("Pearson p-value", pearson_p_diff),
        ("Pearson hierarchy distance", pearson_distance_diff),
    )
    for label, difference in checks:
        if difference > tolerance:
            fail(
                f"Fast {label} failed equivalence microtest: "
                f"max abs diff={difference:.3e}"
            )

    invalid_nan = X.copy()
    invalid_nan[0, 0] = np.nan
    invalid_constant = X.copy()
    invalid_constant[0, :] = 1.0

    for eligibility_name, eligibility in (
        ("Spearman", helper.can_use_fast_spearman),
        ("Pearson", helper.can_use_fast_pearson),
    ):
        if eligibility(invalid_nan):
            fail(f"Fast {eligibility_name} helper did not reject missing values.")
        if eligibility(invalid_constant):
            fail(f"Fast {eligibility_name} helper did not reject a constant feature.")

    print("[OK] Fast Pearson/Spearman numerical equivalence microtests passed.")
    print(f"[INFO] Spearman max similarity difference: {spearman_rho_diff:.3e}")
    print(f"[INFO] Spearman max p-value difference:    {spearman_p_diff:.3e}")
    print(f"[INFO] Spearman max distance difference:   {spearman_distance_diff:.3e}")
    print(f"[INFO] Pearson max similarity difference:  {pearson_rho_diff:.3e}")
    print(f"[INFO] Pearson max p-value difference:     {pearson_p_diff:.3e}")
    print(f"[INFO] Pearson max distance difference:    {pearson_distance_diff:.3e}")


def validate_patch(package_dir: Path) -> None:
    main_file = package_dir / "main.py"
    hierarchy_file = package_dir / "hierarchy.py"
    helper_file = package_dir / "utils" / HELPER_FILENAME

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
        fail("HAllA fast correlation helper marker is missing.")

    expected_sha = sha256_text(HELPER_SOURCE)
    observed_sha = sha256_text(helper_text)
    if observed_sha != expected_sha:
        fail("HAllA fast correlation helper contents do not match patch v2.")

    # A completed v2 upgrade should no longer leave the v1 helper active.
    old_helper = package_dir / "utils" / V1_HELPER_FILENAME
    if old_helper.exists():
        fail(f"Obsolete v1 HAllA performance helper is still present: {old_helper}")

    for path in (main_file, hierarchy_file, helper_file):
        compile_module(path)

    run_numerical_microtests(helper_file)

    print("[OK] HAllA 0.8.20 Pearson/Spearman performance patch v2 validation passed.")
    print(f"[INFO] HAllA package: {package_dir}")


def apply_patch(package_dir: Path) -> None:
    main_file = package_dir / "main.py"
    hierarchy_file = package_dir / "hierarchy.py"
    helper_file = package_dir / "utils" / HELPER_FILENAME

    for path in (main_file, hierarchy_file):
        if not path.is_file():
            fail(f"Expected HAllA 0.8.20 module was not found: {path}")

    upgrade_v1_if_needed(package_dir)

    # Write the helper first; main.py/hierarchy.py will not import it until a
    # future HAllA process starts. Validation occurs after both callers are
    # patched.
    write_helper(helper_file)
    patch_main(main_file)
    patch_hierarchy(hierarchy_file)
    validate_patch(package_dir)

    print("[OK] HAllA 0.8.20 Pearson/Spearman performance patch v2 applied.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Apply or validate the MTD Explorer HAllA 0.8.20 "
            "Pearson/Spearman performance patch v2."
        )
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate patch v2 without modifying HAllA.",
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
