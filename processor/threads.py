"""
Thread count from THREADS (0 = the container's CPU quota, never the host's core
count), applied to the BLAS / OpenMP / numba environment before numpy is imported.
"""
from __future__ import annotations

import math
import os
import sys

_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS")


def cgroup_cpus() -> int | None:
    """CPUs from the cgroup v2 quota (cpu.max) or v1 (cfs quota/period); None without a quota."""
    try:
        quota, period = open("/sys/fs/cgroup/cpu.max").read().split()[:2]
        if quota != "max":
            return max(1, math.ceil(int(quota) / int(period)))
        return None
    except (OSError, ValueError):
        pass
    try:
        quota = int(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
        period = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
        if quota > 0:
            return max(1, math.ceil(quota / period))
    except (OSError, ValueError):
        pass
    return None


def resolve_threads(value: str | int | None) -> int:
    try:
        n = int(value) if value not in (None, "") else 0
    except ValueError:
        n = 0
    if n > 0:
        return n
    n = cgroup_cpus()
    if n is None:
        try:
            n = len(os.sched_getaffinity(0))
        except AttributeError:
            n = os.cpu_count() or 1
    return max(1, min(n, 8))


def configure(value: str | int | None) -> int:
    n = resolve_threads(value)
    numba = sys.modules.get("numba")
    for var in _VARS:
        if var == "NUMBA_NUM_THREADS" and numba is not None:
            continue    # numba already imported (tests): changing the variable now makes numba raise
        os.environ[var] = str(n)
    if numba is not None:
        numba.set_num_threads(min(n, numba.config.NUMBA_NUM_THREADS))
    return n
