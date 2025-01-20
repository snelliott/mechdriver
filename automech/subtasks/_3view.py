"""Visualize results."""

from collections.abc import Callable
from pathlib import Path

import autofile
import automol

from . import fs


def display(
    task_key: str, subtask_key: str, path: str | Path = ".", runlvl: str | None = None
) -> None:
    """Display results for a task.

    :param task_key: Task key
    :param subtask_key: Subtask key
    :param path: Path to AutoMech run directory, containing `inp/` folder
    :param runlvl: Specify the run level (needed if there are duplicate tasks)
    """
    if not task_has_display_function(task_key):
        raise NotImplementedError(f"Display for {task_key} is not yet implemented...")

    disp_ = task_display_function(task_key)

    for task_path in fs.task_paths(task_key, subtask_key, path=path, runlvl=runlvl):
        disp_(task_path)


def _display_find_ts(task_path: str) -> None:
    """Display results for a "find_ts" task.

    :param task_path: Task path
    """
    cnf_fs = autofile.fs.conformer(task_path)
    for cnf_loc in cnf_fs[-1].existing():
        geo = cnf_fs[-1].file.geometry.read(cnf_loc)
        hess = cnf_fs[-1].file.hessian.read(cnf_loc)
        freqs, norm_coos = automol.geom.vibrational_analysis(geo, hess)
        print("Lowest frequency mode:")
        print(f"  frequency: {freqs[0]}")
        automol.geom.display(geo, mode=norm_coos[:, 0])


def _display_rpath_scan(task_path: str) -> None:
    """Display results for a "find_ts" task.

    :param task_path: Task path
    """
    for scan_fs in autofile.fs.iterate_managers(
        task_path, ["CONFORMER", "ZMATRIX"], "SCAN"
    ):
        irc_loc = next(loc for loc in scan_fs[1].existing() if loc[-1][0] == "IRC")
        irc_geos, *_ = zip(*scan_fs[1].file.trajectory.read(irc_loc))
        automol.geom.display_trajectory(irc_geos)


TASK_DISPLAY_FUNCTION = {"find_ts": _display_find_ts, "rpath_scan": _display_rpath_scan}


def task_has_display_function(task_key: str) -> bool:
    """Determine if task has display function.

    :param task_key: Task key
    :return: `True` if it does, otherwise `False`
    """
    return task_key in TASK_DISPLAY_FUNCTION


def task_display_function(task_key: str) -> Callable[[str], None]:
    """Get the appropriate task display function.

    :param task_key: Task key
    :return: Display function
    """
    return TASK_DISPLAY_FUNCTION.get(task_key)
