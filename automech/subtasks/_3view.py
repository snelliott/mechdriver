"""Visualize results."""

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
    disp_ = {"find_ts": _display_find_ts}.get(task_key)

    for task_path in fs.task_paths(task_key, subtask_key, path=path, runlvl=runlvl):
        disp_(task_path)


def _display_find_ts(task_path: str) -> None:
    """Display results for a "find_ts" task.

    :param task_key: Task key
    :param subtask_key: Subtask key
    :param path: Path to AutoMech run directory, containing `inp/` folder
    :param runlvl: Specify the run level (needed if there are duplicate tasks)
    """
    cnf_fs = autofile.fs.conformer(task_path)
    for cnf_loc in cnf_fs[-1].existing():
        geo = cnf_fs[-1].file.geometry.read(cnf_loc)
        hess = cnf_fs[-1].file.hessian.read(cnf_loc)
        freqs, norm_coos = automol.geom.vibrational_analysis(geo, hess)
        print("Lowest frequency mode:")
        print(f"  frequency: {freqs[0]}")
        automol.geom.display(geo, mode=norm_coos[:, 0])
