""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

import itertools
import subprocess
import tarfile
from collections.abc import Sequence
from pathlib import Path

import pandas
import yaml

from ..base import Status
from ._0setup import INFO_FILE, SUBTASK_DIR, SubtasksInfo, Task
from ._1status import log_paths_with_check_results, parse_subtask_status

SCRIPT_DIR = Path(__file__).parent / "scripts"
RUN_SCRIPT = str(SCRIPT_DIR / "run_adhoc.sh")


def run(
    paths: Sequence[str | Path] = (".",),
    nodes: Sequence[str] | None = None,
    dir_name: str = SUBTASK_DIR,
    activation_hook: str | None = None,
    statuses: Sequence[Status] = (Status.TBD,),
    tar: bool = False,
) -> None:
    """Runs subtasks in parallel on Ad Hoc cluster

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param nodes: A comma-separated list of nodes to run on
    :param activation_hook: Shell commands for activating the AutoMech environment on the remote
    :param statuses: A comma-separated list of status to run or re-run
    :param tar: Tar the subtask data and save filesystem after running?
    """
    # Determine paths
    paths = [Path(p) / dir_name for p in paths]
    for path in paths:
        assert (
            path.exists()
        ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    # Read in subtask information
    infos = [SubtasksInfo(**yaml.safe_load((p / INFO_FILE).read_text())) for p in paths]

    # Make sure the run and save directories exist
    for info in infos:
        (info.root_path / info.run_path).mkdir(exist_ok=True)
        (info.root_path / info.save_path).mkdir(exist_ok=True)

    # Resolve paths relative to the current working directory
    infos: list[SubtasksInfo] = [p / i for p, i in zip(paths, infos, strict=True)]

    # Zip tasks together
    task_lsts = list(
        itertools.zip_longest(*(itertools.chain(*i.task_groups) for i in infos))
    )
    print(len(task_lsts))

    import sys

    sys.exit()

    # Make sure the paths and their run and save directories exist
    group_ids = set()
    for path in paths:
        assert (
            path.exists()
        ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

        info_path = path / INFO_FILE
        info = SubtasksInfo(**yaml.safe_load(info_path.read_text()))

        group_ids |= set(info.group_ids)
        run_path = Path(info.root_path) / info.run_path
        save_path = Path(info.root_path) / info.save_path

        run_path.mkdir(exist_ok=True)
        save_path.mkdir(exist_ok=True)

    # for group_id in group_ids:
    #     tasks_lst = list(
    #         itertools.zip_longest(
    #             *(read_task_list(p / f"{group_id}.yaml") for p in paths)
    #         )
    #     )
    #     for tasks in tasks_lst:

    #     print(len(paths))
    #     print(len(tasks_lst))

    import sys

    sys.exit()

    for group_id in group_ids:
        tasks = read_task_list(paths / f"{group_id}.yaml")
        if not tasks:
            continue

        df = pandas.read_csv(paths / f"{group_id}.csv")
        for task_key, row in df.iterrows():
            task: Task = tasks[task_key]
            assert row[TableKey.task] == task.name, f"{row} does not match {task.name}"

            subtask_paths = []
            subtask_logs = []

            for key, nworkers in zip(
                task.subtask_keys, task.subtask_nworkers, strict=True
            ):
                assert key in row, f"Key {key} not present in row:\n{row}"
                subtask_path = row.get(key)
                status = parse_subtask_status(
                    log_paths_with_check_results(subtask_path)
                )
                if status in statuses:
                    subtask_paths.extend([subtask_path] * nworkers)
                    subtask_logs.extend([f"out{i}.log" for i in range(nworkers)])

            if subtask_paths:
                run_args = [
                    RUN_SCRIPT,
                    work_path,
                    f"{task.mem}",
                    f"{task.nprocs}",
                    ",".join(subtask_paths),
                    ",".join(subtask_logs),
                    ",".join(nodes),
                    "" if activation_hook is None else activation_hook,
                ]
                subprocess.run(run_args)

    if tar:
        tar_subtask_data(paths)


def tar_subtask_data(path: str | Path = SUBTASK_DIR) -> None:
    """Tar the save directory for a subtask run

    :param path: The path where the AutoMech subtasks were set up
    """
    path = Path(path)
    assert (
        path.exists()
    ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    info_path = path / INFO_FILE
    info_dct = yaml.safe_load(info_path.read_text())
    save_path = Path(info_dct[InfoKey.save_path])

    tar_directory(path)
    tar_directory(save_path)


def untar_subtask_data(path: str | Path = SUBTASK_DIR) -> None:
    """Un-tar the save directory for a subtask run, if it exists

    :param path: The path where the AutoMech subtasks were set up
    """
    untar_directory(path)

    path = Path(path)
    assert (
        path.exists()
    ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    info_path = path / INFO_FILE
    info_dct = yaml.safe_load(info_path.read_text())
    save_path = Path(info_dct[InfoKey.save_path])

    untar_directory(save_path)


def tar_directory(dir_path: str | Path) -> None:
    """Tar a directory in its current location

    :param path: The path to a directory
    """
    dir_path = Path(dir_path)
    tar_path = dir_path.with_suffix(".tgz")
    print(f"Tarring {dir_path} into {tar_path}...")
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(dir_path, arcname=dir_path.name)


def untar_directory(dir_path: str | Path) -> None:
    """Un-tar a directory in its current location

    :param path: The path where the AutoMech subtasks were set up
    """
    dir_path = Path(dir_path)
    tar_path = dir_path.with_suffix(".tgz")
    if tar_path.exists():
        print(f"Un-tarring {tar_path} into {dir_path}...")
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(dir_path.parent)
