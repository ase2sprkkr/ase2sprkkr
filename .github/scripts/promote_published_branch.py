"""Advance distribution branches only for the latest successful PyPI publication.

The publishing workflow runs this script *after* the PyPI upload job succeeds.
The exact tag commit, rather than the head of develop, becomes the branch tip.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time
import tomllib
from urllib.error import URLError
from urllib.request import Request, urlopen

from packaging.version import InvalidVersion, Version


PYPI_JSON = "https://pypi.org/pypi/ase2sprkkr/json"


def target_branch(version: Version) -> str:
    """Prereleases live on master; final releases live on release."""
    return "master" if version.is_prerelease else "release"


def published_versions(releases: dict) -> set[Version]:
    """Ignore yanked releases and versions without downloadable files."""
    result = set()
    for name, files in releases.items():
        try:
            version = Version(name)
        except InvalidVersion:
            continue
        if files and any(not file.get("yanked", False) for file in files):
            result.add(version)
    return result


def latest_in_channel(versions: set[Version], version: Version) -> Version | None:
    same_channel = (v for v in versions if v.is_prerelease == version.is_prerelease)
    return max(same_channel, default=None)


def can_advance(current: Version, candidate: Version) -> bool:
    """Never replace a target branch already on the same or newer version."""
    return current < candidate


def source_version_and_branch(tag: str, pyproject: str) -> tuple[Version, str]:
    project = tomllib.loads(pyproject)["project"]
    if project["name"] != "ase2sprkkr":
        raise ValueError("Unexpected project name; refusing to move a branch")
    raw_version = project["version"]
    if tag != f"v{raw_version}":
        raise ValueError(f"Tag {tag!r} does not match pyproject.toml version {raw_version!r}")
    version = Version(raw_version)
    return version, target_branch(version)


def index_versions() -> set[Version]:
    request = Request(
        PYPI_JSON,
        headers={"Accept": "application/json", "Cache-Control": "no-cache",
                 "User-Agent": "ase2sprkkr-release-promotion"},
    )
    with urlopen(request, timeout=20) as response:
        return published_versions(json.load(response)["releases"])


def candidate_is_latest(version: Version, *, wait_for_index: bool) -> bool:
    """Do not mistake PyPI index propagation delay for a missing publication."""
    attempts = 12 if wait_for_index else 1
    for attempt in range(attempts):
        try:
            versions = index_versions()
        except (URLError, TimeoutError, OSError) as exc:
            if attempt + 1 == attempts:
                raise RuntimeError("Cannot verify PyPI release index") from exc
            print(f"PyPI index unavailable ({exc}); retrying", flush=True)
        else:
            latest = latest_in_channel(versions, version)
            if latest is not None and latest > version:
                print(f"Skipping {version}: newer published version {latest} exists", flush=True)
                return False
            if version in versions:
                return True
            if attempt + 1 == attempts:
                raise RuntimeError(f"Published version {version} was not found on PyPI")
            print(f"Waiting for PyPI to index {version}", flush=True)
        time.sleep(10)
    raise AssertionError("Unreachable")


def git(*arguments: str) -> str:
    return subprocess.check_output(["git", *arguments], text=True).strip()


def promote() -> None:
    version, branch = source_version_and_branch(
        os.environ["RELEASE_TAG"], Path("pyproject.toml").read_text(encoding="utf-8")
    )
    tagged_commit = git("rev-parse", "HEAD")
    print(f"Published {version}: target {branch} -> {tagged_commit}", flush=True)

    if not candidate_is_latest(version, wait_for_index=True):
        return

    for attempt in range(3):
        # Refresh every attempt. A concurrent publication can move the ref.
        git("fetch", "--no-tags", "origin",
            f"+refs/heads/{branch}:refs/remotes/origin/{branch}")
        remote_ref = f"refs/remotes/origin/{branch}"
        previous_sha = git("rev-parse", remote_ref)
        if previous_sha == tagged_commit:
            print(f"{branch} already points to {tagged_commit}", flush=True)
            return

        branch_pyproject = git("show", f"{remote_ref}:pyproject.toml")
        branch_version = Version(tomllib.loads(branch_pyproject)["project"]["version"])
        if not can_advance(branch_version, version):
            print(f"Skipping: {branch} already has version {branch_version}", flush=True)
            return
        if not candidate_is_latest(version, wait_for_index=False):
            return

        # This may be non-fast-forward (master currently diverges from develop).
        # Lease prevents an older workflow run from overwriting a newly promoted ref.
        command = [
            "git", "push", f"--force-with-lease=refs/heads/{branch}:{previous_sha}",
            "origin", f"HEAD:refs/heads/{branch}",
        ]
        result = subprocess.run(command, text=True, capture_output=True)
        if result.returncode == 0:
            print(f"Updated {branch} to published tag {os.environ['RELEASE_TAG']}", flush=True)
            return
        print(result.stderr.strip(), file=sys.stderr, flush=True)
        if attempt == 2:
            raise RuntimeError(f"Unable to promote {branch}; check branch protection/permissions")
        print("Ref changed or push was rejected; rechecking before retry", flush=True)
        time.sleep(2)


if __name__ == "__main__":
    promote()
