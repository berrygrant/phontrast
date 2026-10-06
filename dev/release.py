#!/usr/bin/env python3
"""Release helper for phontrast.

Bumps the version across the metadata files, builds the tarball, runs the
CRAN checks, triggers R-hub, and tags the release. Standard library only;
needs R on the PATH (and the GitHub CLI `gh` for the `rhub` command).

    python3 dev/release.py bump 2.5.0            # DESCRIPTION, CITATION.cff (+ today), README pin; NEWS heading check
    python3 dev/release.py check                 # R CMD build + R CMD check --as-cran --timings, then a summary
    python3 dev/release.py bump 2.5.0 --check    # both
    python3 dev/release.py report <dir.Rcheck>   # summarize an existing check directory
    python3 dev/release.py rhub [--platforms linux,windows,macos,mkl]
    python3 dev/release.py tag 2.5.0 [--yes]     # verify the metadata agree, then create and push vX.Y.Z
    python3 dev/release.py doi 2.5.0 10.5281/zenodo.NNNNNNN
                                                 # record the version DOI Zenodo minted: inst/CITATION map,
                                                 # test-citation.R, CITATION.cff (doi + url; concept DOI kept
                                                 # under identifiers)

Run it from anywhere inside the repository. `bump` and `doi` are idempotent.

Release order that gets the version-specific DOI into the CRAN tarball:
bump -> check -> rhub -> tag -> GitHub release (Zenodo mints the DOI) -> doi
-> commit -> check -> submit to CRAN.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESCRIPTION = ROOT / "DESCRIPTION"
CFF = ROOT / "CITATION.cff"
README = ROOT / "README.md"
NEWS = ROOT / "NEWS.md"
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+(\.\d+)?$")


def die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    print("$ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=kw.pop("cwd", ROOT), **kw)


def description_version() -> str:
    m = re.search(r"^Version:\s*(\S+)\s*$", DESCRIPTION.read_text(), re.M)
    if not m:
        die("no Version: field in DESCRIPTION")
    return m.group(1)


def cff_version() -> str:
    m = re.search(r"^version:\s*\"?([^\"\n]+)\"?\s*$", CFF.read_text(), re.M)
    return m.group(1).strip() if m else ""


def readme_pin() -> str:
    m = re.search(r"phontrast@v([0-9][0-9.]*)", README.read_text())
    return m.group(1) if m else ""


def news_has_heading(version: str) -> bool:
    return re.search(rf"^# phontrast {re.escape(version)}\s*$", NEWS.read_text(), re.M) is not None


# ---- bump --------------------------------------------------------------------

def substitute(path: Path, pattern: str, replacement: str, dry_run: bool) -> bool:
    text = path.read_text()
    new, n = re.subn(pattern, replacement, text, flags=re.M)
    if n == 0:
        die(f"{path.relative_to(ROOT)}: pattern not found: {pattern}")
    changed = new != text
    label = "unchanged" if not changed else ("would update" if dry_run else "updated")
    print(f"  {label:12s} {path.relative_to(ROOT)}")
    if changed and not dry_run:
        path.write_text(new)
    return changed


def cmd_bump(args: argparse.Namespace) -> None:
    version = args.version
    if not VERSION_RE.match(version):
        die(f"'{version}' is not an X.Y.Z version")
    today = dt.date.today().isoformat()
    print(f"Setting version {version} (release date {today}):")
    substitute(DESCRIPTION, r"^Version:\s*.*$", f"Version: {version}", args.dry_run)
    substitute(CFF, r"^version:\s*.*$", f"version: {version}", args.dry_run)
    substitute(CFF, r"^date-released:\s*.*$", f"date-released: {today}", args.dry_run)
    substitute(README, r"phontrast@v[0-9][0-9.]*", f"phontrast@v{version}", args.dry_run)
    if news_has_heading(version):
        print(f"  ok           NEWS.md has a '# phontrast {version}' heading")
    else:
        top = next((l for l in NEWS.read_text().splitlines() if l.startswith("# ")), "(none)")
        die(f"NEWS.md has no '# phontrast {version}' heading (top heading: {top})")
    print(f"inst/CITATION is version-aware; once Zenodo mints this version's DOI run "
          f"`python3 dev/release.py doi {version} 10.5281/zenodo.NNNNNNN`.")
    if args.check:
        cmd_check(args)


# ---- doi ---------------------------------------------------------------------

INST_CITATION = ROOT / "inst" / "CITATION"
TEST_CITATION = ROOT / "tests" / "testthat" / "test-citation.R"
DOI_RE = re.compile(r"^10\.5281/zenodo\.\d+$")


def concept_doi() -> str:
    m = re.search(r'if \(is\.na\(doi\)\) \{\s*doi <- "([^"]+)"', INST_CITATION.read_text())
    if not m:
        die("could not find the concept-DOI fallback in inst/CITATION")
    return m.group(1)


def update_doi_map(path: Path, var_name: str, version: str, doi: str, dry_run: bool) -> None:
    text = path.read_text()
    # The block may sit inside a test_that() call, so keep its indentation.
    block = re.search(
        rf"^(?P<indent>[ \t]*){re.escape(var_name)} <- c\((?P<body>.*?)^[ \t]*\)",
        text, re.M | re.S,
    )
    if not block:
        die(f"{path.relative_to(ROOT)}: no `{var_name} <- c(...)` block")
    indent = block.group("indent")
    entries = dict(re.findall(r'"(\d+\.\d+\.\d+(?:\.\d+)?)"\s*=\s*"([^"]+)"', block.group("body")))
    entries[version] = doi
    ordered = sorted(entries.items(), key=lambda kv: tuple(int(p) for p in kv[0].split(".")))
    body = ",\n".join(f'{indent}  "{v}" = "{d}"' for v, d in ordered)
    new = text[:block.start()] + f"{indent}{var_name} <- c(\n{body}\n{indent})" + text[block.end():]
    changed = new != text
    label = "unchanged" if not changed else ("would update" if dry_run else "updated")
    print(f"  {label:12s} {path.relative_to(ROOT)}")
    if changed and not dry_run:
        path.write_text(new)


def cmd_doi(args: argparse.Namespace) -> None:
    version, doi = args.version, args.doi.strip()
    if not VERSION_RE.match(version):
        die(f"'{version}' is not an X.Y.Z version")
    if doi.lower().startswith("https://doi.org/"):
        doi = doi[len("https://doi.org/"):]
    if not DOI_RE.match(doi):
        die(f"'{doi}' does not look like a Zenodo DOI (10.5281/zenodo.NNNNNNN)")
    concept = concept_doi()
    if doi == concept:
        die(f"{doi} is the concept (all-versions) DOI; pass the DOI Zenodo minted for {version}")
    print(f"Recording {doi} as the DOI of phontrast {version}:")
    update_doi_map(INST_CITATION, "version_dois", version, doi, args.dry_run)
    update_doi_map(TEST_CITATION, "expected_dois", version, doi, args.dry_run)
    if not args.no_cff:
        if cff_version() != version:
            die(f"CITATION.cff is at version {cff_version() or '(none)'}, not {version}; run `bump {version}` first")
        substitute(CFF, r"^doi:\s*.*$", f"doi: {doi}", args.dry_run)
        substitute(CFF, r"^url:\s*.*$", f'url: "https://doi.org/{doi}"', args.dry_run)
        if "identifiers:" not in CFF.read_text():
            block = ("identifiers:\n"
                     "  - type: doi\n"
                     f"    value: {concept}\n"
                     "    description: Concept DOI for all versions\n")
            label = "would update" if args.dry_run else "updated"
            print(f"  {label:12s} CITATION.cff (identifiers block with the concept DOI)")
            if not args.dry_run:
                CFF.write_text(CFF.read_text().rstrip("\n") + "\n" + block)
    print("Commit the result; if the CRAN tarball has not been built yet, build it after this so "
          "citation(\"phontrast\") in the released package carries the version DOI.")


# ---- check -------------------------------------------------------------------

def build_tarball() -> Path:
    version = description_version()
    tarball = ROOT / f"phontrast_{version}.tar.gz"
    if tarball.exists():
        tarball.unlink()
    res = run(["R", "CMD", "build", ".", "--no-resave-data"])
    if res.returncode != 0 or not tarball.exists():
        die("R CMD build failed")
    return tarball


def cmd_check(args: argparse.Namespace) -> None:
    tarball = build_tarball()
    check_dir = ROOT / ".release-check"
    if check_dir.exists():
        shutil.rmtree(check_dir)
    check_dir.mkdir()
    moved = check_dir / tarball.name
    shutil.move(str(tarball), moved)

    cmd = ["R", "CMD", "check", "--timings"]
    if not getattr(args, "no_as_cran", False):
        cmd.append("--as-cran")
    if getattr(args, "no_manual", False):
        cmd.append("--no-manual")
    cmd.append(moved.name)
    env = dict(os.environ)
    if getattr(args, "no_force_suggests", False):
        env["_R_CHECK_FORCE_SUGGESTS_"] = "false"
    res = run(cmd, cwd=check_dir, env=env)
    rcheck = check_dir / "phontrast.Rcheck"
    status = summarize_check(rcheck)
    print(f"\nTarball: {moved}\nCheck directory: {rcheck}")
    if res.returncode != 0 or status["errors"] or status["warnings"]:
        sys.exit(1)


def summarize_check(rcheck: Path) -> dict:
    log = rcheck / "00check.log"
    if not log.exists():
        die(f"no 00check.log in {rcheck}")
    lines = log.read_text(errors="replace").splitlines()
    flagged: list[tuple[str, list[str]]] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.search(r"\.\.\.\s*(NOTE|WARNING|ERROR)\s*$", line)
        if m:
            detail = []
            j = i + 1
            while j < len(lines) and not lines[j].startswith("* ") and not lines[j].startswith("Status:"):
                detail.append(lines[j])
                j += 1
            flagged.append((line.strip(), detail))
            i = j
            continue
        i += 1
    status_line = next((l for l in lines if l.startswith("Status:")), "Status: (not found)")

    print("\n=== R CMD check summary ===")
    print(status_line)
    for head, detail in flagged:
        print(f"\n{head}")
        for d in detail:
            if d.strip():
                print(f"    {d}")

    timings = rcheck / "phontrast-Ex.timings"
    if timings.exists():
        slow = []
        for row in timings.read_text().splitlines()[1:]:
            parts = row.split("\t")
            if len(parts) >= 4:
                try:
                    user, system, elapsed = (float(parts[1]), float(parts[2]), float(parts[3]))
                except ValueError:
                    continue
                if max(user + system, elapsed) > 5.0:
                    slow.append((parts[0], user + system, elapsed))
        print("\nExamples over 5 s (CRAN flags these):" if slow else "\nNo example over 5 s.")
        for name, cpu, elapsed in slow:
            print(f"    {name:28s} cpu {cpu:6.2f}s  elapsed {elapsed:6.2f}s")

    tests_out = sorted((rcheck / "tests").glob("testthat.Rout*")) if (rcheck / "tests").exists() else []
    for out in tests_out:
        tail = [l for l in out.read_text(errors="replace").splitlines() if l.startswith("[ FAIL")]
        if tail:
            print(f"\nTests ({out.name}): {tail[-1]}")

    kinds = [h.split()[-1] for h, _ in flagged]
    return {"errors": kinds.count("ERROR"), "warnings": kinds.count("WARNING"), "notes": kinds.count("NOTE")}


def cmd_report(args: argparse.Namespace) -> None:
    summarize_check(Path(args.rcheck_dir).resolve())


# ---- rhub --------------------------------------------------------------------

def cmd_rhub(args: argparse.Namespace) -> None:
    if shutil.which("gh") is None:
        die("the GitHub CLI `gh` is not on the PATH. Alternative from R: "
            "rhub::rhub_check(platforms = c(\"linux\", \"windows\", \"macos\", \"mkl\"))")
    name = args.name or f"phontrast {description_version()} release check"
    res = run(["gh", "workflow", "run", "rhub.yaml", "-f", f"config={args.platforms}", "-f", f"name={name}"])
    if res.returncode != 0:
        die("gh workflow run failed (is the rhub.yaml workflow on the default branch, and is RHUB_TOKEN set?)")
    print("Triggered. Follow it with: gh run list --workflow rhub.yaml -L 3   (or on the Actions tab)")


# ---- tag ---------------------------------------------------------------------

def cmd_tag(args: argparse.Namespace) -> None:
    version = args.version
    if not VERSION_RE.match(version):
        die(f"'{version}' is not an X.Y.Z version")
    problems = []
    if description_version() != version:
        problems.append(f"DESCRIPTION has {description_version()}")
    if cff_version() != version:
        problems.append(f"CITATION.cff has {cff_version() or '(none)'}")
    if readme_pin() != version:
        problems.append(f"README install pin is @v{readme_pin() or '(none)'}")
    if not news_has_heading(version):
        problems.append("NEWS.md has no heading for this version")
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    if dirty:
        problems.append("working tree has uncommitted changes")
    if problems:
        die("not tagging:\n  - " + "\n  - ".join(problems) + f"\nRun `python3 dev/release.py bump {version}` and commit first.")
    tag = f"v{version}"
    if subprocess.run(["git", "rev-parse", "-q", "--verify", f"refs/tags/{tag}"], cwd=ROOT, capture_output=True).returncode == 0:
        die(f"tag {tag} already exists locally")
    if args.dry_run:
        print(f"would run: git tag -a {tag} -m 'phontrast {version}' && git push origin {tag}")
        return
    if not args.yes:
        answer = input(f"Create and push tag {tag} at HEAD? [y/N] ").strip().lower()
        if answer not in ("y", "yes"):
            print("aborted")
            return
    if run(["git", "tag", "-a", tag, "-m", f"phontrast {version}"]).returncode != 0:
        die("git tag failed")
    if run(["git", "push", "origin", tag]).returncode != 0:
        die("git push failed")
    print(f"Pushed {tag}. The R-CMD-check workflow verifies the tag against DESCRIPTION and CITATION.cff.")


# ---- main --------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("bump", help="set the version in DESCRIPTION, CITATION.cff, README")
    p.add_argument("version")
    p.add_argument("--check", action="store_true", help="then build and check")
    p.add_argument("--dry-run", action="store_true")
    for flag, help_text in [("--no-as-cran", "omit --as-cran"), ("--no-manual", "skip the PDF manual (no TeX)"),
                            ("--no-force-suggests", "do not require every Suggests package")]:
        p.add_argument(flag, action="store_true", help=help_text)
    p.set_defaults(func=cmd_bump)

    p = sub.add_parser("check", help="R CMD build + R CMD check --as-cran --timings")
    for flag, help_text in [("--no-as-cran", "omit --as-cran"), ("--no-manual", "skip the PDF manual (no TeX)"),
                            ("--no-force-suggests", "do not require every Suggests package")]:
        p.add_argument(flag, action="store_true", help=help_text)
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("report", help="summarize an existing .Rcheck directory")
    p.add_argument("rcheck_dir")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("rhub", help="trigger the R-hub GitHub Actions workflow")
    p.add_argument("--platforms", default="linux,windows,macos", help="comma-separated R-hub platforms (add mkl for CRAN's MKL flavor)")
    p.add_argument("--name", default=None)
    p.set_defaults(func=cmd_rhub)

    p = sub.add_parser("tag", help="verify the metadata and create/push the vX.Y.Z tag")
    p.add_argument("version")
    p.add_argument("--yes", action="store_true", help="do not ask for confirmation")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=cmd_tag)

    p = sub.add_parser("doi", help="record the version DOI Zenodo minted for a release")
    p.add_argument("version")
    p.add_argument("doi", help="10.5281/zenodo.NNNNNNN (a https://doi.org/ prefix is accepted)")
    p.add_argument("--no-cff", action="store_true", help="leave CITATION.cff alone")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=cmd_doi)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
