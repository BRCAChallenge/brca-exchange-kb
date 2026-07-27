#!/usr/bin/env python3
"""
Generate a new BRCA Exchange data release.

This script sets up the environment and kicks off the pipeline to create
a new data release. It handles:
- Creating a working directory
- Cloning the code repository
- Generating the pipeline configuration from template
- Running the build-release Make target

Python 3.13+ compatible.
"""

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    from jinja2 import Environment, FileSystemLoader, select_autoescape
except ImportError:
    print("Error: jinja2 is required. Install with: pip install jinja2", file=sys.stderr)
    sys.exit(1)


def resolve_path(path: str) -> Path:
    """Resolve a path to its absolute form."""
    return Path(path).resolve()


def get_current_branch(script_dir: Path) -> str:
    """
    Determine the git branch of the checkout this script is running from,
    so running from a feature branch deploys that branch by default rather
    than silently falling back to master.

    Falls back to "master" if the branch can't be determined (e.g. the
    script isn't inside a git checkout, or HEAD is detached).
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=script_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        branch = result.stdout.strip()
        if branch and branch != "HEAD":
            return branch
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return "master"


def run_command(
    cmd: list[str],
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = False
) -> subprocess.CompletedProcess:
    """
    Run a shell command with proper error handling.

    Args:
        cmd: Command and arguments as a list
        cwd: Working directory for the command
        check: Whether to raise exception on non-zero exit
        capture_output: Whether to capture stdout/stderr

    Returns:
        CompletedProcess instance
    """
    print(f"Running: {' '.join(cmd)}")
    return subprocess.run(
        cmd,
        cwd=cwd,
        check=check,
        capture_output=capture_output,
        text=True
    )


def clone_or_update_repo(code_base: Path, git_commit: str) -> None:
    """
    Clone the repository if it doesn't exist, then check out the specified
    branch/commit/tag and fast-forward it to match origin.

    Args:
        code_base: Path where the code should be cloned
        git_commit: Git commit/branch/tag to checkout
    """
    repo_url = "https://github.com/BRCAChallenge/brca-exchange-kb.git"

    if not code_base.exists():
        print(f"Cloning repository to {code_base}...")
        run_command(["git", "clone", repo_url, str(code_base)])
    else:
        print(f"Repository already exists at {code_base}")
        print("Fetching latest changes from remote...")
        run_command(["git", "fetch", "origin"], cwd=code_base)

    print(f"Checking out {git_commit}...")
    # Try to checkout directly first
    try:
        run_command(["git", "checkout", git_commit], cwd=code_base)
    except subprocess.CalledProcessError:
        # If direct checkout fails, try as a remote branch
        print(f"Direct checkout failed, trying origin/{git_commit}...")
        run_command(["git", "checkout", "-b", git_commit, f"origin/{git_commit}"], cwd=code_base)
        return

    # `git checkout` above is a no-op if code_base was already on this
    # branch, so a reused work directory would otherwise keep running
    # whatever commit it happened to be on -- even after the `git fetch`
    # above -- silently ignoring newer commits on origin. Fast-forward
    # explicitly to pick those up.
    try:
        run_command(["git", "merge", "--ff-only", f"origin/{git_commit}"], cwd=code_base)
    except subprocess.CalledProcessError:
        print(f"Note: no origin/{git_commit} to fast-forward from "
              "(may be a tag or a specific commit) -- leaving as checked out.")


def generate_config(
    template_path: Path,
    output_path: Path,
    context: dict[str, str]
) -> None:
    """
    Generate configuration file from Jinja2 template.

    Args:
        template_path: Path to the Jinja2 template file
        output_path: Path where the generated config should be written
        context: Dictionary of template variables
    """
    print(f"Generating configuration file: {output_path}")

    # Set up Jinja2 environment
    template_dir = template_path.parent
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(),
        keep_trailing_newline=True
    )

    # Load and render template
    template = env.get_template(template_path.name)
    rendered = template.render(**context)

    # Write output
    output_path.write_text(rendered)
    print(f"Configuration written to {output_path}")


def spawn_pipeline(code_base: Path, work_dir: Path) -> Path:
    """
    Launch the pipeline build-release target in the background and return
    immediately -- this does not wait for the pipeline to finish.

    The gene configuration is baked into the generated brca_pipeline_cfg.mk
    (via GENE_CONFIG_FILENAME in the template context), not passed as a make
    variable here -- the Makefile no longer reads a GENE_CONFIG_FILENAME
    override directly.

    Both logs are now written by the Makefile itself: run-pipeline's Luigi
    output goes to PIPELINE_LOG, and the rest of build-release's own output
    (checkout, resource downloads, docker service startup, ...) goes to
    BUILD_RELEASE_LOG -- that's where any early, fast-failing setup problems
    would show up. Both paths are computed here rather than left to the
    Makefile's own defaults, so the caller can report them immediately
    without waiting on or parsing `make`'s output.

    Args:
        code_base: Path to the code repository
        work_dir: This release's working directory (logs live under a tmp/
            subdirectory of it)

    Returns:
        Path to the log file the running Luigi pipeline is writing to.
    """
    pipeline_dir = code_base / "pipeline"

    log_dir = work_dir / "tmp"
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pipeline_log = log_dir / f"luigi_run_{timestamp}.log"
    build_release_log = log_dir / f"build_release_{timestamp}.log"

    print(f"\nKicking off pipeline! (spawned in the background, not waiting)")
    cmd = [
        "make",
        f"PIPELINE_LOG={pipeline_log}",
        f"BUILD_RELEASE_LOG={build_release_log}",
        "build-release",
    ]
    print(f"Running: {' '.join(cmd)}")
    print(f"build-release log: {build_release_log}")
    # Both logs above are opened by the Makefile's own recipes, not here --
    # nothing meaningful is left on this process's stdout/stderr to capture.
    subprocess.Popen(
        cmd,
        cwd=pipeline_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        stdin=subprocess.DEVNULL,
        start_new_session=True,
    )

    return pipeline_log


def main() -> int:
    """Main entry point."""
    # The branch to build always comes from the checkout this script is
    # itself running from -- no CLI override, so there's no way to
    # accidentally deploy a different branch than the one you're standing in.
    git_commit = get_current_branch(Path(__file__).resolve().parent)

    parser = argparse.ArgumentParser(
        description="Generate a new BRCA Exchange data release",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Builds from the current branch of this checkout ('{git_commit}').

Example usage:
  %(prog)s /data/releases /data/credentials
  %(prog)s /data/releases /data/credentials --previous-release-dir /data/previous_releases
  %(prog)s /data/releases /data/credentials gene_config_brca_hbop.txt --previous-release-dir /data/previous_releases
        """
    )

    parser.add_argument(
        "root_dir",
        type=str,
        help="Root directory for the release (working directory will be created here)"
    )

    parser.add_argument(
        "credentials_path",
        type=str,
        help="Path to Luigi credentials configuration file"
    )

    parser.add_argument(
        "--previous-release-dir",
        dest="previous_release_dir",
        type=str,
        default=None,
        help="Directory containing previous release for comparison "
             "(optional; defaults to the working directory itself, i.e. no separate comparison dir)"
    )

    parser.add_argument(
        "gene_config_filename",
        type=str,
        nargs='?',
        default="gene_config_brca_only.txt",
        help="Gene configuration filename (default: gene_config_brca_only.txt)"
    )

    args = parser.parse_args()

    # Resolve all paths to absolute
    root_dir = resolve_path(args.root_dir)
    credentials_path = resolve_path(args.credentials_path)
    previous_release_dir = resolve_path(args.previous_release_dir) if args.previous_release_dir else None

    # Generate data date and working directory. The working directory name
    # doubles as the pipeline's PostgreSQL schema name (see
    # brca_pipeline_cfg.mk.j2's VCF_ASSEMBLY_DB_SCHEMA), so it's built from an
    # underscore-separated date rather than DATA_DATE's hyphenated ISO form --
    # unquoted Postgres identifiers can't contain hyphens.
    data_date = datetime.now().strftime("%Y-%m-%d")
    work_dir = root_dir / f"data_release_{data_date.replace('-', '_')}"

    print(f"=== BRCA Exchange Release Generator ===")
    print(f"Data Date: {data_date}")
    print(f"Working Directory: {work_dir}")
    print(f"Gene Configuration: {args.gene_config_filename}")
    print(f"Git Commit: {git_commit}")
    print(f"Previous Release Dir: {previous_release_dir or '(none -- defaulting to working directory)'}")
    print("=" * 40)

    # Create working directory
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created working directory: {work_dir}")

    # Set up code base
    code_base = work_dir / "code"
    clone_or_update_repo(code_base, git_commit)

    # Prepare template context
    template_path = code_base / "pipeline" / "pipeline_running" / "brca_pipeline_cfg.mk.j2"
    config_path = code_base / "pipeline" / "brca_pipeline_cfg.mk"

    context = {
        "DATA_DATE": data_date,
        "WORK_DIR": str(work_dir),
        "CODE_BASE": str(code_base),
        "CREDENTIALS_PATH": str(credentials_path),
        "GIT_COMMIT": git_commit,
        "GENE_CONFIG_FILENAME": args.gene_config_filename,
    }
    if previous_release_dir is not None:
        context["PREVIOUS_RELEASE_DIR"] = str(previous_release_dir)

    # Generate configuration
    generate_config(template_path, config_path, context)

    # Print usage information
    print("\n" + "=" * 40)
    print("Configuration generated successfully!")
    print("\nYou can issue pipeline commands using:")
    print(f"  make CONFIG_PATH={config_path} [cmd]")
    print("-- or --")
    print(f"  cd {code_base / 'pipeline'} && make [cmd]")
    print("=" * 40)

    # Launch the pipeline in the background. build-release (particularly the
    # actual Luigi run within it) can take hours, so we don't wait for it --
    # just report where to watch it.
    try:
        pipeline_log = spawn_pipeline(code_base, work_dir)
    except OSError as e:
        print(f"\nError: could not launch the pipeline: {e}", file=sys.stderr)
        return 1

    print("\n" + "=" * 40)
    print("Pipeline launched in the background (not waiting for it to finish).")
    print(f"Follow its progress with:\n  tail -f {pipeline_log}")
    print("=" * 40)
    return 0


if __name__ == "__main__":
    sys.exit(main())
