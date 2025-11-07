#!/usr/bin/env -S uv --quiet run --script
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pyyaml",
#     "requests",
#     "rich",
# ]
# ///

"""
packrat - Reference genome database manager

Like a packrat collecting and organizing treasures,
packrat gathers and curates your genomic references.
"""

import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import requests
import yaml
from rich.console import Console
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

console = Console()

# Embedded YAML configuration for reference data sources
CONFIG_YAML = """
genomes:
  hg38:
    name: "Human GRCh38/hg38"
    fasta:
      url: "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
      md5: "e3f10d2c1e2f2f5b5c5a5c5a5c5a5c5a"  # Placeholder - will skip verification if empty
      has_chr_prefix: true

  hg19:
    name: "Human GRCh37/hg19"
    fasta:
      url: "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
      md5: ""
      has_chr_prefix: true

  mm10:
    name: "Mouse GRCm38/mm10"
    fasta:
      url: "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
      md5: ""
      has_chr_prefix: true
"""


def load_configuration() -> dict[str, Any]:
    """Load and parse the embedded YAML configuration."""
    return yaml.safe_load(CONFIG_YAML)


def check_tool_available(tool_name: str) -> bool:
    """Check if a required tool is available in PATH."""
    return shutil.which(tool_name) is not None


def get_tool_version(tool_name: str) -> str | None:
    """Get the version of a tool if available."""
    try:
        result = subprocess.run(
            [tool_name, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        version_output = result.stdout + result.stderr
        return version_output.strip().split("\n")[0]
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return None


def calculate_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Calculate MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(chunk_size), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def download_file(url: str, output_path: Path, description: str = "Downloading") -> bool:
    """
    Download a file from a URL with progress bar.

    Args:
        url: URL to download from
        output_path: Path to save the downloaded file
        description: Description for the progress bar

    Returns:
        True if download successful, False otherwise
    """
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)

        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))

        with Progress(
            TextColumn("[bold blue]{task.description}"),
            BarColumn(),
            DownloadColumn(),
            TransferSpeedColumn(),
            TimeRemainingColumn(),
            console=console,
        ) as progress:
            task = progress.add_task(description, total=total_size)

            with open(output_path, "wb") as output_file:
                for chunk in response.iter_content(chunk_size=8192):
                    output_file.write(chunk)
                    progress.update(task, advance=len(chunk))

        return True
    except requests.exceptions.RequestException as error:
        console.print(f"[red]Error downloading file: {error}[/red]")
        return False


def decompress_gzip(input_path: Path, output_path: Path) -> bool:
    """
    Decompress a gzipped file.

    Args:
        input_path: Path to the gzipped file
        output_path: Path to save the decompressed file

    Returns:
        True if decompression successful, False otherwise
    """
    try:
        console.print(f"[cyan]Decompressing {input_path.name}...[/cyan]")

        result = subprocess.run(
            ["gunzip", "-c", str(input_path)],
            stdout=open(output_path, "wb"),
            stderr=subprocess.PIPE,
            check=True,
        )

        console.print(f"[green]✓ Decompressed to {output_path}[/green]")
        return True
    except subprocess.CalledProcessError as error:
        console.print(f"[red]Error decompressing file: {error.stderr.decode()}[/red]")
        return False


def index_fasta(fasta_path: Path) -> bool:
    """
    Create FASTA index using samtools faidx.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        True if indexing successful, False otherwise
    """
    if not check_tool_available("samtools"):
        console.print("[red]Error: samtools not found in PATH[/red]")
        return False

    try:
        console.print(f"[cyan]Indexing {fasta_path.name}...[/cyan]")

        result = subprocess.run(
            ["samtools", "faidx", str(fasta_path)],
            capture_output=True,
            text=True,
            check=True,
        )

        index_path = Path(f"{fasta_path}.fai")
        if index_path.exists():
            console.print(f"[green]✓ Created index: {index_path}[/green]")
            return True
        else:
            console.print("[red]Error: Index file was not created[/red]")
            return False
    except subprocess.CalledProcessError as error:
        console.print(f"[red]Error indexing FASTA: {error.stderr}[/red]")
        return False


def verify_fasta_file(fasta_path: Path) -> bool:
    """
    Perform basic sanity checks on a FASTA file.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        True if file appears valid, False otherwise
    """
    try:
        with open(fasta_path, "r") as fasta_file:
            first_line = fasta_file.readline().strip()

            if not first_line.startswith(">"):
                console.print("[red]Error: File does not appear to be a valid FASTA (missing '>' header)[/red]")
                return False

            console.print(f"[green]✓ FASTA validation passed[/green]")
            console.print(f"  First sequence: {first_line[:80]}")
            return True
    except Exception as error:
        console.print(f"[red]Error validating FASTA: {error}[/red]")
        return False


def download_genome_fasta(genome_id: str, config: dict[str, Any], base_output_dir: Path) -> bool:
    """
    Download and prepare a genome FASTA file.

    Args:
        genome_id: Genome identifier (e.g., 'hg38')
        config: Configuration dictionary
        base_output_dir: Base directory for output files

    Returns:
        True if successful, False otherwise
    """
    if genome_id not in config["genomes"]:
        console.print(f"[red]Error: Unknown genome '{genome_id}'[/red]")
        return False

    genome_config = config["genomes"][genome_id]
    genome_name = genome_config["name"]
    fasta_config = genome_config["fasta"]

    console.print(f"\n[bold cyan]Processing {genome_name} ({genome_id})[/bold cyan]")

    # Setup paths
    genome_output_dir = base_output_dir / genome_id / "fasta"
    genome_output_dir.mkdir(parents=True, exist_ok=True)

    compressed_fasta_path = genome_output_dir / f"{genome_id}.fa.gz"
    final_fasta_path = genome_output_dir / "genome.fa"

    # Download
    if final_fasta_path.exists():
        console.print(f"[yellow]FASTA already exists at {final_fasta_path}[/yellow]")

        response = console.input("Overwrite? [y/N]: ").strip().lower()
        if response != "y":
            console.print("[yellow]Skipping download[/yellow]")

            # Still check if index exists
            if not Path(f"{final_fasta_path}.fai").exists():
                console.print("[cyan]Index missing, creating...[/cyan]")
                return index_fasta(final_fasta_path)
            return True

    console.print(f"[cyan]Downloading from: {fasta_config['url']}[/cyan]")

    if not download_file(fasta_config["url"], compressed_fasta_path, f"Downloading {genome_id}"):
        return False

    # Verify MD5 if provided
    if fasta_config.get("md5"):
        console.print("[cyan]Verifying MD5 checksum...[/cyan]")
        calculated_md5 = calculate_md5(compressed_fasta_path)
        expected_md5 = fasta_config["md5"]

        if calculated_md5 != expected_md5:
            console.print(f"[red]MD5 mismatch![/red]")
            console.print(f"  Expected: {expected_md5}")
            console.print(f"  Got:      {calculated_md5}")
            return False
        console.print("[green]✓ MD5 checksum verified[/green]")

    # Decompress
    if not decompress_gzip(compressed_fasta_path, final_fasta_path):
        return False

    # Clean up compressed file
    console.print(f"[cyan]Removing compressed file...[/cyan]")
    compressed_fasta_path.unlink()

    # Verify FASTA
    if not verify_fasta_file(final_fasta_path):
        return False

    # Create index
    if not index_fasta(final_fasta_path):
        return False

    console.print(f"[bold green]✓ Successfully prepared {genome_name}[/bold green]")
    console.print(f"  FASTA: {final_fasta_path}")
    console.print(f"  Index: {final_fasta_path}.fai")

    return True


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="packrat - Reference genome database manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--genome",
        type=str,
        help="Genome to download (e.g., hg38, hg19, mm10)",
    )

    parser.add_argument(
        "--list",
        action="store_true",
        help="List available genomes",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help="Base output directory (default: current directory)",
    )

    args = parser.parse_args()

    # Load configuration
    config = load_configuration()

    # List genomes
    if args.list:
        console.print("\n[bold cyan]Available genomes:[/bold cyan]")
        for genome_id, genome_info in config["genomes"].items():
            console.print(f"  {genome_id:10s} - {genome_info['name']}")
        return 0

    # Check if genome specified
    if not args.genome:
        console.print("[yellow]No genome specified. Use --genome <genome_id> or --list[/yellow]")
        parser.print_help()
        return 1

    # Check for required tools
    if not check_tool_available("samtools"):
        console.print("[red]Error: samtools is required but not found in PATH[/red]")
        console.print("[yellow]Please install samtools and ensure it's in your PATH[/yellow]")
        return 1

    samtools_version = get_tool_version("samtools")
    if samtools_version:
        console.print(f"[green]✓ Found samtools: {samtools_version}[/green]")

    # Download genome
    success = download_genome_fasta(args.genome, config, args.output_dir)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
