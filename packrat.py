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
      md5: ""  # Will skip verification if empty
      has_chr_prefix: true
    annotation:
      gencode:
        name: "GENCODE v47"
        url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
        md5: ""
        has_chr_prefix: true

  hg19:
    name: "Human GRCh37/hg19"
    fasta:
      url: "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
      md5: ""
      has_chr_prefix: true
    annotation:
      gencode:
        name: "GENCODE v47lift37"
        url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gtf.gz"
        md5: ""
        has_chr_prefix: true

  mm10:
    name: "Mouse GRCm38/mm10"
    fasta:
      url: "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
      md5: ""
      has_chr_prefix: true
    annotation:
      gencode:
        name: "GENCODE vM35"
        url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz"
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


def verify_gtf_file(gtf_path: Path) -> bool:
    """
    Perform basic sanity checks on a GTF file.

    Args:
        gtf_path: Path to the GTF file

    Returns:
        True if file appears valid, False otherwise
    """
    try:
        with open(gtf_path, "r") as gtf_file:
            # Skip comments
            for line in gtf_file:
                if line.startswith("#"):
                    continue

                # Check first data line
                fields = line.strip().split("\t")
                if len(fields) != 9:
                    console.print(f"[red]Error: GTF file has invalid format (expected 9 fields, got {len(fields)})[/red]")
                    return False

                console.print(f"[green]✓ GTF validation passed[/green]")
                console.print(f"  First feature: {fields[2]} on {fields[0]}")
                return True

            console.print("[red]Error: GTF file appears to be empty[/red]")
            return False
    except Exception as error:
        console.print(f"[red]Error validating GTF: {error}[/red]")
        return False


def standardize_chromosome_names(input_gtf_path: Path, output_gtf_path: Path, add_chr_prefix: bool = True) -> bool:
    """
    Standardize chromosome names in a GTF file by adding/removing 'chr' prefix.

    Args:
        input_gtf_path: Path to the input GTF file
        output_gtf_path: Path to the output GTF file
        add_chr_prefix: If True, add 'chr' prefix; if False, remove it

    Returns:
        True if successful, False otherwise
    """
    try:
        console.print(f"[cyan]Standardizing chromosome names (add_chr_prefix={add_chr_prefix})...[/cyan]")

        modifications_made = False
        with open(input_gtf_path, "r") as input_file, open(output_gtf_path, "w") as output_file:
            for line in input_file:
                if line.startswith("#"):
                    output_file.write(line)
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    output_file.write(line)
                    continue

                chromosome_name = fields[0]

                if add_chr_prefix:
                    # Add chr prefix if not present
                    if not chromosome_name.startswith("chr"):
                        fields[0] = f"chr{chromosome_name}"
                        modifications_made = True
                else:
                    # Remove chr prefix if present
                    if chromosome_name.startswith("chr"):
                        fields[0] = chromosome_name[3:]
                        modifications_made = True

                output_file.write("\t".join(fields) + "\n")

        if modifications_made:
            console.print(f"[green]✓ Chromosome names standardized[/green]")
        else:
            console.print(f"[yellow]Chromosome names already in correct format[/yellow]")

        return True
    except Exception as error:
        console.print(f"[red]Error standardizing chromosome names: {error}[/red]")
        return False


def sort_gtf_file(input_gtf_path: Path, output_gtf_path: Path) -> bool:
    """
    Sort a GTF file by chromosome and position (required for tabix indexing).

    Args:
        input_gtf_path: Path to the input GTF file
        output_gtf_path: Path to the output sorted GTF file

    Returns:
        True if sorting successful, False otherwise
    """
    try:
        console.print(f"[cyan]Sorting GTF file for tabix indexing...[/cyan]")

        # Separate header and data lines
        header_lines = []
        data_lines = []

        with open(input_gtf_path, "r") as input_file:
            for line in input_file:
                if line.startswith("#"):
                    header_lines.append(line)
                else:
                    data_lines.append(line)

        # Sort data lines by chromosome (col 1) and start position (col 4)
        # Using a temporary file for large datasets
        temp_data_file = input_gtf_path.parent / f"{input_gtf_path.name}.temp_data"

        with open(temp_data_file, "w") as temp_file:
            for line in data_lines:
                temp_file.write(line)

        # Use system sort command for efficiency
        # -t$'\t': tab delimiter
        # -k1,1: sort by chromosome (column 1)
        # -k4,4n: sort by start position (column 4, numeric)
        # -k5,5n: sort by end position (column 5, numeric)
        result = subprocess.run(
            ["sort", "-t", "\t", "-k1,1", "-k4,4n", "-k5,5n", str(temp_data_file)],
            capture_output=True,
            text=True,
            check=True,
        )

        # Write sorted output with headers first
        with open(output_gtf_path, "w") as output_file:
            # Write headers
            for header_line in header_lines:
                output_file.write(header_line)

            # Write sorted data
            output_file.write(result.stdout)

        # Clean up temp file
        temp_data_file.unlink()

        console.print(f"[green]✓ GTF file sorted[/green]")
        return True
    except subprocess.CalledProcessError as error:
        console.print(f"[red]Error sorting GTF: {error.stderr}[/red]")
        return False
    except Exception as error:
        console.print(f"[red]Error sorting GTF file: {error}[/red]")
        return False


def compress_with_bgzip(input_file_path: Path, output_file_path: Path) -> bool:
    """
    Compress a file using bgzip (block gzip for tabix indexing).

    Args:
        input_file_path: Path to the input file
        output_file_path: Path to the output .gz file

    Returns:
        True if compression successful, False otherwise
    """
    if not check_tool_available("bgzip"):
        console.print("[red]Error: bgzip not found in PATH[/red]")
        console.print("[yellow]bgzip is part of htslib/tabix. Please install it.[/yellow]")
        return False

    try:
        console.print(f"[cyan]Compressing with bgzip...[/cyan]")

        result = subprocess.run(
            ["bgzip", "-c", str(input_file_path)],
            stdout=open(output_file_path, "wb"),
            stderr=subprocess.PIPE,
            check=True,
        )

        console.print(f"[green]✓ Compressed to {output_file_path}[/green]")
        return True
    except subprocess.CalledProcessError as error:
        console.print(f"[red]Error compressing with bgzip: {error.stderr.decode()}[/red]")
        return False


def index_gtf_with_tabix(gtf_gz_path: Path) -> bool:
    """
    Create tabix index for a bgzipped GTF file.

    Args:
        gtf_gz_path: Path to the bgzipped GTF file

    Returns:
        True if indexing successful, False otherwise
    """
    if not check_tool_available("tabix"):
        console.print("[red]Error: tabix not found in PATH[/red]")
        console.print("[yellow]tabix is part of htslib. Please install it.[/yellow]")
        return False

    try:
        console.print(f"[cyan]Indexing with tabix...[/cyan]")

        result = subprocess.run(
            ["tabix", "-p", "gff", str(gtf_gz_path)],
            capture_output=True,
            text=True,
            check=True,
        )

        index_path = Path(f"{gtf_gz_path}.tbi")
        if index_path.exists():
            console.print(f"[green]✓ Created index: {index_path}[/green]")
            return True
        else:
            console.print("[red]Error: Index file was not created[/red]")
            return False
    except subprocess.CalledProcessError as error:
        console.print(f"[red]Error indexing with tabix: {error.stderr}[/red]")
        return False


def download_genome_annotation(
    genome_id: str,
    annotation_source: str,
    config: dict[str, Any],
    base_output_dir: Path,
) -> bool:
    """
    Download and prepare a genome annotation file.

    Args:
        genome_id: Genome identifier (e.g., 'hg38')
        annotation_source: Annotation source (e.g., 'gencode')
        config: Configuration dictionary
        base_output_dir: Base directory for output files

    Returns:
        True if successful, False otherwise
    """
    if genome_id not in config["genomes"]:
        console.print(f"[red]Error: Unknown genome '{genome_id}'[/red]")
        return False

    genome_config = config["genomes"][genome_id]

    if "annotation" not in genome_config:
        console.print(f"[red]Error: No annotations configured for '{genome_id}'[/red]")
        return False

    if annotation_source not in genome_config["annotation"]:
        console.print(f"[red]Error: Unknown annotation source '{annotation_source}' for '{genome_id}'[/red]")
        available_sources = ", ".join(genome_config["annotation"].keys())
        console.print(f"[yellow]Available sources: {available_sources}[/yellow]")
        return False

    genome_name = genome_config["name"]
    annotation_config = genome_config["annotation"][annotation_source]
    annotation_name = annotation_config["name"]

    console.print(f"\n[bold cyan]Processing {annotation_name} for {genome_name}[/bold cyan]")

    # Setup paths
    annotation_output_dir = base_output_dir / genome_id / "annotation" / annotation_source
    annotation_output_dir.mkdir(parents=True, exist_ok=True)

    compressed_download_path = annotation_output_dir / f"{annotation_source}.gtf.gz"
    uncompressed_gtf_path = annotation_output_dir / f"{annotation_source}.temp.gtf"
    standardized_gtf_path = annotation_output_dir / "genes.gtf"
    final_compressed_path = annotation_output_dir / "genes.gtf.gz"

    # Check if already exists
    if final_compressed_path.exists() and Path(f"{final_compressed_path}.tbi").exists():
        console.print(f"[yellow]Annotation already exists at {final_compressed_path}[/yellow]")

        response = console.input("Overwrite? [y/N]: ").strip().lower()
        if response != "y":
            console.print("[yellow]Skipping download[/yellow]")
            return True

    # Download
    console.print(f"[cyan]Downloading from: {annotation_config['url']}[/cyan]")

    if not download_file(
        annotation_config["url"],
        compressed_download_path,
        f"Downloading {annotation_source}",
    ):
        return False

    # Verify MD5 if provided
    expected_md5 = annotation_config.get("md5", "").strip()
    if expected_md5:
        console.print("[cyan]Verifying MD5 checksum...[/cyan]")
        calculated_md5 = calculate_md5(compressed_download_path)

        if calculated_md5 != expected_md5:
            console.print(f"[red]MD5 mismatch![/red]")
            console.print(f"  Expected: {expected_md5}")
            console.print(f"  Got:      {calculated_md5}")
            return False
        console.print("[green]✓ MD5 checksum verified[/green]")
    else:
        console.print("[yellow]Skipping MD5 verification (no checksum provided)[/yellow]")

    # Decompress
    if not decompress_gzip(compressed_download_path, uncompressed_gtf_path):
        return False

    # Clean up downloaded compressed file
    console.print(f"[cyan]Removing downloaded compressed file...[/cyan]")
    compressed_download_path.unlink()

    # Verify GTF
    if not verify_gtf_file(uncompressed_gtf_path):
        return False

    # Standardize chromosome names if needed
    needs_chr_prefix = not annotation_config.get("has_chr_prefix", True)

    if needs_chr_prefix:
        if not standardize_chromosome_names(
            uncompressed_gtf_path,
            standardized_gtf_path,
            add_chr_prefix=True,
        ):
            return False
        # Clean up temp file
        uncompressed_gtf_path.unlink()
    else:
        # Just rename if no standardization needed
        uncompressed_gtf_path.rename(standardized_gtf_path)

    # Sort GTF file (required for tabix indexing)
    sorted_gtf_path = annotation_output_dir / "genes.sorted.gtf"
    if not sort_gtf_file(standardized_gtf_path, sorted_gtf_path):
        return False

    # Clean up unsorted file
    console.print(f"[cyan]Removing unsorted file...[/cyan]")
    standardized_gtf_path.unlink()

    # Compress with bgzip for tabix indexing
    if not compress_with_bgzip(sorted_gtf_path, final_compressed_path):
        return False

    # Clean up sorted uncompressed file
    console.print(f"[cyan]Removing uncompressed sorted file...[/cyan]")
    sorted_gtf_path.unlink()

    # Create tabix index
    if not index_gtf_with_tabix(final_compressed_path):
        return False

    console.print(f"[bold green]✓ Successfully prepared {annotation_name}[/bold green]")
    console.print(f"  GTF: {final_compressed_path}")
    console.print(f"  Index: {final_compressed_path}.tbi")

    return True


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
    expected_md5 = fasta_config.get("md5", "").strip()
    if expected_md5:
        console.print("[cyan]Verifying MD5 checksum...[/cyan]")
        calculated_md5 = calculate_md5(compressed_fasta_path)

        if calculated_md5 != expected_md5:
            console.print(f"[red]MD5 mismatch![/red]")
            console.print(f"  Expected: {expected_md5}")
            console.print(f"  Got:      {calculated_md5}")
            return False
        console.print("[green]✓ MD5 checksum verified[/green]")
    else:
        console.print("[yellow]Skipping MD5 verification (no checksum provided)[/yellow]")

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
        "--fasta",
        action="store_true",
        help="Download genome FASTA file (default if no other options specified)",
    )

    parser.add_argument(
        "--annotation",
        action="store_true",
        help="Download genome annotation file",
    )

    parser.add_argument(
        "--annotation-source",
        type=str,
        default="gencode",
        help="Annotation source to download (default: gencode)",
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
            if "annotation" in genome_info:
                annotation_sources = ", ".join(genome_info["annotation"].keys())
                console.print(f"    Annotations: {annotation_sources}")
        return 0

    # Check if genome specified
    if not args.genome:
        console.print("[yellow]No genome specified. Use --genome <genome_id> or --list[/yellow]")
        parser.print_help()
        return 1

    # Determine what to download
    download_fasta_flag = args.fasta
    download_annotation_flag = args.annotation

    # If neither specified, default to downloading FASTA only
    if not download_fasta_flag and not download_annotation_flag:
        download_fasta_flag = True

    # Track overall success
    all_successful = True

    # Download FASTA
    if download_fasta_flag:
        # Check for required tools
        if not check_tool_available("samtools"):
            console.print("[red]Error: samtools is required but not found in PATH[/red]")
            console.print("[yellow]Please install samtools and ensure it's in your PATH[/yellow]")
            return 1

        samtools_version = get_tool_version("samtools")
        if samtools_version:
            console.print(f"[green]✓ Found samtools: {samtools_version}[/green]")

        success = download_genome_fasta(args.genome, config, args.output_dir)
        if not success:
            all_successful = False

    # Download annotation
    if download_annotation_flag:
        # Check for required tools
        tools_required = ["bgzip", "tabix"]
        missing_tools = []

        for tool in tools_required:
            if not check_tool_available(tool):
                missing_tools.append(tool)

        if missing_tools:
            console.print(f"[red]Error: Required tools not found in PATH: {', '.join(missing_tools)}[/red]")
            console.print("[yellow]Please install htslib (provides bgzip and tabix)[/yellow]")
            return 1

        bgzip_version = get_tool_version("bgzip")
        if bgzip_version:
            console.print(f"[green]✓ Found bgzip: {bgzip_version}[/green]")

        tabix_version = get_tool_version("tabix")
        if tabix_version:
            console.print(f"[green]✓ Found tabix: {tabix_version}[/green]")

        success = download_genome_annotation(
            args.genome,
            args.annotation_source,
            config,
            args.output_dir,
        )
        if not success:
            all_successful = False

    return 0 if all_successful else 1


if __name__ == "__main__":
    sys.exit(main())
