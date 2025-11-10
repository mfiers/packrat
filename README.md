# packrat

A single-file Python script for automated management of genomic reference data and tool-specific databases.

> Like a packrat collecting and organizing treasures, packrat gathers and curates your genomic references.

<img width="512" height="512" alt="image" src="https://github.com/user-attachments/assets/324552b1-1072-4083-a08b-0767a109660a" />


## Overview

This tool automates the download, indexing, and preparation of genomic reference data for bioinformatics pipelines. It handles:

- **Genome references**: Downloads FASTA files and creates indexes (.fai)
- **Gene annotations**: Downloads and indexes annotation files (GTF format)
- **Tool databases**: Builds indexes for STAR (Salmon, BWA support coming soon)

## Features

### Organized File Structure

Files are organized by genome and data type, or by tool and version:

```
hg38/
├── fasta/
│   ├── genome.fa
│   └── genome.fa.fai
├── annotation/
│   └── gencode/
│       ├── genes.gtf       (decompressed for STAR)
│       ├── genes.gtf.gz    (bgzipped for tabix)
│       └── genes.gtf.gz.tbi
└── STAR/
    └── 2.7.11b/
        └── <index files>
```

### Modular Design

- **Extensible**: Easy to add support for additional tools
- **Self-contained**: Python dependencies managed inline with uv script metadata (PEP 723)
- **Environment-aware**: Detects bioinformatics tools (STAR, Salmon, BWA, samtools) from your PATH
- **Version tracking**: Automatically detects and uses tool versions from your environment
- **Configurable**: Download paths defined in an embedded YAML configuration at the top of the script

### Data Management

**Quality assurance:**
- Sanity checks for all downloads to verify file integrity
- Proper indexing for all reference files
- Correct handling of coordinate systems (0-based vs 1-based)

**Standardization:**
- All gene annotations converted to GTF format
- Chromosome names standardized with "chr" prefix (chr1, chr2, etc.)

## Requirements

- Python 3.13+
- [uv](https://github.com/astral-sh/uv) - Python package manager
- Bioinformatics tools must be available in your PATH:
  - **samtools** - For FASTA indexing
  - **bgzip** and **tabix** (from htslib) - For GTF compression and indexing
  - **STAR** - Optional, for building STAR alignment indexes
  - **Salmon, BWA** - Coming soon

**Note**: Python dependencies are managed inline using uv's script metadata format (PEP 723). No separate virtual environment setup needed - uv handles everything automatically when you run the script.

## Usage

The script uses uv's inline dependency management, so just run it directly:

```bash
# Make executable (first time only)
chmod +x packrat.py

# List available genomes and annotations
./packrat.py --list

# Download genome FASTA (default behavior)
./packrat.py --genome hg38

# Download T2T-CHM13 (complete telomere-to-telomere human genome)
./packrat.py --genome hs1

# Download mouse genome (mm39 is the latest)
./packrat.py --genome mm39

# Download genome annotation
./packrat.py --genome hg38 --annotation

# Download both FASTA and annotation
./packrat.py --genome hg38 --fasta --annotation

# Specify custom output directory
./packrat.py --genome hg38 --annotation --output-dir /data/references

# Download specific annotation source
./packrat.py --genome hg38 --annotation --annotation-source gencode

# Force re-download even if files exist
./packrat.py --genome hg38 --annotation --force

# Build STAR index (auto-detects annotation if only one exists)
./packrat.py --genome hg38 --star

# Build everything: download FASTA, annotation, and build STAR index
./packrat.py --genome hg38 --fasta --annotation --star
```

**Note**: By default, packrat will skip downloads/builds if the output files already exist. Use `--force` to re-download or rebuild.

Alternatively, run explicitly with uv:

```bash
uv run packrat.py --genome hg38 --annotation
```

## Configuration

Download paths and reference sources are configured in the YAML section at the top of the script. Edit this section to customize data sources.

## Contributing

To add support for a new tool:
1. Ensure the tool is available in your PATH
2. Add configuration to the YAML section
3. Implement the tool-specific indexing logic

---

**Note**: This README assumes familiarity with bioinformatics workflows and genomic reference data.
