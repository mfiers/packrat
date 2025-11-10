# packrat

A single-file Python script for automated management of genomic reference data and tool-specific databases.

> Like a packrat collecting and organizing treasures, packrat gathers and curates your genomic references.

<img width="512" height="512" alt="image" src="https://github.com/user-attachments/assets/324552b1-1072-4083-a08b-0767a109660a" />


## Overview

This tool automates the download, indexing, and preparation of genomic reference data for bioinformatics pipelines. It handles:

- **Genome references**: Downloads FASTA files and creates indexes (.fai)
- **Gene annotations**: Downloads and indexes annotation files (GTF format)
- **Tool databases**: Builds indexes for STAR and Salmon (BWA support coming soon)

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
│       ├── genes.gtf       (decompressed for STAR/Salmon)
│       ├── genes.gtf.gz    (bgzipped for tabix)
│       └── genes.gtf.gz.tbi
├── STAR/
│   └── 2.7.11b/
│       └── <index files>
└── salmon/
    └── 1.10.3/
        ├── transcripts.fa
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
  - **Salmon** - Optional, for building transcript indexes (requires gffread)
  - **gffread** - Optional, for extracting transcript sequences (required by Salmon)
  - **BWA** - Coming soon

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

# Download RefSeq annotation (alternative to GENCODE)
./packrat.py --genome mm39 --annotation --annotation-source refseq

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

# Build Salmon transcript index
./packrat.py --genome hg38 --salmon

# Build everything: download FASTA, annotation, and build both STAR and Salmon indexes
./packrat.py --genome hg38 --fasta --annotation --star --salmon

# Merge genomes for hybrid/xenograft analysis
./packrat.py merge hs1mm39 hs1 mm39

# Specify annotation source for merging
./packrat.py merge hs1mm39 hs1 mm39 --annotation-source refseq
```

### Genome Merging

For analyzing hybrid samples (e.g., human xenografts in mouse, PDX models), packrat can merge multiple genomes:

```bash
# First, download the source genomes
./packrat.py --genome hs1 --fasta --annotation
./packrat.py --genome mm39 --fasta --annotation

# Merge them into a hybrid genome
./packrat.py merge hs1mm39 hs1 mm39

# Or specify which annotation source to use (e.g., 'gencode' or 'refseq')
./packrat.py merge hs1mm39 hs1 mm39 --annotation-source refseq

# Build STAR index for the merged genome
./packrat.py --genome hs1mm39 --star
```

**What happens during merge:**
- FASTA chromosomes are prefixed: `chr1` → `hs1__chr1`, `mm39__chr1`
- Gene IDs are prefixed: `ENSG00000...` → `hs1__ENSG00000...`
- Gene names are prefixed: `GAPDH` → `hs1__GAPDH`, `mm39__Gapdh`
- Creates properly sorted and indexed GTF
- Ready for STAR alignment to distinguish human vs mouse reads
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
