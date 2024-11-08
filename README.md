# Ribo-seq Metadata Standardization Pipeline
A standardization pipeline for RiboCrypt and Riboseq.org metadata processing

## Overview
This pipeline standardizes and enriches Ribo-seq metadata from SRA (Sequence Read Archive) for use in RiboCrypt and Riboseq.org. It performs three main functions:

1. **Column Name Standardization**: Harmonizes variant column names into standard formats
   - Example: `CELL_LINE`, `CELL LINE`, `celllines` → `CELL_LINE`

2. **Value Standardization**: Normalizes different representations of the same value
   - Example: `Ribo-seq`, `Riboseq`, `RIBOSEQ` → `Ribo-seq`

3. **Metadata Enrichment**: Adds standardized annotations 
   - Example: Cell line sex annotation (HeLa → female, HEK → male)

## Prerequisites

This was tested and works using:
- Docker version 27.3.1

## Installation

### Using Docker (Recommended)
```bash
# Build the Docker image
docker build -t riboseq-collection .

# Run the container
docker run --rm \
    -v $(pwd)/resources:/usr/local/share/riboseq/resources \
    -v $(pwd):/data \
    riboseq-collection Rscript \
    /usr/local/bin/metadata_main_script.R \
    --output /data/output \
    --temp /data/temp_files \
    --sra /data/SraRunInfo \
    --scripts /data
```