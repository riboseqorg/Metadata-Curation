# Dockerfile
FROM rocker/tidyverse:latest

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and Bioconductor packages
RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")' && \
    R -e 'BiocManager::install(version = "3.20", ask = FALSE)' && \
    R -e 'BiocManager::install(c("S4Vectors", "AnnotationDbi", "Biostrings", "BSgenome", "KEGGREST", "ORFik"), ask = FALSE)' && \
    R -e 'if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")'

RUN R -e 'devtools::install_github("m-swirski/RiboCrypt")'
RUN R -e 'devtools::install_github("rc-biotech/massive_NGS_pipe", upgrade = "always")'

# Install other required packages
RUN R -e 'install.packages(c("data.table", "googlesheets4", "readr", "dplyr", "optparse", "argparse"))'

# Copy scripts to /usr/local/bin
COPY *.R /usr/local/bin/
COPY *.py /usr/local/bin/

# Make scripts executable
RUN chmod +x /usr/local/bin/*.R /usr/local/bin/*.py