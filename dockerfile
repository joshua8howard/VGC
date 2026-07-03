FROM condaforge/miniforge3:latest

# Install bioinformatics tools
RUN conda install -y \
    -c conda-forge \
    -c bioconda \
    python=3.11 \
    fastqc \
    fastq-screen \
    multiqc \
    trim-galore \
    cutadapt \
    star=2.7.11b \
    bowtie2 \
    samtools \
    subread \
    sra-tools \
    && conda clean -afy

WORKDIR /workspace

CMD ["/bin/bash"]