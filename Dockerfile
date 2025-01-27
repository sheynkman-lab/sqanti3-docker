# Initially found here: https://github.com/mlorthiois/bioinformatics-docker/tree/e9141202b407ea420a4213840a425f7031a7cd47/SQANTI3
FROM continuumio/miniconda3:4.8.2
LABEL description="Basic docker image containing util software requirements for the sheynkman-lab/Long-Read-Proteogenomics pipeline's SQANTI3 module"
ENV CONDA_ENV="sqanti3"



# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
RUN apt-get update && \
    apt-get install -y \
        procps \
        libcurl4-openssl-dev \
        r-base \
        build-essential && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ${CONDA_ENV} > ${CONDA_ENV}.yml

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${CONDA_ENV}/bin:$PATH

# Clone SQANTI v1.3 from release tag using the --branch flag into new folder /opt/sqanti3
RUN git clone https://github.com/ConesaLab/SQANTI3.git --branch v3.3 /opt/sqanti3 && \
    rm -rf /opt/sqanti3/.git && \
    rm -rf /opt/sqanti3/example

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P /opt/sqanti3/utilities/ && \
    chmod +x /opt/sqanti3/utilities/gtfToGenePred
# Clone cDNA_Cupcake v15.1.0 from release tag using the --branch flag into new folder /opt/cDNA_Cupcake
RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git --branch v15.1.0 /opt/cDNA_Cupcake && \
    cd /opt/cDNA_Cupcake && \
    python setup.py build && \
    python setup.py install && \
    rm -rf /opt/cDNA_Cupcake/.git && \
    rm -rf /opt/cDNA_Cupcake/.git

ENV PYTHONPATH /opt/cDNA_Cupcake/sequence/
ENV PATH "$PATH:/opt/sqanti3/"


# Install R packages necessary to run
RUN Rscript -e "install.packages('BiocManager', repos='http://cran.rstudio.com'); \
        BiocManager::install( \
            c( \
                'ggpubr', \
                'ggthemes', \
                'kableExtra', \
                'NOISeq', \
                'optparse', \
                'paletteer', \
                'qs', \
                'rstatix', \
                'scales', \
                'stringr', \
                'targets', \
                'tarchetypes', \
                'tidyverse' \
            ) \
        )"
RUN Rscript -e "install.packages( \
                    c( \
                    'ggplot2', \
                    'scales', \
                    'reshape', \
                    'gridExtra', \
                    'grid', \
                    'dplyr', \
                    'NOISeq', \
                    'ggplotify'), repo='http://cran.rstudio.com/')"


# Copy additional scripts from bin and add to PATH

COPY src/*.py /opt/sqanti3
RUN chmod +x /opt/sqanti3/*
ENV PATH="$PATH:/opt/sqanti3/"
