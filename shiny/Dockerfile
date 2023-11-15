FROM rocker/shiny:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev

# Install R packages required
RUN R -e 'install.packages("remotes", repos="https://packagemanager.posit.co/cran/__linux__/focal/2023-07-31")'

# Install CRAN packages with specified versions
RUN R -e ' \
    remotes::install_version("shiny", version = "1.7.5"); \
    remotes::install_version("shinyBS", version = "0.61.1"); \
    remotes::install_version("DT", version = "0.28"); \
    remotes::install_version("tidyverse", version = "2.0.0"); \
    remotes::install_version("plotly", version = "4.10.3"); \
'

# Install Bioconductor package
RUN R -e 'install.packages("BiocManager", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-07-31"); \
          BiocManager::install("DESeq2", version = "1.40.2")'

WORKDIR /home/shinyusr
COPY app.R app.R
COPY input_data input_data
CMD Rscript app.R
