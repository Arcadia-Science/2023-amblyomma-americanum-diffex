FROM rocker/shiny:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev

# Install R packages required
RUN R -e 'install.packages(c(\
              "shiny", \
              "shinyBS", \
              "DT", \
              "tidyverse", \
              "plotly"), \
            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-11-13"\
          )'

# Install Bioconductor package
RUN R -e 'install.packages("BiocManager", repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-11-13"); \
          BiocManager::install("DESeq2")'

WORKDIR /home/shinyusr
COPY app.R app.R
COPY input_data input_data
CMD Rscript app.R

