## This file builds a container with R 
## for the epigenetics manuscript

FROM rocker/r-ver:latest
LABEL maintainer="atyler"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y libglpk-dev \
	libgmp-dev \
	libssl-dev \
	libxml2-dev \
	pandoc \
	pandoc-citeproc \
	zlib1g-dev

RUN ["install2.r", "abind", "here", "igraph", "RColorBrewer", "rmarkdown", "shape", "VennDiagram", "ape", "gprofiler2", "e1071", "rticles", "tinytex", "knitr", "pheatmap", "hexbin", "MASS", "gridExtra", "ggplotify", "BiocManager", "remotes"]
RUN R -e "BiocManager::install('grimbough/biomaRt')"
RUN R -e "BiocManager::install('DESeq2')"

WORKDIR /chrom_ms/
CMD ["R"]
