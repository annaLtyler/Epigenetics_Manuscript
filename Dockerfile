## This file builds a container with R 
## and ChromHMM for the epigenetics manuscript

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

RUN ["install2.r", "abind", "here", "igraph", "RColorBrewer", "rmarkdown", "shape"]
RUN R -e "install.packages('qtl2convert', repos='http://cran.r-project.org')"
RUN R -e "install.packages('qtl2', repos = 'http://cran.r-project.org')"

FROM biocontainers/chromhmm

WORKDIR /chrom_ms/
CMD ["R"]
