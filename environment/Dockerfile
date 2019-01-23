FROM registry.codeocean.com/codeocean/r-base:3.4.4-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'devtools::install_version("VennDiagram", \
      version = "1.6.20", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("cowplot", \
      version = "0.9.4", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("dplyr", \
      version = "0.7.4", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("ggplot2", \
      version = "3.1.0", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("ggpubr", \
      version = "0.2", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("ggsci", \
      version = "2.9", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("ggthemes", \
      version = "4.0.1", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("grid", \
      version = "3.4.4", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("papeR", \
      version = "1.0-4", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("plyr", \
      version = "1.8.4", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("survival", \
      version = "2.41-3", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("survminer", \
      version = "0.4.3", \
      dependencies = TRUE)'
RUN Rscript -e 'devtools::install_version("tidyr", \
      version = "0.8.0", \
      dependencies = TRUE)'
