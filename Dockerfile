FROM rocker/shiny:4.4

RUN R -e "install.packages('remotes'); \
    remotes::install_version('bslib', version = '0.9.0'); \
    remotes::install_version('circlize', version = '0.4.16'); \
    remotes::install_version('ggridges', version = '0.5.6'); \
    remotes::install_version('magrittr', version = '2.0.3'); \
    remotes::install_version('patchwork', version = '1.3.0'); \
    remotes::install_version('plotly', version = '4.10.4'); \
    remotes::install_version('shiny', version = '1.10.0'); \
    remotes::install_version('shinyFiles', version = '0.9.3'); \
    remotes::install_version('shinyWidgets', version = '0.9.0'); \
    remotes::install_version('tidyverse', version = '2.0.0'); \
    remotes::install_version('tinytex', version = '0.56.0'); \
    remotes::install_version('scales', version = '1.3.0'); \
    remotes::install_version('fmsb', version = '0.7.6'); \
    remotes::install_version('ggrepel', version = '0.9.6'); \
    remotes::install_version('png', version = '0.1-8'); \
    remotes::install_version('hexbin', version = '1.28.5'); \
    "

EXPOSE 3838