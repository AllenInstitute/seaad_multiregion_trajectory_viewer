# hash:sha256:f92761f3765d856c2dcc89c64e76f6579b1f0837696cdc5a4e670b993c1fb536
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.12.4-R4.3.3-JupyterLab4.2.4-RStudiorstudio-server-2024.04.2-764-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

ENV NODE_PATH=/opt/conda/lib/node_modules

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/shiny-server-jammy/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=2.9.2.1-3ubuntu2 \
        shiny-server=1.5.23.1030 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        libcurl==8.18.0 \
        libgcc-ng==16.1.0 \
        libstdcxx-ng==16.1.0 \
        nodejs==25.2.1 \
        r-arrow==21.0.0 \
    && conda clean -ya

RUN --mount=type=secret,id=git-access-token \
    GIT_ACCESS_TOKEN="$(cat /run/secrets/git-access-token)" \
    pip install -U --no-cache-dir \
    anndata==0.13.2 \
    h5py==3.16.0 \
    numpy==2.5.2 \
    pandas==3.0.5 \
    pyarrow==25.0.1 \
    pybiomart==0.2.0 \
    scipy==1.18.0 \
    -e git+https://github.com/AllenInstitute/abc_atlas_access.git@b16ee223d3466a57b30066973413ccc948c2e983#egg=abc-atlas-access

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("shiny", "1.14.0")' \
    && Rscript -e 'remotes::install_version("bslib", "0.12.0")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.2.1")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.2")' \
    && Rscript -e 'remotes::install_version("scales", "1.4.0")' \
    && Rscript -e 'remotes::install_version("shinycssloaders", "1.1.0")' \
    && Rscript -e 'remotes::install_version("cachem", "1.1.0")' \
    && Rscript -e 'remotes::install_version("DT", "0.34.0")' \
    && Rscript -e 'remotes::install_version("renv", "1.2.4")' \
    && Rscript -e 'remotes::install_version("rsconnect", "1.10.1")' \
    && Rscript -e 'remotes::install_version("curl", "8.0.0")'

ADD "https://github.com/coder/code-server/releases/download/v4.100.3/code-server-4.100.3-linux-amd64.tar.gz" /.code-server/code-server.tar.gz

RUN cd /.code-server \
    && tar -xvf code-server.tar.gz \
    && rm code-server.tar.gz \
    && ln -s /.code-server/code-server-4.100.3-linux-amd64/bin/code-server /usr/bin/code-server \
    && mkdir -p /.vscode/extensions

RUN code-server --extensions-dir="/.vscode/extensions" --install-extension ms-python.python@2025.14.0 \
    && code-server --extensions-dir="/.vscode/extensions" --install-extension ms-toolsai.jupyter@2025.4.1 \
    && code-server --extensions-dir="/.vscode/extensions" --install-extension reditorsupport.r@2.8.6 \
    && code-server --extensions-dir="/.vscode/extensions" --install-extension saoudrizwan.claude-dev@3.35.1

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
