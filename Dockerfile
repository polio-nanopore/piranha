# start with an image with conda installed
FROM condaforge/mambaforge AS compile-image

# set working directory
WORKDIR /data

# check for updates
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt install build-essential -y --no-install-recommends && \
  apt-get clean && apt-get autoclean

# copy in piranha
RUN git clone https://github.com/polio-nanopore/piranha.git && \
  cd /data/piranha && \
  mamba install conda -n base -c conda-forge -c defaults && \
  mamba env create -f /data/piranha/environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "piranha", "/bin/bash", "-c"]
RUN mamba install -c conda-forge -n piranha python=3.10 conda-pack && \
  cd /data/piranha && \
  pip install .

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda list && \
  conda-pack -n piranha -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar && /venv/bin/conda-unpack

SHELL ["/bin/bash", "-c"]

RUN conda clean --all &&\
  conda remove --name piranha --all

# build piranha
WORKDIR /data/piranha
RUN source /venv/bin/activate && pip install --user --no-cache-dir . \ 
  && pip uninstall -y tensorflow tensorflow-estimator \
  && mamba install -c anaconda -c defaults tensorflow<2.15.0 tensorflow-estimator

# build image
FROM debian:bullseye-slim AS runtime-image

COPY --from=compile-image /root/.local /root/.local
ENV PATH=/root/.local/bin:$PATH

# Copy /venv from the previous stage:
COPY --from=compile-image /venv /venv

# create directory to mount the basecalled directory and output directory
RUN mkdir -p /data/run_data/basecalled && mkdir -p /data/run_data/output

# check for updates
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt install build-essential -y --no-install-recommends && \
  apt install -y procps && \
  apt-get clean && apt-get autoclean

WORKDIR /data/run_data/analysis

SHELL ["/bin/bash", "-c"]

# to allow streamed log output
ENV PYTHONUNBUFFERED=1
ENV PATH=/venv/bin:$PATH

CMD   source /venv/bin/activate && \
      piranha -b /data/run_data/analysis/barcodes.csv -i /data/run_data/basecalled --outdir /data/run_data/output/piranha_output -t ${THREADS}
