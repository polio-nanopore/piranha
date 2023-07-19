# start with an image with conda installed
FROM mambaorg/micromamba:latest AS compile-image

# set working directory
WORKDIR /data

# check for updates
RUN apt-get update -y && apt-get upgrade -y

# install gcc
RUN apt install build-essential -y --no-install-recommends

# copy in piranha
RUN git clone https://github.com/polio-nanopore/piranha.git

# TEMP change branch
RUN cd /data/piranha

# install env
RUN /opt/conda/bin/mamba env create -f /data/piranha/environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "piranha", "/bin/bash", "-c"]
RUN cd /data/piranha  && pip install .

# Install conda-pack:
RUN conda install -c conda-forge conda-pack

# fix conda conflicts for conda-pack
RUN conda list

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n piranha -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack


SHELL ["/bin/bash", "-c"]
# build piranha
WORKDIR /data/piranha
RUN source /venv/bin/activate && pip install --user --no-cache-dir .
#--user .


# build image
FROM debian:bullseye-slim AS runtime-image

COPY --from=compile-image /root/.local /root/.local

# Copy /venv from the previous stage:
COPY --from=compile-image /venv /venv
#COPY --from=compile-image /opt/conda/ /opt/conda

# create directory to mount the basecalled directory
RUN mkdir -p /data/run_data/basecalled

# create directory to mount the output directory
RUN mkdir -p /data/run_data/output

WORKDIR /data/run_data/analysis

SHELL ["/bin/bash", "-c"]

# to allow streamed log output
ENV PYTHONUNBUFFERED=1

ENV PATH=/root/.local/bin:$PATH

ENTRYPOINT   source /venv/bin/activate
