FROM continuumio/miniconda3:4.9.2 as base
USER root

RUN conda install pandas=1.3.1 numpy=1.20.3 scipy\
        && apt-get update --allow-releaseinfo-change \
        && apt-get install -y --no-install-recommends procps=2:3.3.15-2 unattended-upgrades=1.11.2 \
        && unattended-upgrade -d -v \
        && apt-get remove -yq unattended-upgrades \
        && apt-get autoremove -yq \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/* \
        && groupadd -r -g 1000 ubuntu \
        && useradd -r -g ubuntu -u 1000 ubuntu

USER ubuntu

WORKDIR /

COPY association_mutation_models.py association_mutation_models.py

#USER ubuntu

FROM base as test
# Testing

COPY --chown=ubuntu:ubuntu test test

WORKDIR /test

RUN bash test.sh

FROM base
