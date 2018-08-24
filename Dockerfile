FROM continuumio/miniconda3:4.5.4

RUN apt-get update
RUN apt-get install -y valgrind procps

ADD env.yml /code/env.yml
RUN conda update -n base conda
RUN conda env create -f /code/env.yml

ADD experiment /code/experiment
ADD src /code/src
