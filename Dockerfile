FROM continuumio/miniconda:4.5.4

RUN apt-get update

RUN apt-get install -y valgrind

ADD src /code/src
ADD experiment /code/experiment
ADD env.yml /code

RUN conda env create -f /code/env.yml
