FROM python:3.4-jessie

RUN apt-get update

RUN apt-get install -y valgrind

ADD src /code/src
ADD experiment /code/experiment
ADD requirements.txt /code
WORKDIR /code

RUN pip install -r requirements.txt
