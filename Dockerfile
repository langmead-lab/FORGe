FROM python:3.4-jessie

RUN apt-get update

# for bx-python
RUN apt-get install -y liblzo2-2
RUN apt-get install -y liblzo2-dev

RUN apt-get install -y valgrind

ADD src /code/src
ADD experiment /code/experiment
ADD requirements.txt /code
WORKDIR /code

RUN pip install -r requirements.txt
