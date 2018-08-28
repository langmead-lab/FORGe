FROM pypy:3-6-jessie

RUN apt-get update -y && apt-get install -y valgrind procps sysstat

ADD requirements.txt /code/requirements.txt
RUN pypy3 -m pip install --no-cache-dir -r /code/requirements.txt

RUN wget https://github.com/refresh-bio/KMC/archive/v3.1.1rc1.tar.gz && \
    tar zxf v3.1.1rc1.tar.gz && \
    cd KMC-3.1.1rc1 && \
    make  -f makefile

RUN cp KMC-3.1.1rc1/bin/kmc_tools \
       KMC-3.1.1rc1/bin/kmc_dump \
       KMC-3.1.1rc1/bin/kmc /usr/local/bin

ADD experiment /code/experiment
ADD src /code/src
