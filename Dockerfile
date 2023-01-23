FROM ubuntu:jammy-20220531

WORKDIR /home/

RUN apt -y update && apt install -y \
    cmake \
    g++ \
    wget \
    zip \
    vim

RUN wget https://github.com/swacisko/ALGA/archive/refs/heads/master.zip && unzip master.zip

RUN mkdir ALGA-master/build && cd ALGA-master/build && cmake .. && make
