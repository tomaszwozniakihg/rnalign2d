# TO USE THIS DOCKERFILE firstly download MUSCLE and VIENNA RNA package
# https://www.tbi.univie.ac.at/RNA/
# The MUSCLE is available here:
# http://www.drive5.com/muscle/
# Please adjust version numbers as found suitable

# To run this dockerfile:
# docker build -t trnalign2d .
# docker run -it trnalign2d /bin/bash

FROM python:3.6.9-stretch
ENV PYTHONUNBUFFERED 1

COPY ./muscle3.8.31_i86linux64 /muscle/muscle
ENV PATH "${PATH}:/muscle"

WORKDIR /
RUN mkdir /vienna
RUN apt update && apt install -y libgsl2
COPY ./ViennaRNA-2.4.14.tar.gz /vienna/ViennaRNA-2.4.14.tar.gz
WORKDIR /vienna
RUN tar -xvzf ViennaRNA-2.4.14.tar.gz
WORKDIR /vienna/ViennaRNA-2.4.14
RUN yes '' | ./configure --without-perl
RUN yes '' | make
RUN yes '' | make install

COPY ./rnalign2d/ /rnalign2d/
WORKDIR /rnalign2d
RUN yes '' | python setup.py install
RUN pip3 install --upgrade pip && pip install -r requirements.txt

RUN mkdir /data
WORKDIR /data
