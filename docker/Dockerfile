FROM nvcr.io/nvidia/pytorch:24.05-py3

COPY requirements.txt  .
RUN  pip3 install -r requirements.txt 
RUN  apt-get update && apt-get install -y swig4.0
