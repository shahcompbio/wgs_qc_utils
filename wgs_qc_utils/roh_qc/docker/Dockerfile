FROM python:3.8

# docker build -t quay.io/soymintc/roh_qc . 
# docker push quay.io/soymintc/roh_qc

RUN apt-get update
RUN apt-get install git build-essential -y
RUN apt-get install samtools bcftools -y
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt
