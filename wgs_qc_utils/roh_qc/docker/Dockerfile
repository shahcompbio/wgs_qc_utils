FROM python:3.8

# docker build -t shahcompbio/roh_qc:0.0.2 .
# docker push shahcompbio/roh_qc:0.0.2

RUN apt-get update
RUN apt-get install git build-essential -y
RUN apt-get install samtools bcftools -y
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt
