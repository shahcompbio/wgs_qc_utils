FROM python:3.8

# docker build -t shahlab/roh_qc:latest .
# docker push shahlab/roh_qc:latest

RUN apt-get update
RUN apt-get install git build-essential -y
RUN apt-get install samtools bcftools -y
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt
COPY plot_roh.py /plot_roh.py
ENTRYPOINT ["python","/plot_roh.py"]
