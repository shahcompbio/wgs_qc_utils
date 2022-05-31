FROM python:3.8

# docker build -t shahcompbio/breakpoint_qc:latest .
# docker push shahcompbio/breakpoint_qc:latest

RUN apt-get update
RUN apt-get install git build-essential -y
COPY requirements.txt /requirements.txt
COPY breakpoint_qc.py /usr/bin/breakpoint_qc.py
RUN chmod +x /usr/bin/breakpoint_qc.py
RUN pip install -r requirements.txt

