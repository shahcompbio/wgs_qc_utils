FROM python:3.8

# docker build -t shahlab/classifycopynumber:latest .
# docker push shahlab/classifycopynumber:latest

RUN apt-get update
RUN apt-get install git build-essential -y
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt
