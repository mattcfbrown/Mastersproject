FROM --platform=linux/amd64 python:3.9.2

RUN apt-get update
RUN apt-get -y install apt-transport-https ca-certificates curl gnupg2 software-properties-common
RUN curl -fsSL https://download.docker.com/linux/debian/gpg | apt-key add -
RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"

RUN apt-get update
RUN apt-get -y install docker-ce


RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir geneci==2.0.1 numpy==1.26.4 requests==2.29.0