FROM rocker/r-base

# install kubeque deps
RUN apt-get update
RUN apt-get install -y time python3-pip && pip3 install gcloud attrs
RUN mkdir /install
ADD kubeque/kubeque-0.1.tar.gz /install/
RUN cd /install/kubeque-0.1 && python3 setup.py install

# install r deps
COPY setup-r /tmp/setup-r
RUN Rscript /tmp/setup-r
