FROM python:3.8-slim

WORKDIR /app

# install osprey and biopython
RUN  apt-get update \
  && apt-get install -y wget unzip \
  && wget -O osprey.zip https://github.com/donaldlab/OSPREY3/releases/download/3.2.213/osprey-linux-python3-3.2.zip \
  && unzip osprey.zip \
  && pip install  --no-cache-dir biopython==1.78 osprey --find-link=osprey-linux-python3-3.2/wheelhouse --pre \
  && rm -rf osprey-linux-python3-3.2 osprey.zip

