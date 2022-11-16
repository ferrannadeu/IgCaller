FROM debian:buster-slim

ENV DEBIAN_FRONTEND noninteractive

RUN set -e \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
      ca-certificates curl python3 samtools python3-pip git 

RUN set -e \
      && apt-get -y update \
      && apt-get install -y python3-setuptools \
      && pip3 install statistics regex argparse==1.1 numpy==1.16.3 scipy==1.3.0 \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*


RUN set -e \
      && cd / \
      && git clone https://github.com/ferrannadeu/IgCaller \
      && chmod +x IgCaller/IgCaller 

ENV PATH=/IgCaller/:${PATH}

ENTRYPOINT ["IgCaller"]
