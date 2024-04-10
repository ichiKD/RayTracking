FROM ubuntu:22.04

RUN apt-get update 
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    cmake \
    build-essential \
    gcc g++ \
    ninja-build \
    imagemagick \
    git \
    gdb \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /cg1

COPY ./scripts/run_all_tests.sh /cg1/scripts/run_all_tests.sh
RUN chmod +x /cg1/scripts/run_all_tests.sh
ENTRYPOINT ["/cg1/scripts/run_all_tests.sh"]


CMD /bin/bash