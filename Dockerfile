FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y && apt-get install -y make rsync git gcc g++ bzip2 hdf5-tools unzip gfortran curl software-properties-common python3
WORKDIR /test
RUN mkdir -p /opt/julia-1.8.5 && \
    curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz | tar -C /opt/julia-1.8.5 -x -z --strip-components=1 -f -

RUN add-apt-repository ppa:fenics-packages/fenics -y 
RUN apt-get update -q 
RUN apt-get install fenics -y
#link to python3
ENV PYTHON /usr/bin/python3
RUN apt-get install libpython3.8 -y
RUN apt-get install python3-pip -y
RUN pip3 install matplotlib
ADD setup.jl . 
RUN echo "\nPATH=/opt/julia-1.8.5/bin:\$PATH\n" >> /root/.bashrc
RUN ln -s /opt/julia-1.8.5/bin/julia /usr/local/bin/
