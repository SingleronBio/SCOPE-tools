FROM continuumio/miniconda3

RUN git clone -b master https://github.com/SingleronBio/SCOPE-tools.git /opt/SCOPE-tools
RUN conda install -c conda-forge -c bioconda python=3.8.6 star=2.7.6a fastqc=0.11.9 picard=2.23.8 ucsc-gtftogenepred=377 subread=2.0.1 samtools=1.11 leidenalg=0.8.2 louvain=0.6.1
RUN cd /opt/SCOPE-tools && python3 setup.py install
