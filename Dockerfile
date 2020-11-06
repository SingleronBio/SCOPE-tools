FROM continuumio/miniconda3

RUN git clone -b feature/visible_params https://github.com/SingleronBio/SCOPE-tools.git /opt/SCOPE-tools
RUN conda install -c conda-forge -c bioconda python star fastqc picard ucsc-gtftogenepred subread samtools cutadapt pysam scipy numpy pandas jinja2 matplotlib click scanpy leidenalg louvain
RUN cd /opt/SCOPE-tools && python3 setup.py install
