安装
=======================================

.. toctree::
   :maxdepth: 2
   :caption: 目录


#. 要求

    * 32G内存或者更高 (根据参考基因组大小确定)
    * Linux
    * conda环境管理器


#. 使用conda安装

    * 创建环境

        .. code-block:: bash

           conda create -n scope
    * 激活环境

        .. code-block:: bash

           conda activate scope
    * 安装scopetools

        .. code-block:: bash

           conda install -c singleronbio scope-tools


    * 验证安装完成

        .. code-block:: bash

            scope -h

            Usage: scope [OPTIONS] COMMAND1 [ARGS]... [COMMAND2 [ARGS]...]...

              Single Cell Omics Preparation Entity Tools

            Options:
              --version   Show the version and exit.
              -h, --help  Show this message and exit.

            Commands:
              STAR           STAR short help
              barcode        extract barcode and umi short help
              cluster        cluster short help
              count          count short help
              cutadapt       cutadapt short help
              featureCounts  featureCounts short help
              run            run short help

#. 准备参照基因组

    * Homo sapiens

        .. code-block:: bash

            wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
            wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

            mkdir -p references/Homo_sapiens/Ensembl/GRCh38
            gzip -c -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.fa
            gzip -c -d Homo_sapiens.GRCh38.99.gtf.gz > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf

            gtfToGenePred -genePredExt -geneNameAsName2 references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf /dev/stdout | \
                awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.refFlat

            conda activate SCOPE-tools
            STAR \
                --runMode genomeGenerate \
                --runThreadN 6 \
                --genomeDir references/Homo_sapiens/Ensembl/GRCh38 \
                --genomeFastaFiles references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.fa \
                --sjdbGTFfile references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf \
                --sjdbOverhang 100


    * Mus musculus

        .. code-block:: bash

            wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
            wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

            mkdir -p references/Mus_musculus/Ensembl/GRCm38
            gzip -c -d Mus_musculus.GRCm38.dna.primary_assembly.fa.gz > references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.fa
            gzip -c -d Mus_musculus.GRCm38.99.gtf.gz > references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.gtf

            gtfToGenePred -genePredExt -geneNameAsName2 references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.gtf /dev/stdout | \
                awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.refFlat

            conda activate SCOPE-tools
            STAR \
                --runMode genomeGenerate \
                --runThreadN 6 \
                --genomeDir references/Mus_musculus/Ensembl/GRCm38 \
                --genomeFastaFiles references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.fa \
                --sjdbGTFfile references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.gtf \
                --sjdbOverhang 100

