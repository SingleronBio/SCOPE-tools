使用
=======================================

.. toctree::
   :maxdepth: 2
   :caption: 目录

#. 示例数据


    示例的样本数据存储于 `Open Science Framework`_. 目前已经上传SCOPEv2样本数据, 其余类型样本SCOPEv1, 10X, dropseq, indrop, BD Rhapsody待测试上传.

    .. _Open Science Framework: https://osf.io/792wz/?view_only=f0345853d9b141708be8a533d09ea93c

#. 快速使用

    .. code-block:: bash

        scope run \
            --fq1 \
            --fq2 \
            --outdir


    * 参数说明

        * --fq1：
        * --fq2：

    * 示例报告

    * 输出文件


#. 使用说明

    SCOPE-tools包含7个子命令, 分别是barcode, cutadapt, STAR, featureCounts, count, cluster和run.

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

    * barcode

        基于read1序列信息 (barcode序列, linker序列, 质量值和polyT长度) 过滤, 提取并矫正barcode, 将矫正后的barcode和原始的UMI序列添加到read2的ID中.

        过滤规则:
            #. 过滤polyT碱基数目<10的reads
            #. 过滤barcode和UMI碱基质量值低于14的个数>2
            #. 过滤两段接头中任何一段中错配碱基数>1
            #. 过滤三段barcode中错配碱基数和>1

        barcode矫正规则:
            将未出现在whitelist中的barcode, 矫正为与whitelist中 汉明距离_ 为1的barcode.
                .. _汉明距离: https://en.wikipedia.org/wiki/Hamming_distance

        * 示例
            .. code-block:: bash

                scope barcode \
                    --fq1 ./rawdata/R2005073_L1_1.fq.gz \
                    --fq2 ./rawdata/R2005073_L1_2.fq.gz \
                    --sample scopev2 \
                    --outdir ./scopev2 \
                    --bctype SCOPEv2

        * 参数说明
            * --fq1: read1 fastq文件路径
            * --fq2: read2 fastq文件路径
            * --sample: 样本名称
            * --outdir: 输出路径
            * --bctype: 预置的接头类型
            * --pattern: 自定义的接头结构, 字母C, L, U, T分别表示cell barcode、linker、UMI、T碱基, 数字表示碱基长度. C8L16C8L16C8L1U8T18即表示以下结构：

                CCCCCCCCLLLLLLLLLLLLLLLLCCCCCCCCLLLLLLLLLLLLLLLLCCCCCCCCLUUUUUUUUTTTTTTTTTTTTTTTTTT

                从第一个碱基开始，为C位置碱基构成cell barcode, 为U位置的碱基构成UMI, 第一段为L的碱基为linker1, 第二段为L的碱基为linker2, 末尾T为polyT.

            * --whitelist 自定义的cell barcode白名单文件
            * --linker: 自定义的linker白名单文件
            * --lowQual: 定义为低质量碱基的质量值, 默认: 14
            * --lowNum: 允许的低质量碱基的个数, 默认: 2

    * cutadapt

        调用 cutadapt_ 对read2进行质控.
            * trim接头序列
            * trim两端低质量碱基

            .. _cutadapt: https://en.wikipedia.org/wiki/Hamming_distance

        * 示例
            .. code-block:: bash

                scope cutadapt \
                    --fq ./scopev2/01.barcode/scopev2_2.fq.gz \
                    --sample scopev2 \
                    --outdir ./scopev2 \
                    --adapter p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC polyT=A{18} \
                    --overlap 5

        * 参数说明
            * --fq: barcode处理后的read2 fastq文件
            * --sample: 样本名称
            * --outdir: 输出路径
            * --adapter: 接头序列, 可多次使用, 默认: p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC polyT=A{18}
            * --overlap: 认为检测到接头时重叠碱基数, 默认: 5
            * --minimum-length: 允许的最短序列长度, 默认: 20
            * --nextseq-trim: trim使用的质量值 (忽略G, 针对双色试剂, 如 NextSeq_), 默认: 20
            * --thread: 线程数, 默认: 2

            .. _NextSeq: https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/

    * STAR

        调用 STAR_ 将read2序列定位到基因组上.

        .. _STAR: https://github.com/alexdobin/STAR

        * 示例
            .. code-block:: bash

                scope STAR \
                    --fq ./scopev2/02.cutadapt/scopev2_2.fq.gz \
                    --sample scopev2 \
                    --outdir ./scopev2 \
                    --refFlat ./references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.refFlat \
                    --genomeDir ./references/Homo_sapiens/Ensembl/GRCh38

        * 参数说明
            * --fq: cutadapt处理后的read2 fastq文件
            * --sample: 样本名称
            * --readFilesCommand: STAR读取输入文件的命令, 默认: zcat
            * --genomeDir: 参考基因组路径, 包含STAR所建索引
            * --runThreadN: 线程数, 默认: 2
            * --outdir: 输出路径
            * --refFlat: refFlat文件路径

    * featureCounts

        调用 featureCounts_ 将定位到基因组上的reads，进一步定位到基因上

        .. _featureCounts: http://bioinf.wehi.edu.au/featureCounts/

        * 示例
            .. code-block:: bash

                scope featureCounts \
                    --bam ./scopev2/03.STAR/scopev2_Aligned.sortedByCoord.out.bam \
                    --annot .references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf \
                    --sample scopev2 \
                    --outdir ./scopev2

        * 参数说明
            * --bam: STAR比对并排序后的bam文件
            * --sample: 样本名称
            * --annot: 基因组注释文件, gtf格式
            * --nthreads: 线程数, 默认: 2
            * --format: 输入文件的格式, 默认: BAM
            * --outdir: 输出路径

    * count

        描述

        * 示例

        * 参数说明

    * cluster

        描述

        * 示例

        * 参数说明

    * run

        描述

        * 示例

        * 参数说明
