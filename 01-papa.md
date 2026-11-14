# 使用PAPA检测和定量新发现的最后一个外显子/多聚腺苷酸化事件

作者 · 不加酵母的面团（DoughWithoutYeast）

原发表于 2026-08-25 · [博客园原文](https://www.cnblogs.com/DoughWithoutYeast/p/22672503)

原文更新于 2026-10-26

## 0.0-introduction

`PAPA：用于从批量RNA测序数据中检测和定量新型末端外显子/多聚腺苷化事件的Snakemake流程`这是在github内对于该工具的介绍。\
2020年该方法开发并发布在预印版内。`https://www.biorxiv.org/content/10.1101/2020.08.21.261644v2.full.pdf+html`\
2025年作者使用该方法在Nature Neuroscience发表了相关文献。`https://www.nature.com/articles/s41593-025-02050-w`\
该方法从StringTie组装的转录本中提取最后一个外显子，并根据其与3’端测序所得poly(A)位点的距离进行筛选。可以对全新末端外显子和多聚腺苷酸化事件的检测。

* 使用StringTie从比对后的RNA-seq读段中组装转录本
* 过滤组装的转录本，保留在相同条件下各样本中平均表达量达到最低阈值的转录本
* 筛选新发现的末外显子，这些外显子可延伸已知外显子（可选择性检查5'端是否匹配），或包含新的末剪接位点（可选择性检查新的5'剪接位点是否与已知5'剪接位点匹配）
* 过滤组装的转录本，保留末外显子附近存在参考多聚A位点（如PolyASite）或保守多聚A信号基序的转录本
* 将新发现的末外显子异构体与参考注释合并，并使用Salmon对转录本进行定量
* 通过tximport输出计数/TPM矩阵，供下游差异转录本使用分析软件使用
* 使用DEXSeq对两种条件之间的转录本使用情况进行差异分析\
  ![image](3805523-20260825093358432-869008610.png)

## 0-installation

### 0.1-download from git

`https://github.com/frattalab/PAPA`在PAPA内下载对应的zip文件，解压缩。可以看到`papa-main`目录：\
![image](3805523-20260825093717079-1256115372.png)\
**记住！运行snakemake需要在该目录下运行`/path/to/your/PAPA-main`下运行**

### 0.2-install conda env

```bash
<conda/mamba> env create -f envs/snakemake_papa.yaml
```

以上是官方给出的安装代码。但是我测试过了(踩了很多坑)，最好的办法其实是：

1. 最好按照`snakemake_papa.yaml`安装对应的环境`papa_snakemake`，如果出现报错，AI选择对应的版本，记住！`Pyranges`这个版本最好不要修改，需要使用`0.0.*`。因为新版的整体语法都发生了改变，而PAPA主要还需要对位置进行处理，如果这里版本出现问题需要进行整体后续脚本的大范围修改，可以说是换血的程度了！
2. 如果后面的`papa.yaml`,`papa_r.yaml`装不上没有关系。我们可以在`papa_snakemake`分别安装其中需要的包、库。尤其是`papa_r.yaml`,主要是后续的dexseq进行差异APA的时候需要，我们可以在服务器的R环境内直接由输出的定量结果运行即可。

### 0.3-test run

```bash
cd /path/to/your/PAPA-main
conda activate papa_snakemake
snakemake -n -p --use-conda # 如果--use-conda参数的话，就会自动安装前面提到的papa.yaml和papa_r.yaml。如果你按照我上面的建议都安装在papa_snakemake环境内的话，这里不需要这个参数
```

## 1-config

### 1.1-sample file

可以替换自己的数据，参考我生成sample file的r代码

```r
sample_name <- c('KD_1','KD_3','NC_1','NC_2')
condition <- c('KD','KD','NC','NC')
path <- c(
  '/path/to/your/KD_1_sorted.bam',
  '/path/to/your/KD_3_sorted.bam',
  '/path/to/your/NC_1_sorted.bam',
  '/path/to/your/NC_2_sorted.bam'
)
fastq1 <- c(
  '/path/to/your/KD_1_1.fq.gz',
  '/path/to/your/KD_3_1.fq.gz',
  '/path/to/your/NC_1_1.fq.gz',
  '/path/to/your/NC_2_1.fq.gz'
)
fastq2 <- c(
  '/path/to/your/KD_1_2.fq.gz',
  '/path/to/your/KD_1_2.fq.gz',
  '/path/to/your/NC_1_2.fq.gz',
  '/path/to/your/NC_-2_2.fq.gz'
)
config <- data.frame(
  sample_name,condition,path,fastq1,fastq2
)
write.csv(config,'/path/to/your/sample.csv',quote = F,row.names = F)
```

**这里的bam最好sort \+ index先**

### 1.2-set config

需要修改的config：

1. working dir `main_output_dir`
2. sample file path `sample_tbl`
3. annotete GTF file path `annotation_gtf` 如果下载的gtf不是gencode需要修改`annotation_source`参数
4. genome fasta file path `genome_fasta`
5. polyA site database file path `polya_site_atlas`\
   基本需要修改的参数就这些，大多数都是一些file path，我们尽量先保持默认设置。其中一些需要的files都是后面我们需要下载的。运行前，你先下载相关的文件，在检索关键词，去对应位置设置自己的file path即可。个性化设置可以参考config内的说明以及作者在git内的说明。

## 2-download ref data

### 2.1-gtf

这里最好下载gencode的GTF，因为gencode相对来说比较全，且很多都是经过实验验证。\
我的下载链接:`https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz`\
**注意！这里有个！gencode内下载的GTF中，chromosome列是带有`chr`前缀的，如果bam文件内的chromosome是不带的话需要在这里进行统一。我的做法是将GTF修改成和bam一致的格式！**

### 2.2 ref.fa

后续PAPA需要对motif进行scan，所以需要序列信息。我的下载链接`https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz`\
这里需要使用`samtools faidx`先对fasta文件进行索引。

### 2.3-APA.bed

在数据库内下载对用的bed文件，我的下载链接:`https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz`需要进行解压缩先。

## 3-run

设置好对应的config，且在之前已经dry run过的话可以直接使用

```bash
cd PAPA-main
snakemake -p --core 8 # 如果你也是将papa.yaml内的包装到papa_snakemake环境内，不要使用--use-conda参数，反则反之
```

运行到最后的一部分，需要使用dexseq进行差异的APA分析时，可能会报错。这里我们可以直接打开RStudio server，进入`/path/to/yout/PAPA-main/scripts/run_dexseq.R`,将输入参数修改成设置参数，将需要的参数全部改成路径输入，并安装需要的包，就可以一run到底了！\
这里可能会出现的报错，但是AI一下就可以很简单找到解决办法！\
随后就可以开开心心去看结果啦\~

---

**这里总结一下易错点：**

1. 环境的安装部分很容易报错。耐心一点，如果太慢了，可以使用mamba，会比conda快一点。如果用的镜像太慢了可以换一个。
2. 一定要主要pyranges的版本，我因为这个版本问题，跑到后面一直报错！全部更新的话，工作量太大！
3. 注意gtf和bam对chromosome的写法要一致
4. 可以dry run或者使用作者的test数据run一下，确保pipeline可以跑起来先
5. snakemake的做事风格就是，做过的事情不再重复。在snakefile内会记录之前做过的事情，只要是config文件不修改就可以继续运行
6. dexseq部分我喜欢直接在RStudio内运行，我也建议这样做，更加直观一些

---

## 4-result show

IGV是个好东西啊！根据我们实验室的好习惯，得到这样类似APA或者是chip的peak，我们都喜欢在IGV内看一下，screenshot结果拿给老板看！\
具体的输出文件可以看：`https://github.com/frattalab/PAPA/blob/main/docs/output_docs.md`.其中有作者详细的说明！(作者是真的想让人用他的方法！好作者)

## 4.1-show in IGV

**Tips**\
可以先将`/path/to/differential_apa/dexseq_apa.results.processed.tsv`转成bed文件，同bam文件一并导入IGV内看

```r
library(tidyr)
library(dplyr)
puad <- fread('/path/to/differential_apa/dexseq_apa.results.processed.tsv')
pau_bed <- puad %>%
   dplyr::select(chromosome, start, end, gene_name)
pau_bed_split <- pau_bed %>%
  mutate(start_list = strsplit(start, ","),
         end_list   = strsplit(end, ",")) %>%
  rowwise() %>%
  mutate(coords = list(
    data.frame(
      start = as.numeric(start_list), 
      end   = as.numeric(end_list),
      stringsAsFactors = FALSE
    )
  )) %>%
  dplyr::select(-start, -end, -start_list, -end_list) %>%
  unnest(coords) %>%
  group_by(gene_name) %>%
  mutate(segment_id = row_number()) %>%
  ungroup() %>%
  mutate(name = paste0(gene_name, "_", segment_id))

bed_write <- pau_bed_split %>%
  mutate(start_bed = start - 1) %>%
  dplyr::select(chromosome, start_bed, end, name)

write.table(bed_write,
            file = "/data/home/hhl/03task/26_08/3_papa_new/split_pau.bed",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
```

这是我的结果中，在处理组和对照组中差异最大的位点：\
![image](3805523-20260825105711540-1092993153.png)

确实是存在一些差异的。\
但是，也是存在一些假阳性的。比如我主要关注的novel的APA，但是找出来的结果中看起来既不是novel，差异也不大。\
![image](3805523-20260825110618603-1853786337.png)

## 5-后续测试

1. 下载一直的与加尾非常相关的基因的KD数据，关注是不是会广泛引起加尾的改变
2. 找加尾改变的阳参，并下载相关数据，看PAPA是不是可以找到该基因的APA。

---

[合集目录](README.md) · [下一篇](02-dapars.md)
