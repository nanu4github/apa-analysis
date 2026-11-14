# Darpas1/2分析APA

作者 · 不加酵母的面团（DoughWithoutYeast）

原发表于 2026-08-25 · [博客园原文](https://www.cnblogs.com/DoughWithoutYeast/p/22681595)

原文更新于 2026-10-26

## Darpars 1/2 介绍

dapars家族主要原理是依据：\
**利用标准 RNA-seq 在 3′UTR 上的覆盖度变化，来推断可变多聚腺苷酸化（APA）事件。**\
![image](3805523-20260825162532241-1552659639.png)\
目前，有两个版本。两个版本各有各的特点。

* 对于Dapars1：可以进行处理组和对照组的差异APA的分析
* 对于Dapars2：可以进行多线程大批量数据并行处理\
  大家可以根据自己的需求进行选择。

## Dapars1

### 0-installation

`https://github.com/ZhengXia/dapars`\
下载git内的源码压缩包👆\
解压之后主要函数在src内

```
ls -l src
-rw-rw-r-- 1 hhl hhl  5499 10月 16  2022 DaPars_Extract_Anno.py
-rw-rw-r-- 1 hhl hhl 23962 10月 16  2022 DaPars_main.py
```

### 1-preparation

#### 1.1 下载UCSC内需要的注释文件并注释3'-UTR

* 下载全基因组的bed文件\
  ![Clipboard_2025-11-18-10-49-28](3805523-20260825160717130-1491228162.png)
* 下载id转换文件\
  ![Clipboard_2025-11-18-10-50-15](3805523-20260825160728166-2015198762.png)
* 选择两个name\
  ![Clipboard_2025-11-18-10-50-31](3805523-20260825160738012-505109725.png)

```bash
cd /path/to/dapars2_db
python /path/to/DaPars2-2.1/src/DaPars_Extract_Anno.py -b hg38_wholeGene_annotation.bed -s hg38_Refseq_id_from_UCSC.txt -o hg38_refseq_extracted_3UTR.bed
```

**所有的dapars2需要的文件都保存在了上面的文件夹内**

#### 1.2 生成wig

```bash
for b in /path/to/bam/*bam; do
    sample=$(basename $b .bam)
    bedtools genomecov -ibam $b -bga -split \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' \
    | sort -k1,1 -k2,2n \
    > ${sample}.wig
done
ls *.wig | while read a; do
    sed -i '1d' $a        # 删除 trackline
    sed -i 's/^/chr/' $a  # 对每一行行首加 chr
done

```

**这里注意有坑：要看wig文件内的chr的格式，是不是和3'UTR的一致，不一致会报错的**\
查看染色体及覆盖度

```bash
cut -f1 test.sorted.bam.out.wig  | uniq  -c
```

挑选出需要分析的染色体

```bash
cat /path/to/dapars2_db/hg38_refseq_extracted_3UTR.bed |cut -f1|sort|uniq > chrList.txt
```

#### 1.3 拿到测序深度文件

```bash
ls *.wig | sed 's/.wig$//' | while read prefix; do
    bam="/path/to/bam/${prefix%.bam.out}.bam"

    if [[ ! -f "$bam" ]]; then
        echo "WARNING: 找不到 BAM: $bam" >&2
        continue
    fi

    # samtools view -c 得到 mapped reads（默认不会算 unmapped）
    depth=$(samtools view -c -F 4 "$bam")

    echo -e "${prefix}.wig\t${depth}"
done > sequencing_depth_file.txt
```

### 2-config

对将上面的准备文件替换到对应的config文件内的位置上就可以啦\~

```
#The following file is the result of step 1.

Annotated_3UTR=/path/to/dapars2_db/hg38_refseq_extracted_3UTR.bed

#A comma-separated list of BedGraph files of samples from condition 1

Group1_Tophat_aligned_Wig=/path/to/kd_1_sorted.wig,/path/to/kd_2_sorted.wig
#Group1_Tophat_aligned_Wig=Condition_A_chrX_r1.wig,Condition_A_chrX_r2.wig if multiple files in one group

#A comma-separated list of BedGraph files of samples from condition 2

Group2_Tophat_aligned_Wig=/path/to/nc_1_sorted.wig,/path/to/nc_2_sorted.wig

Output_directory=DaPars1/

Output_result_file=DaPars1

#At least how many samples passing the coverage threshold in two conditions
Num_least_in_group1=2

Num_least_in_group2=2

Coverage_cutoff=30

#Cutoff for FDR of P-values from Fisher exact test.

FDR_cutoff=0.05


PDUI_cutoff=0.5

Fold_change_cutoff=0.59
```

### 3-run

```bash
python /path/to/dapars-1.0.0/src/DaPars_main.py dapars1_config
```

### 4-result

![image](3805523-20260825165701774-1840075907.png)\
大致的结果是这样的。\
这里我们主要需要关注以下几列：

* "Group\_A\_Mean\_PDUI"：config内A组的平均PDUI,这里我们是KD组的平均PDUI
* "Group\_B\_Mean\_PDUI"：config内B组的平均PDUI
* "PDUI\_Group\_diff"：两组之间差异PDUI
* "P\_val"：差异PDUI检验的pvalue
* "adjusted.P\_val"：pvalue的校正
* "Pass\_Filter"：是否通过config内设置的cutoff

## Dapars2

### 0-installation

* 直接从github上下载dapars2的源码之后加压缩使用。`https://github.com/3UTR/DaPars2`
  + 基本的函数在scr目录下`/DataM/hhl/software/DaPars2-2.1/src`
  + python scr/function.py 参数s
* 参考教程`https://blog.csdn.net/qq_61543229/article/details/128186673` `https://www.jianshu.com/p/a75ceb3adbcb`
* 需要bedtools的环境
* 官方教程：`https://github.com/3UTR/DaPars2/wiki`

### 1-preparation

**参照dapars1**\
拿到需要的文件有

* 所有sample的wig文件
* chrList.txt
* annotation 3'UTR bed

### 2-config

将上面处理得到的文件在config内进行替换就可以啦\~

```
# Specify the reference of 3'UTR region
 
Annotated_3UTR=/path/to/hg38_refseq_extracted_3UTR.bed
 
# A comma separated list of wig files of all samples
 
Aligned_Wig_files=/path/to/sample_1_sorted.wig,/path/to/sample_2_sorted.wig,/path/to/sample_3_sorted.wig,/path/to/sample_4_sorted.wig
 
Output_directory=Dapars2/
 
Output_result_file=Dapars2
 
# Specify Coverage threshold
 
Coverage_threshold=10
 
# Specify the number of threads to process the analysis
 
Num_Threads=8
 
# Provide sequencing depth file for normalization
 
sequencing_depth_file=sequencing_depth_file.txt
```

### 3-run

```bash
python /path/to/DaPars2-2.1/src/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure_file chrList.txt
```

**整理所有的chr成为一个all文件**

```bash
#!/bin/bash
out=DaPars2_merged.txt
> $out

for chr in chr{1..22} chrX chrY chrM; do
    f=Dapars2_$chr/Dapars2_result_temp.$chr.txt
    [[ -f $f ]] || continue
    [[ ! -s $out ]] && cat $f >> $out || tail -n +2 $f >> $out
done
```

### 4-result

结果大致的样子：\
![image](3805523-20260825163414952-974445227.png)\
相对比起Dapars1的版本少了差异APA检验结果的部分

# 结果可视化

* igv展示那些APA位点，show这些结果给老板看
* 将所有的Dapars1分析出来的位点绘制火山图

---

[合集目录](README.md) · [上一篇](01-papa.md)
