# BAM&SAM文件格式, 分析工具, 具体案例(BD单细胞测序)

### SAM格式

`sam`文件是Sequence Alignment Format的简写，记录了序列比对的具体情况。`sam`文件是记录比对信息的标准结构化文件，从总体上可以分为两个部分

* SAM文件头部

以@键开头。

主要涵盖了文件标准格式版本（`VN`）、比对中使用的参考序列信息（`SQ`）、测序数据分组信息（`RG`）、比对或后期处理使用的程序信息（`PG`）上面的信息是必须存在的

还可以存在一些自定义的信息，最为常见的UMI标签信息。

* 文件比对信息

每一行是由固定的11列的组成

| 1    | QNAME | String | query序列名称                                       |
| ---- | ----- | ------ | --------------------------------------------------- |
| 2    | FLAG  | Int    | FLAG标签，主要记录比对的基本情况，取值为2的整数次幂 |
| 3    | RNAME | String | 比对至参考序列的名称，例如：chr1                    |
| 4    | POS   | Int    | 比对至参考序列的位置                                |
| 5    | MAPQ  | Int    | 比对质量，不同的比对软件对其定义不尽相同            |
| 6    | CIGAR | String | 比对的CIGAR字符串，下文有详细说明                   |
| 7    | MRNM  | String | Mate Reference NaMe                                 |
| 8    | MPOS  | Int    | Mate read比对到的参考序列最左侧的位置坐标           |
| 9    | ISIZE | Int    | 比对序列对应的模板（Template）长度                  |
| 10   | SEQ   | String | query序列                                           |
| 11   | QUAL  | String | query序列的碱基质量Phred值                          |

上述11列信息，是SAM文件必要的组成部分。除此之外，每行可以追加可选信息。该部分信息以`TAG:TYPE:VALUE` 形式存储

案例：
![sam文件](https://github.com/JaxonHello/BISN/blob/main/src/bam_sam/%E6%88%AA%E5%B1%8F2024-01-22%2013.55.17.png)

无效或者没有的字段一般用“0”或者“*”表示。

该文件来自于以下更原始的测序格式：
![origin sam文件](/Users/jaxonhe/Library/Application Support/typora-user-images/截屏2024-01-22 14.01.25.png)

*注：以上图片均来自于官方文档。*https://samtools.github.io/hts-specs/SAMv1.pdf

对于具体的字段解释：
`Qname`: r001 : read name

`FLAG`: 状态参数

![截屏2024-01-22 14.06.24](/Users/jaxonhe/Library/Application Support/typora-user-images/截屏2024-01-22 14.06.24.png)

`RNAME`: 比对上的参考序列的名字，该名字出现在Header section的@SQ行的SN标识中

`MAPQ`: 比对的质量值，计算方法为比对错误率的-10*log10的值，一般是四舍五入的整数值，如果是255，说明该比对值无效。

`CIGAR`: 8M2I4M1D3M(8 match+ 2 insection + 4 match + 1 deletion + 3 match)。除了以上的字母，其他的请参考文档

`MRNM`: 该read的mate read比对上的参考序列的名字，该名字出现在Header section的@SQ行的SN标识中，

- 如果和该read所在行的第三列“RNAME”一样，则用“=”表示，说明这对read比对到了同一条参考序列上；
- 如果mate read没有比对上，第七列则用“*”表示；
- 如果这对read没有比对到同一条参考序列，那么这一列则是mate read所在行第三列的“RNAME”。

`MPOS`: mate read比对到的参考序列“RNAME”最左侧的位置坐标

`ISIZE`: 比对长度。如果是情况1，长度是read2最右侧减去read1最左侧。情况2的话暂未达成共识，它的长度可以是read1，也可以是read2。

![截屏2024-01-22 14.36.25](/Users/jaxonhe/Library/Application Support/typora-user-images/截屏2024-01-22 14.36.25.png)

`SEQ` & `QUAL`：序列和测序质量(ASCII编码)



### BAM文件

`bam`文件作为`sam`文件的压缩版本，记录的信息本质上是一样的。但是对其进行了BGZF压缩。

该压缩的目的，在减少存储消耗的前提下，提供在有索引的前提下，提供快速的随机访问的功能。

`bam`文件的基本信息和`sam`信息存在对应关系

* Bam.bai文件

`BAM` 索引文件包含了` BAM` 文件中每个比对块（chunk）的元数据，以及其在文件中的位置。这使得对于某个给定的 genomic 区域，可以快速地找到对应的比对记录，而无需遍历整个 文件。这对于一些操作，比如可视化、子集提取和快速的查询非常有用



---

由于BAM文件是压缩文件，因此一般的linux字符处理工具（如cat、head）是无法直接处理bam文件的读取、修改、查询的。在实际处理过程中，经常使用`samtools` 和`htslib`提供的接口对`bam`文件进行处理。目前生信分析常用的主流语言中基本上都涵盖了处理`bam`文件的扩展，例如python语言中的`pysam`，R中的`Rsamtools`等。



### samtools工具(Linux)

#### samtools安装

访问 https://www.htslib.org/download/ 下载并解压`samtools-1.x`，例如`samtools-1.19.1`。

```bash
# 打开Linux
# 进入下载之后的文件夹
cd samtools-1.x
# 为防止出现权限问题，如果./configure出现问题，首先运行以下代码
chmod +x configure

./configure --prefix=/home/user/absolute/install/path
# --prefix=之后需要输入下载路径的绝对路径
# 这里使用路径就是 samtools-1.x 文件夹
make
make install
# 安装成功之后进行测试
./samtools --help
# 添加到系统路径
export PATH=/home/user/absolute/install/path/bin:$PATH
# 注意这里的路径是install路径+/bin
# 测试
samtools --help
# 或许以后在使用的时候会出现以下的问题
# Command 'samtools' not found
# 继续添加系统路径可以解决
```

#### samtools使用





---

### 案例：以BD单细胞测序结果BAM文件为例 

对于BD平台来说：对齐程序将only R2 (read2)读数与参考文件并输出与对齐质量相关的标签。

```bash
# 观察前10条比对信息
samtools view your_file.bam | head -n 10
# 结果如下：
E150025748L1C012R02200176853    0       chr1    10171   255     14M619455N135M1S        *       0       0      ACCCTAACCTAACCATTAATCCCCTGGCCCAACCCGTCATCTACTCTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCTCGCACTGATTTTTTACCTGAGTAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCC  C>AACBC8CCCCA=CCC@AC>ACCCB>CAC@;CCCA;?CC=CCCCBB@BBB=CCCCC@BA=@CC7CCCB1@C;B0CCCBB?CC@@9CCCC<CCCCCCC;B9CAA=CCB8BCCA3AC?A@?CCCCCCC?CACBCCCC??CCCCCCC?CB?C  NH:i:1  HI:i:1  AS:i:136        nM:i:0  XF:Z:__intergenic      TR:Z:*   TF:Z:*  CB:Z:3711785    MR:Z:TCGGGGAC   MA:Z:TCGGGGAC   CN:Z:T  ST:Z:03
```

| Key   | Value                        |
| ----- | ---------------------------- |
| QNAME | E150025748L1C012R02200176853 |
| FLAG  | 0                            |
| RNAME | chr1                         |
| POS   | 10171                        |
| MAPQ  | 255                          |
| CIGAR | 14M619455N135M1S             |
| MRNM  | *                            |
| MPOS  | 0                            |
| ISIZE | 0                            |
| SEQ   | ACCCTAA...                   |
| QUAL  | C>AACBC...                   |

其他信息(由BD以及序列比对工具STAR注释)

| Key                                           | Value                                  |
| --------------------------------------------- | -------------------------------------- |
| CB(cell barcode)                              | Z: 3711785 (cell barcode = 3711785)    |
| ST(sample tag)                                | Z: 03 (sampletag = 3)                  |
| NH(number of loci the reads maps to)          | Z：1(有的时候一条序列很比对到多个基因) |
| HI(multiple alignment index)                  | Z:1                                    |
| AS(alignment score: +1/-1 for match/mismatch) | Z:136                                  |
| XF(potential gene name)                       | XF:Z:IGHG2                             |
| MR(UMI)                                       | Z:TAATCCTG                             |
| CN(经过BD的算法是否认为该barcode是cell/noise) | CN:Z:T(T为细胞/X为noise)               |

*注：以上参考资料来自于BD pipeline和序列比对工具STAR的官方pdf*

[BAM and BAM Index - BD Rhapsody™ Sequence Analysis Pipeline](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/outputs/outputs_bam.html)

[STAR/doc/STARmanual.pdf at master · alexdobin/STAR (github.com)](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

```bash
# 找到细胞71的所有比对信息
samtools view ./Combined_HJY-v01.BAM | grep -E "\bCB:Z:71\b" 
# 查询条目（统计CB为798680的细胞的比对数量）
samtools view ./Combined_HJY-v01.BAM | grep -E "\bCB:Z:7986880\b" | wc -l
14351
# 该结果与seurat对象中的nCount_RNA数量不相符(推测BD存在filter过程)

# 这里的chr14:...是某个具体的基因（IGHG1 here）
samtools view ./Combined_HJY-v01.BAM chr14:105741473-105743070 | grep -E "\bCB:Z:7986880\b" | head -n 10
# XF注释中含有基因的信息
# 统计细胞7986880比对到IGHG1基因上的数量
samtools view ./Combined_HJY-v01.BAM chr14:105741473-105743070 | grep -E "\bCB:Z:7986880\b" | wc -l
41
# 而在seurat对象中，7086880细胞的IGHG1表达量为4
raw@assays$RNA@counts["IGHG1", '7986880']

# 观察到该序列的NH:i:3，比对到了3个结果
samtools view -h ./Combined_HJY-v01.BAM | grep -E "E150025748L1C032R01601632864"
# !!!!在BD中观察到了一个现象，比如说同一个细胞会同时高表达IGHG1基因和IGHG2基因，这在biology是不make sense的
# 推测到是IGHG1和IGHG2基因过于相似，因此一条序列同时比对到了IGHG1和IGHG2
# 查询结果如下
------------------------------------------------------------------
E150025748L1C032R01601632864    16      chr14   105643220       1       149M1S  *       0       0       GGGAGAGGCTCTTCTGCGTGTAGTGGTTGTGCAGAGCCTCATGCATCACGGAGCATGAGAAGACGTTCCCCTGCTGCCACCTGCTCTTGTCCACGGTGAGCTTGCTGTAGAGGAAGAAGGAGCCGTCGGAGTCCAGCATGGGAGGCGTGA      CC<D@C8BC1CB>CC+CCCCCCCCCDBDCC7BC5CD1CC8B?D?<CDD=@DB=CD;B=BCCCCC;C?@C4DB(6D+<DCCDCCCCCBC?CDCCDB/6BDBC=@1C<CDCBDC@CCCDCCACB6CCC6.?B.D(BDC=DC4B5D1:@%C<3      NH:i:3  HI:i:1  AS:i:143        nM:i:2  XF:Z:IGHG2|IGH_BReverseCONSTANTVDJ      CB:Z:7986880    MR:Z:TAATCCTG       MA:Z:TAATCCTG   CN:Z:T  ST:Z:03
-------------------------------------------------------------------
E150025748L1C032R01601632864    272     chr14   105741490       1       149M1S  *       0       0       GGGAGAGGCTCTTCTGCGTGTAGTGGTTGTGCAGAGCCTCATGCATCACGGAGCATGAGAAGACGTTCCCCTGCTGCCACCTGCTCTTGTCCACGGTGAGCTTGCTGTAGAGGAAGAAGGAGCCGTCGGAGTCCAGCATGGGAGGCGTGA      CC<D@C8BC1CB>CC+CCCCCCCCCDBDCC7BC5CD1CC8B?D?<CDD=@DB=CD;B=BCCCCC;C?@C4DB(6D+<DCCDCCCCCBC?CDCCDB/6BDBC=@1C<CDCBDC@CCCDCCACB6CCC6.?B.D(BDC=DC4B5D1:@%C<3      NH:i:3  HI:i:2  AS:i:143        nM:i:2  XF:Z:IGHG1|IGH_BReverseCONSTANTVDJ      CB:Z:7986880    MR:Z:TAATCCTG
-------------------------------------------------------------------
E150025748L1C032R01601632864    272     chr14   105741490       1       80M27764N69M1S  *       0       0       GGGAGAGGCTCTTCTGCGTGTAGTGGTTGTGCAGAGCCTCATGCATCACGGAGCATGAGAAGACGTTCCCCTGCTGCCACCTGCTCTTGTCCACGGTGAGCTTGCTGTAGAGGAAGAAGGAGCCGTCGGAGTCCAGCATGGGAGGCGTGA      CC<D@C8BC1CB>CC+CCCCCCCCCDBDCC7BC5CD1CC8B?D?<CDD=@DB=CD;B=BCCCCC;C?@C4DB(6D+<DCCDCCCCCBC?CDCCDB/6BDBC=@1C<CDCBDC@CCCDCCACB6CCC6.?B.D(BDC=DC4B5D1:@%C<3      NH:i:3  HI:i:3  AS:i:143        nM:i:1  XF:Z:__ambiguous[IGHG1|IGH_BReverseCONSTANTVDJ+IGHG3|IGH_BReverseCONSTANTVDJ]       CB:Z:7986880    MR:Z:TAATCCTG
-------------------------------------------------------------------
# 从XF结果来看确实一条序列比对到了IGHG1和IGHG2
# !!!!如何根据转录组分配没有查询到isotype的免疫细胞的isotype是一个重要的问题
```



**Last Updata: 01/24/2024**



参考教程：
http://www.xtaohub.com/BI-solutions/bam-file-format.html

https://www.jianshu.com/p/ff6187c97155


