# RNA-seq downstream analysis 

This is a shiny app for RNA-seq downstream analysis for _Arabidopsis thaliana_.

It can do the following analysis:

- Differetial gene analysis using DESeq2
- GO and KEGG enrichment analysis using clusterProfiler
- GSEA analysis using clusterProfiler

The only data you need provided is expression matrix. 

I use `read.table("data.txt",sep="\t",header = TRUE)` to parse the data.

|Geneid | treat1 | tread2 | untreat1 | untreat2 |
| --- | --- | ----| --- | --- |
|AT1G01010|	163|	168|	233	|219|
|AT1G01020|	289|	304|	390|	451|
|AT1G01030|	32|	27|	37	|25|
|AT1G01040|	1430|	1350|	1706	|1593|
|AT1G01046|	14|	13|	18|	19|
|AT1G01050|	1646|	1637|	1904|	1823|
|AT1G01060|	132	|143|	152|	103|
|AT1G01070|	28	|30	|31|	16|
|AT1G01073|	0|	0|0	|0|

Requirement:

- R > 3.5.0
- R packages:
  - shiny 
  - DESeq2
  - stringr
  - tidyr
  - ggplot2
  - plotly
  - shinythemes
  - clusterProfiler
  - org.At.tair.db

## 功能升级

- 修改界面, 以shinydashboard 作为主题
- 增加文件下载功能
- 增加更多的注释信息