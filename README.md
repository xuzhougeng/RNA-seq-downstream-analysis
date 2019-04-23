第一步: 数据上传并进行解析

第二步:



挑选出背景和对照组，设定阈值





## 变量命名

InputId

上传

- upload_vcf : 上传VCF
- EMS: 用于选择是否为EMS诱变群体，（目前仅仅是做EMS过滤，未来可能会单独写个脚本做滑窗的index计算）
- submit: 确定上传文件，进行处理

候选基因筛选

- input_bg: 选择背景样本
- input_fg: 选择突变呀根本
- input_delta: 筛选的差值index计算
- submit2: 开始进行候选位点筛选

画图部分

- input_chr



输出部分

- ratioplot
- download_ratioplot
- deltaplot
- download_deltaplot
- dataset1
- download_dataset1
- dataset2
- download_dataset2



全部变量(global_value)

- output_vcf
- chr_selection
- samples
- candidate_sites