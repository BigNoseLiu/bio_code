[ PIP v1.0.0rc ]
database: PIDB202101_Simplified_v27


[ PIP v1.0.0rc1 ]
database: 
	DB version: PIDB202101_Simplified_v27 (即PIDB_v1.0.1)
fixed:
	1: 加入病毒blast验证
	2: 输出复合群的统计结果，过滤复合群中个别blast验证为假的物种，矫正复合群检出reads（矫正reads = 原复合群检出reads数 - 复合群中假阳物种的reads数），复合群统计结果输出至*xls，用于LIMS系统展示。


[ PIP v1.0.0rc2 ]
database:
	DB version: PIDB_v1.0.3
	Note: 
		1: 更新物种注释信息表 【all.Taxonomy.txt，all.Taxonomy.other.txt】
		2: 新增病毒的seqid（基因组taxid）与S1taxid(种株水平taxid)的映射表 【seqTaxid2S1LevelTaxid】
fixed:
	1: 修复病毒因分类水平、注释文件的差异导致的输出结果中“blast_maxscore”列被标记为“Flase”的bug
	2: QC统计结果中展示superkingdom水平的物种检出reads的统计结果
	3: 去掉覆盖度模块bwa mem的-k4参数，减少比对时间；bwa比对用的序列由“属及种水平reads”改为“仅种水平reads”
	4: plugin_out/Result/底下的样本结果文件，$sample.merge.xls内容替换为Total_Detail.xlsx中的内容


[ PIP v1.0.0rc3 ]
database:
	DB version: PIDB_v1.1.0
	Note:
		1: 更新物种注释信息
		2: 新增物种序列

fixed:
	1: 优化病毒的展示方式
	2: 优化by2*、plugin_out、basecaller_results等的生成顺序，适配LIMS抓取
	3: 优化LY的结果展示、序列提取等
	4: 优化覆盖度模块，更新覆盖度bwa数据库为BWA.v4，新增病毒亚种水平覆盖度图
date:
	2022/3/31


[ PIP v1.0.0rc4 ]
database:
        DB version: PIDB_v1.1.1
        Note:
                1: 更新知识库物种注释信息至病原数据库中的注释文件

fixed:
	1: 在v1.0.0rc3的基础上，更新度覆盖度模块相关的脚本；
	2: 更新覆盖度bwa数据库为BWA.v5， 补充来自PIDB的所有病毒序列、大肠杆菌基因组12个、金葡基因组11个、铜绿基因组11个;
	3: 更新覆盖度统计模块，适配BWA.v5的相关注释文件以便统计多个参考基因组上的比对信息，减少统计步骤的耗时。
date:
	2022/4/15

fixed:
	1: 设置ggplot2画图用的字体为" DejaVu Sans", 修复覆盖图字体在不同节点分析出现差异的问题;
	2: Python重写KoutScore.pl, 优化kmer得分的计算、物种检出统计步骤，修复个别物种, 如鼻病毒78，因kmer得分、blast验证统计有误而漏检的问题
date:
	2022/4/25

fixed:
	1: 增加PIDB_v1.2.0的调用
	2: 修改plugin out中单样本的merge.xls的格式，将分析结果分为四个sheet：Data, ARG, VF, Pathogen表
date:
	2022/4/26

fixed:
	1: 修改覆盖度画图模块，修复含有特殊字符的物种名的覆盖图无法绘制问题
	2: raw_summary.txt添加去宿主后reads数
	3: plugin_out中结果序列文件名加入样本编号


[ PIP v1.1.0rc ]
database:
        DB version: PIDB_v1.2.0
       
fixed:
	1: 在v1.0.0rc4的基础上，根据云康LIMS开发需求，修改分析结果输出格式，适配云康LIMS结果抓取；
date:
	2022/5/27


[ PIP v1.1.0rc2 ]
database:
        DB version: PIDB_v1.2.1
        Note:
        	1. 新增背景菌数据库，包含物种历史检出结果统计；
       
fixed:
	1: 在v1.1.0rc的基础上，根据相关人员提出的开发需求，修改分析结果输出格式，适配云康LIMS结果抓取；
	2: 使用新的覆盖度模块替换原有的覆盖度模块，新增新版覆盖度数据库的调用，更新检出序列分类脚本。
date:
	2022/5/31


[ PIP v1.1.0rc3 ]
database:
        DB version: PIDB_v1.2.1
       
fixed:
	1: 在v1.1.0rc2的基础上，根据相关人员提出的开发需求，修改分析结果输出格式，适配云康LIMS结果抓取；
	2: 个样本输出结果表中新增“Comment”列，用于比较该样本的物种在批次中的检出频率；以“|”分隔，左边为频率，右边为排位，数字越小排位越前；
	3: 补充常见菌检出矩阵
date:
	2022/6/1

fixed:
	1: 新增批次3,4级致病菌统计列表
	2: 覆盖度统计修改图片文件后缀，增加平均深度统计
date:
	2022/6/2

fixed:
	1: 更新毒力结果表输出格式
	2: Stat.py中将均一化步骤转移至Supplement.py


[ PIP v1.1.0rc4]
database:
        DB version: PIDB_v1.2.1

fixed:
	1: 新增复合群blast验证模块
	2: 耐药基因部分, 替换的人Resfinder为CARD
	3: 覆盖度美化, 加入AverageDepth
	4: 统计部分加入均一化步骤
	5: 加入流程分析进度反馈模块
date:
	2022/6/13

[ PIP v1.1.0rc5]
database:
	 DB version: PIDB_v1.2.3

fixed:
	1: 新增数据库PIDB_v1.2.3版本
	2: 脚本适配新版注释表
	3: 毒力、耐药修改文件输出格式, txt, xls, xlsx
	4: 耐药更新版本为1.0.2
	5: 优化流程分析进度反馈模块
	6: 深度图美化, 补充复合群图深度图
	7: 加入病毒分析模块
	8: 加入强阳物种（reads>=1000）批次统计图, 分析结果表过滤及排序, 即分别按属、种reads由大到小排序, 每个属筛出包含3,4级致病菌在内的前10个种水平检出结果
