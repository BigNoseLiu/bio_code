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
