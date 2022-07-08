from django.db import models

# Create your models here.
class Gene(models.Model):
    GETYPE = [
        ('onco', '原癌基因'),
        ('anti', '抑癌基因'),
        ('oncoandanti', '原癌&抑癌基因'),
        ('other', '其他'),
        ('na', '未定义'),
    ]
    type = models.CharField(blank=True, null=True, max_length=200, choices=GETYPE, verbose_name='基因类型')
    name = models.CharField(blank=True, max_length=200, verbose_name='基因名称')
    cancer_descr = models.TextField( blank=True, null=True, verbose_name='肿瘤基因描述' )
    bak_descr = models.TextField( blank=True, null=True, verbose_name='备注说明' )
    class Meta:
        verbose_name = '10 基因库'
    def __str__(self):
        return self.name
        

class Mutation(models.Model):
    MUTYPE = [
        ('group', '变异类'),
        ('combine_mut', '合并变异'),
        ('SNV', 'SNV'),
        ('InDel', 'InDel'),
        ('Fusion', 'Fusion'),
        ('CNV', 'CNV'),
        ('other', '其他'),
        ('na', '未定义'),
    ]

    PREFF = [
        ('LOF', 'Loss_of_function'),
        ('LLOF', 'Likely_Loss_of_function'),
        ('GOF', 'Gain_of_function'),
        ('LGOF', 'Likely_Gain_of_function'),
        ('NEU', 'Neutral'),
        ('LNEU', 'Likely_Neutral'),
        ('SWITCH', 'Switch_of_function'),
        ('LSWITCH', 'Likely_Switch_of_function'),
        ('NO', 'NO_effect'),
        ('unknown', 'Unknown'),
        ('na', '未定义'),
    ]

    type = models.CharField(blank=True, null=True, max_length=200, choices=MUTYPE, verbose_name='变异类型')
    mut_function = models.CharField(blank=True, null=True, max_length=200, choices=PREFF, verbose_name='变异功能')
    mut_name = models.CharField(blank=True, max_length=200, verbose_name='变异名称')
    descr = models.TextField( blank=True, null=True, verbose_name='变异位点描述' )
    bak_descr = models.TextField( blank=True, null=True, verbose_name='备注说明' )
    hotspot_mut = models.ManyToManyField('Mutation', blank=True, related_name="hotspot_muts", verbose_name='热点突变')
    contain_mut = models.ManyToManyField('Mutation', blank=True, related_name="contain_muts", verbose_name='包含突变')

    #SNV&InDel 变异信息
    chr_hg19 = models.CharField(blank=True, null=True, max_length=50, verbose_name='【HG19】染色体')
    pos1_hg19 = models.IntegerField(blank=True, null=True, verbose_name='【HG19】位置1')
    pos2_hg19 = models.IntegerField(blank=True, null=True, verbose_name='【HG19】位置2')
    ref_hg19 = models.CharField(blank=True, null=True, max_length=200, verbose_name='【HG19】ref碱基')
    alt_hg19 = models.CharField(blank=True, null=True, max_length=200, verbose_name='【HG19】alt碱基')
    chr_hg38 = models.CharField(blank=True, null=True, max_length=50, verbose_name='【HG38】染色体')
    pos1_hg38 = models.IntegerField(blank=True, null=True, verbose_name='【HG38】位置1')
    pos2_hg38 = models.IntegerField(blank=True, null=True, verbose_name='【H38】位置2')
    ref_hg38 = models.CharField(blank=True, null=True, max_length=200, verbose_name='【HG38】ref碱基')
    alt_hg38 = models.CharField(blank=True, null=True, max_length=200, verbose_name='【HG38】alt碱基')
    gene = models.CharField(blank=True, null=True, max_length=200, verbose_name='基因')
    nm_id = models.CharField(blank=True, null=True, max_length=200, verbose_name='标准转录本（NM_id）')
    c_point = models.CharField(blank=True, null=True, max_length=200, verbose_name='核苷酸变异(c.)')
    p_point = models.CharField(blank=True, null=True, max_length=200, verbose_name='氨基酸变化(p.)')
    #基因扩增
    cnv_gene = models.CharField(blank=True, null=True, max_length=200, verbose_name='基因扩增')

    #Fusion变异
    fusion_gene1 = models.CharField(blank=True, null=True, max_length=200, verbose_name='融合基因1')
    fusion_region1 = models.CharField(blank=True, null=True, max_length=200, verbose_name='区域1')
    fusion_gene2 = models.CharField(blank=True, null=True, max_length=200, verbose_name='融合基因2')
    fusion_region2 = models.CharField(blank=True, null=True, max_length=200, verbose_name='区域2')
    
    #蛋白影响
    protein_effect = models.CharField(blank=True, null=True, max_length=200, choices=MUTYPE, verbose_name='变异类型')

    class Meta:
        verbose_name = '13 变异库'
    def __str__(self):
        return self.mut_name
        
class Alias(models.Model):
    type = models.CharField(blank=True, max_length=50, verbose_name='类型')
    name = models.CharField(blank=True, max_length=200, verbose_name='别名')
    class Meta:
        verbose_name = '91 别名库'
    def __str__(self):
        return self.name
        

class Disease(models.Model):
    DITYPE = [
        ('cancer', '癌症'),
        ('sub_cancer', '癌症亚型'),
        ('other', '其他'),
        ('na', '未定义'),
    ]
    type = models.CharField(blank=True, null=True, max_length=50, choices=DITYPE, verbose_name='类型')

    father_dis = models.ForeignKey('Disease', blank=True, null=True, on_delete=models.RESTRICT, related_name="father_disease", verbose_name='所属父疾病')
    ch_name = models.CharField(blank=True, max_length=200, verbose_name='中文疾病名')
    en_name = models.CharField(blank=True, null=True, max_length=200, verbose_name='英文疾病名')
    ch_alias = models.ManyToManyField('Alias', blank=True, related_name="dis_alias_ch", verbose_name='中文别名')
    en_alias = models.ManyToManyField('Alias', blank=True, related_name="dis_alias_en", verbose_name='英文别名')
    do_id = models.CharField(blank=True,null=True, max_length=200, verbose_name='DOID编号')
    omim_id = models.CharField(blank=True,null=True, max_length=200, verbose_name='OMIM编号')
    orphanet_id = models.CharField(blank=True,null=True, max_length=200, verbose_name='Orphanet编号')
    inher_mode = models.CharField(blank=True,null=True, max_length=200, verbose_name='遗传模式')
    descr = models.TextField( blank=True,null=True, verbose_name='疾病描述' )

    class Meta:
        verbose_name = '11 疾病库'
    def __str__(self):
        return self.ch_name


class Drug(models.Model):
    DRTYPE = [
        ('drug_group', '药物类'),
        ('combine', '联合用药'),
        ('molecular', '小分子药物'),
        ('mAb', '单克隆抗体'),
        ('fusion_protein', '融合蛋白'),
        ('alkylating', '烷基化药物'),
        ('chemo', '化疗药'),
        ('靶向治疗药物', '靶向治疗药物'),
        ('免疫治疗药物', '免疫治疗药物'),
        ('内分泌治疗药物', '内分泌治疗药物'),
        ('other', '其他'),
        ('na', '未定义'),
    ]
    type = models.CharField(blank=True,null=True, max_length=50, choices=DRTYPE, verbose_name='类型')
    brand_name = models.CharField(blank=True,null=True, max_length=200, verbose_name='药物商品名')
    ch_name = models.CharField(blank=True, max_length=200, verbose_name='中文药物名')
    en_name = models.CharField(blank=True, null=True, max_length=200, verbose_name='英文药物名')
    ch_alias = models.ManyToManyField('Alias', blank=True, related_name="drug_alias_ch", verbose_name='中文别名')
    en_alias = models.ManyToManyField('Alias', blank=True, related_name="drug_alias_en", verbose_name='英文别名')
    related_drugs = models.ManyToManyField('Drug', blank=True, related_name="drugs_in_drug", verbose_name='包含药物')
    drugbank_id = models.CharField(blank=True, null=True, max_length=200, verbose_name='drugbank编号')
    descr = models.TextField( blank=True, null=True, verbose_name='药物描述' )

    class Meta:
        verbose_name = '12 药物库'
    def __str__(self):
        return self.ch_name


REQSTATUS = [
        ('yes', '用药前要求基因检测'),
        ('no', '不必须基因检测'),
        ('na', '未定义'),
]
class Disease2drug(models.Model):
    PVSTATUS = [
        ('yes', '批准'),
        ('no', '未批准'),
        ('na', '未定义'),
    ]
    disease = models.ForeignKey('Disease', blank=True, null=True, on_delete=models.RESTRICT, related_name="drug_relate_disease", verbose_name='相关疾病')
    drug = models.ForeignKey('Drug', blank=True, null=True, on_delete=models.RESTRICT, related_name="disease_relate_drug", verbose_name='相关药物')
    FDA_approve = models.CharField(blank=True,null=True,  max_length=200, choices=PVSTATUS, verbose_name='FDA批准状态')
    NMPA_approve = models.CharField(blank=True,null=True,  max_length=200, choices=PVSTATUS, verbose_name='NMPA批准状态')
    test_required = models.CharField(blank=True,null=True,  max_length=200, choices=REQSTATUS, verbose_name='是否要求基因检测')
    descr = models.TextField( blank=True,null=True,  verbose_name='药物说明' )
    fda_descr = models.TextField( blank=True,null=True,  verbose_name='FDA批文说明' )
    nmpa_descr = models.TextField( blank=True,null=True,  verbose_name='NMPA批文说明' )

    class Meta:
        verbose_name = '21 疾病药物关系'
    def __str__(self):
        return self.disease.ch_name + ' vs ' + self.drug.ch_name


#    IMP = [
#        ('sensitive', '敏感'),
#        ('resistance', '耐药'),
#        ('na', '未定义'),
#    ]


EVLEVEL = [
    ('A', 'A级'),
    ('B', 'B级'),
    ('C', 'C级'),
    ('D', 'D级'),
    ('default', '默认（同证据类型）'),
    ('na', '未定义'),
]
class Drug2Mut_summary(models.Model):
    DRUGIMP = [
        ('sensitive', '敏感'),
        ('resistance', '耐药'),
        ('decreased_response', '疗效降低'),
        ('no_response', '用药无效'),
        ('increase_dose', '增加用量'),
        ('decrease_dose', '降低用量'),
        ('increase_toxicity', '毒副作用增强'),
        ('decrease_toxicity', '毒副作用降低'),
        ('prog_better', '预后好'),
        ('prog_worse', '预后差'),
        ('diagnosis', '诊断'),
        ('no_effect', '无影响'),
        ('na', '未定义'),
    ]
    ch_name = models.CharField(blank=True, max_length=500, verbose_name='名称')
    support_ev_level = models.CharField(blank=True,null=True, max_length=50, choices=EVLEVEL, verbose_name='最高证据等级')
    drug_impact = models.CharField(blank=True,null=True, max_length=50, choices=DRUGIMP, verbose_name='临床影响')
    descr = models.TextField(blank=True,null=True, verbose_name="临床意义描述")
    evidences = models.ManyToManyField('Mut2evidence', blank=True, related_name="ev2summary", verbose_name='参考原始证据')
    evidences_new = models.ManyToManyField('Evidence', blank=True, related_name="ev_new2summary", verbose_name='参考证据')
    class Meta:
        verbose_name = '31 临床意义小结'
    def __str__(self):
        t_descr = '-'
        t_ev = '-'
        t_impact = '-'
        if self.descr is not None:
            t_descr = self.descr
        if self.support_ev_level is not None:
            t_ev = self.support_ev_level
        if self.drug_impact is not None:
            t_impact = self.drug_impact
        
        return '【' + t_impact + '|' + t_ev + '】' + t_descr

MUT_LEVEL = [
        ('1', 'I级'),
        ('2', 'II级'),
        ('3', 'III级'),
        ('4', 'IV级'),
        ('na', '未定义'),
]
class Mut2drug(models.Model):
    disease2drug = models.ForeignKey('Disease2drug', blank=True, null=True, on_delete=models.RESTRICT, related_name="mut2drug_disease2drug", verbose_name='相关疾病药物')
    mut = models.ForeignKey('Mutation', blank=True, null=True, on_delete=models.RESTRICT, related_name="mut2drug_mut", verbose_name='相关变异')
    mut_level = models.CharField(blank=True, null=True,  max_length=50, choices=MUT_LEVEL, verbose_name='治疗分级')
    test_required = models.CharField(blank=True, null=True,  max_length=200, choices=REQSTATUS, verbose_name='用药前是否要求该基因检测')
    drug_impact = models.ForeignKey('Drug2Mut_summary', blank=True, null=True, on_delete=models.RESTRICT, related_name="drug_impact2summary", verbose_name='药效影响')
    dose_impact = models.ForeignKey('Drug2Mut_summary', blank=True, null=True, on_delete=models.RESTRICT, related_name="drug_dose", verbose_name='药量影响')
    toxi_impact = models.ForeignKey('Drug2Mut_summary', blank=True, null=True, on_delete=models.RESTRICT, related_name="drug_toxi2summary", verbose_name='毒副作用')

    class Meta:
        verbose_name = '22 药物变异关系'
    def __str__(self):
        t_mut_name = '-'
        if self.mut is not None:
            t_mut_name = self.mut.mut_name
        t_dis_name = '-'
        if self.disease2drug is not None:
            t_dis_name = self.disease2drug.__str__()
        return t_mut_name + '( ' + t_dis_name + ' )'

class Mut2evidence(models.Model):
    IMP = [
        ('sensitive', '敏感'),
        ('resistance', '耐药'),
        ('decreased_response', '疗效降低'),
        ('no_response', '用药无效'),
        ('increase_dose', '增加用量'),
        ('decrease_dose', '降低用量'),
        ('increase_toxicity', '毒副作用增强'),
        ('decrease_toxicity', '毒副作用降低'),
        ('na', '未定义'),
    ]
    SUP = [
        ('support', '支持'),
        ('not_support', '不支持'),
        ('na', '未定义'),
    ]
    #disease2drug = models.ForeignKey('Disease2drug', blank=True, null=True, on_delete=models.RESTRICT, related_name="mu2evidence_disease2drug", verbose_name='相关疾病药物')
    mut2drug = models.ForeignKey('Mut2drug', blank=True, null=True, on_delete=models.RESTRICT, related_name="Mut2drug_evidence", verbose_name='药物变异')
    drug_impact = models.CharField(blank=True, max_length=100, choices=IMP, verbose_name='用药影响')
    support_direct = models.CharField(blank=True, max_length=100, choices=SUP, verbose_name='证据方向')
    ev_level= models.CharField(blank=True, max_length=50, choices=EVLEVEL, default='na', verbose_name='证据等级')
    evidence = models.ForeignKey('Evidence', blank=True, null=True, on_delete=models.RESTRICT, related_name="mu2evidence_evidence", verbose_name='相关证据')
    descr = models.TextField(blank=True, verbose_name="解释说明")

    class Meta:
        verbose_name = '41 药物变异相关原始证据'
 #   def __str__(self):
  #      return self.mut2drug.mut.mut_name



EVTYPE = [
    ('101_NCCN', 'A级_NCCN临床指南'),
    ('102_CSCO', 'A级_CSCO临床指南'),
    ('103_NMPA', 'A级_NMPA批准'),
    ('104_FDA', 'A级_FDA批准'),
    ('105_otherGuidline', 'A级_其他指南'),
    ('201_consensus', 'B级_专家共识'),
    ('301_CTPhase3', 'B级_III期临床'),
    ('302_CTBig', 'B级_其他大规模临床'),
    ('401_CTPhase2', 'B/C级_II期临床'),
    ('501_CTPhase1', 'B/C级_I期临床'),
    ('502_CTSmall', 'C级_其他小规模临床'),
    ('601_Preclinical', 'D级_临床前研究'),
    ('602_Case', 'D级_病例报道'),
    ('603_Other', 'D级_其他研究'),
    ('na', '未定义'),
]
class Evidence(models.Model):
    ev_type = models.CharField(blank=True, null=True, max_length=50, choices=EVTYPE,  default='na', verbose_name='证据类型')
    ev_name= models.CharField(blank=True, null=True, max_length=500, verbose_name='名称')
    pub_version = models.CharField(blank=True, null=True, max_length=200, verbose_name='版本/发表日期')
    pubmed_id = models.CharField(blank=True, null=True, max_length=50, verbose_name='pubMed编号')
    nct_id = models.CharField(blank=True, null=True, max_length=50, verbose_name='NCT编号')
    ev_url = models.CharField(blank=True, null=True, max_length=500, verbose_name='链接')
    ev_file = models.ManyToManyField('File', blank=True, related_name="ev2file", verbose_name='资料文件')
    ev_abstract = models.TextField( blank=True, null=True, verbose_name='摘要' )

    class Meta:
        verbose_name = '61 证据资料'
    def __str__(self):
        if self.ev_name is not None:
            if self.pub_version is not None:
                return self.ev_name + '(' + self.pub_version + ')'
            return self.ev_name
        if self.ev_abstract is not None:
            return self.ev_abstract
        return self.id

class File(models.Model):
    file = models.FileField(upload_to="files", max_length=200, verbose_name="文件上传")
    file_descr = models.TextField(blank=True, verbose_name="文件说明")
    class Meta:
        verbose_name = '92 文件管理'

