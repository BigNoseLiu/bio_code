from django.db import models

# Create your models here.

class SNV(models.Model):
    chr_hg19 = models.CharField(blank=True, max_length=50, verbose_name='【hg19】染色体')
    pos1_hg19 = models.IntegerField(blank=True, null=True, verbose_name='【HG19】位置1')
    pos2_hg19 = models.IntegerField(blank=True, null=True, verbose_name='【HG19】位置2')
    ref_hg19 = models.CharField(blank=True, max_length=200, verbose_name='【HG19】ref碱基')
    alt_hg19 = models.CharField(blank=True, max_length=200, verbose_name='【HG19】alt碱基')
    chr_hg38 = models.CharField(blank=True, max_length=50, verbose_name='【hg38】染色体')
    pos1_hg38 = models.IntegerField(blank=True, null=True, verbose_name='【HG38】位置1')
    pos2_hg38 = models.IntegerField(blank=True, null=True, verbose_name='【H38】位置2')
    ref_hg38 = models.CharField(blank=True, max_length=200, verbose_name='【HG38】ref碱基')
    alt_hg38 = models.CharField(blank=True, max_length=200, verbose_name='【HG38】alt碱基')
    mut_hgvs = models.CharField(blank=True, max_length=200, verbose_name='变异HGVS名称')
    gene = models.CharField(blank=True, max_length=200, verbose_name='基因')
    nm_id = models.CharField(blank=True, max_length=200, verbose_name='标准转录本（NM）')
    c_point = models.CharField(blank=True, max_length=200, verbose_name='核苷酸变异')
    p_point = models.CharField(blank=True, max_length=200, verbose_name='氨基酸变化')
    def __str__(self):
        return self.mut_hgvs
        
class Alias(models.Model):
    type = models.CharField(blank=True, max_length=50, verbose_name='类型')
    name = models.CharField(blank=True, max_length=200, verbose_name='别名')
    class Meta:
        verbose_name = '别名库'
        
class DiseaseInher(models.Model):
    type = models.CharField(blank=True, max_length=50, verbose_name='类型')
    ch_name = models.CharField(blank=True, max_length=200, verbose_name='中文疾病名')
    en_name = models.CharField(blank=True, max_length=200, verbose_name='英文疾病名')
    ch_alias = models.ManyToManyField('Alias', blank=True, related_name="dis_alias_ch", verbose_name='中文别名')
    en_alias = models.ManyToManyField('Alias', blank=True, related_name="dis_alias_en", verbose_name='英文别名')
    omim_id = models.CharField(blank=True, max_length=200, verbose_name='OMIM编号')
    orphanet_id = models.CharField(blank=True, max_length=200, verbose_name='Orphanet编号')
    inher_mode = models.CharField(blank=True, max_length=200, verbose_name='遗传模式')
    
class Phenotype(models.Model):
    type = models.CharField(blank=True, max_length=50, verbose_name='类型')
    ch_name = models.CharField(blank=True, max_length=200, verbose_name='中文名称')
    en_name = models.CharField(blank=True, max_length=200, verbose_name='英文名称')
    ch_alias = models.ManyToManyField('Alias', blank=True, related_name="pheno_alias_ch", verbose_name='中文别名')
    en_alias = models.ManyToManyField('Alias', blank=True, related_name="pheno_alias_en", verbose_name='英文别名')
    hpo_id = models.CharField(blank=True, max_length=200, verbose_name='HPO编号')
    father_hpo_id = models.CharField(blank=True, max_length=200, verbose_name='父HPO编号')
    
class TestItem(models.Model):
    test_name = models.CharField(blank=True, max_length=200, verbose_name='检测项目名称')
    test_descr = models.TextField(blank=True, verbose_name='项目描述')

    def __str__(self):
        return self.test_name

class SampleTest(models.Model):
    sample = models.ForeignKey('Sample', null=True, on_delete=models.RESTRICT, verbose_name='样本')
    test_item = models.ForeignKey('TestItem',blank=True, null=True, on_delete=models.RESTRICT, verbose_name='检测项目')
    class Meta:
        verbose_name = "02-样本检测"
    def sample_id(self):
        return_str = "abc"
        if self.sample:
            return_str = self.sample.sample_id
        return return_str
        
    def __str__(self):
        return_str = ""
        if self.sample:
            return_str = self.sample.sample_id
        if self.test_item:
            return_str = return_str + "-" + self.test_item.test_name
        return return_str
    
class Sample(models.Model):
    sample_id = models.CharField(blank=True, unique=True, max_length=50, verbose_name='样本编号')
    sample_type = models.CharField(blank=True, max_length=50, verbose_name='样本类型')
    patient = models.ForeignKey('Patient',blank=True, null=True, on_delete=models.RESTRICT, verbose_name='患者')
    sample_descr = models.TextField(blank=True, verbose_name="样本说明")
    
    def __str__(self):
        return_str = self.sample_id
        if self.patient:
            return_str = return_str + "-" + self.patient.patient_name
        return return_str
class Patient(models.Model):
    PATYPE = [
        ('patient', '患者'),
        ('normal', '正常人'),
    ]
    patient_type = models.CharField(blank=True, max_length=20, choices=PATYPE, verbose_name='是否患病')
    patient_name = models.CharField(blank=True, max_length=50, verbose_name='姓名')
    patient_gender = models.CharField(blank=True, max_length=50, verbose_name='性别')
    patient_age = models.CharField(blank=True, max_length=50, verbose_name='年龄')
    main_pheno = models.ManyToManyField('Phenotype', blank=True, related_name="patient2pheno", verbose_name='主要表型')
    patient_clinical = models.TextField(blank=True, verbose_name='患者描述')
    files = models.ManyToManyField('File', blank=True, related_name="patient2file", verbose_name='资料')
    family = models.ForeignKey('Family',blank=True, null=True, on_delete=models.RESTRICT, verbose_name='家族')
    
class Family(models.Model):
    FATYPE = [
        ('patient', '明显'),
        ('normal', '疑似'),
        ('normal', '无'),
    ]
    family_type = models.CharField(blank=True, max_length=20, choices=FATYPE, verbose_name='家族聚集')
    family_id = models.CharField(blank=True, unique=True, max_length=50, verbose_name='家系编号')
    files = models.ManyToManyField('File', blank=True, related_name="family2file", verbose_name='资料')
    family_descr = models.TextField(blank=True, verbose_name="家系说明")

class File(models.Model):
    file = models.FileField(upload_to="files", max_length=200, verbose_name="文件上传")
    file_descr = models.TextField(blank=True, verbose_name="文件说明")

class PathSNV(models.Model):

    class Meta:
        verbose_name = '01-致病变异'
    PYTHOTYPE = [
        ('path', '致病'),
        ('likelypath', '可能致病'),
        ('benign', '良性'),
        ('likelybenign', '可能良性'),
        ('vus', 'VUS'),
    ]
    MUTTYPE = [
        ('germline', 'germline'),
        ('somatic', 'somatic'),
        ('mito', '线粒体'),
    ]
    MUTRES = [
        ('hom', '纯合'),
        ('het', '杂合'),
        ('comhet', '复合杂合'),
        ('hemi', '半合子'),
    ]
    MUTORG = [
        ('both', '父母'),
        ('mot', '母源'),
        ('fat', '父源'),
        ('novo', '新发'),
    ]
    mut_type = models.CharField(blank=True, max_length=50, choices=MUTTYPE, verbose_name='类型')
    mut = models.ForeignKey('SNV', blank=True, null=True, on_delete=models.RESTRICT, verbose_name='变异')
    sample_test = models.ForeignKey('SampleTest',blank=True, null=True, on_delete=models.RESTRICT, verbose_name='样本检测')
    mut_result = models.CharField(blank=True, max_length=50, choices=MUTRES, verbose_name='变异结果')
    father_result = models.CharField(blank=True, max_length=50, verbose_name='父亲结果')
    mother_result = models.CharField(blank=True, max_length=50, verbose_name='母亲结果')
    mut_freq = models.FloatField(blank=True, null=True, verbose_name='变异频率（%）')
    fat_mut_freq = models.FloatField(blank=True, null=True, verbose_name='父亲变异频率（%）')
    mot_mut_freq = models.FloatField(blank=True, null=True, verbose_name='母亲变异频率（%）')
    mut_origin = models.CharField(blank=True, max_length=50, choices=MUTORG, verbose_name='突变来源')
    path_type = models.CharField(blank=True, max_length=50, choices=PYTHOTYPE, verbose_name='致病性')
    #ACMG = models.ManyToManyField('ACMGEvidence', blank=True, related_name="pathsnv2acmg", verbose_name='ACMG证据')
    guidline = models.ManyToManyField('Guidline', blank=True, related_name="pathsnv2guidline", verbose_name='指南依据')
    main_pheno = models.ManyToManyField('Phenotype', blank=True, related_name="pathsnv2pheno", verbose_name='主要表型')
    main_disease = models.ManyToManyField('DiseaseInher', blank=True, related_name="pathsnv2disease", verbose_name='疾病')
 

class ACMGEvidence(models.Model):
    ACMGEV = [
        ('PVS1', 'PVS1'),
        ('PS1', 'PS1'),
        ('PS2', 'PS2'),
        ('PS3', 'PS3'),
        ('PS4', 'PS4'),
        ('PM1', 'PM1'),
        ('PM2', 'PM2'),
        ('PM3', 'PM3'),
        ('PM4', 'PM4'),
        ('PM5', 'PM5'),
        ('PM6', 'PM6'),
        ('PP1', 'PP1'),
        ('PP2', 'PP2'),
        ('PP3', 'PP3'),
        ('PP4', 'PP4'),
    ]
    ACMGLEVEL = [
        ('verystrong', 'VeryStrong'),
        ('strong', 'Strong'),
        ('moderate', 'Moderate'),
        ('suporting', 'Suporting'),
    ]
    path_snv = models.ForeignKey('PathSNV', blank=True, null=True, on_delete=models.RESTRICT, related_name="acmg2pathsnv", verbose_name='致病变异')
    acmg_evidence = models.CharField(blank=True, max_length=50, choices=ACMGEV, verbose_name='ACMG证据等级' )
    acmg_level = models.CharField(blank=True, max_length=50, choices=ACMGLEVEL, verbose_name='证据强度' )
    evidence= models.ManyToManyField('Evidence', blank=True, related_name="acmg2ev", verbose_name='参考依据')
    abstract = models.TextField(blank=True, verbose_name='描述')
    def __str__(self):
        return self.acmg_evidence + "-" + self.acmg_level
    class Meta:
        verbose_name = "ACMG证据项"
    
class Guidline(models.Model):
    guid_type = models.CharField(blank=True, max_length=50, verbose_name='类型')
    guid_name = models.CharField(blank=True, max_length=500, verbose_name='指南名称')
    guid_url = models.CharField(blank=True, max_length=500, verbose_name='链接')
    guid_file = models.ManyToManyField('File', blank=True, related_name="guidline2file", verbose_name='资料')
    guid_abstract = models.TextField(blank=True, verbose_name='指南描述')
    
class Evidence(models.Model):
    EVTYPE = [
        ('pubmed', 'PubMed'),
        ('other', '其他'),
    ]
    ev_type = models.CharField(blank=True, max_length=50, choices=EVTYPE, verbose_name='类型')
    ev_id = models.CharField(blank=True, max_length=50, verbose_name='编号')
    ev_url = models.CharField(blank=True, max_length=500, verbose_name='链接')
    ev_file = models.ManyToManyField('File', blank=True, related_name="ev2file", verbose_name='资料')
    ev_abstract = models.TextField(blank=True, verbose_name='描述')
