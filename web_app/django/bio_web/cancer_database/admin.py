from django.contrib import admin
from import_export import resources
from import_export.admin import ImportExportActionModelAdmin
from .models import *

# Register your models here.

#class MutationAdmin(admin.ModelAdmin):
class MutationAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    list_filter  = ('type',)
    search_fields  = ('mut_name','gene')
    ordering = ['chr_hg19','pos1_hg19']
    list_display = ['id', 'gene','mut_name', 'type', 'mut_function', 'descr']
    fieldsets = (
        (None,{
              'fields':(('type','mut_name', 'mut_function'),('hotspot_mut','contain_mut'),('descr','bak_descr')),
              'classes':('extrapretty',),
        }),
        ('SNV&InDel变异信息',{
              'fields':( ('chr_hg19','pos1_hg19','pos2_hg19'), ('ref_hg19','alt_hg19'), ('gene','nm_id'), ('c_point','p_point'), ),
              'classes':('extrapretty',),
        }),
        ('Fusion变异',{
              'fields':( ('fusion_gene1','fusion_region1'), ('fusion_gene2','fusion_region2') ),
              'classes':('extrapretty',),
        }),
        ('基因扩增',{
              'fields':( ('cnv_gene',), ),
              'classes':('extrapretty',),
        }),
    )
    autocomplete_fields = ('hotspot_mut','contain_mut')
admin.site.register(Mutation,MutationAdmin)


class FileAdmin(ImportExportActionModelAdmin):
    search_fields  = ('file_descr',)
admin.site.register(File,FileAdmin)

class AliasAdmin(ImportExportActionModelAdmin):
    search_fields  = ('name',)
    list_display = [ 'name', 'type' ]
admin.site.register(Alias,AliasAdmin)


class GeneAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    list_filter  = ('type',)
    search_fields  = ('name','cancer_descr')
    list_display = [ 'id', 'name', 'type','cancer_descr','bak_descr' ]
admin.site.register(Gene,GeneAdmin)



class DiseaseAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    list_filter  = ('type',)
    search_fields  = ('ch_name','en_name')
    autocomplete_fields = ('ch_alias','en_alias','father_dis')
    list_display = [ 'id', 'ch_name', 'en_name', 'type','father_dis' ]
    #raw_id_fields = ('ch_alias','en_alias')
admin.site.register(Disease,DiseaseAdmin)

class DrugAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    list_filter  = ('type',)
    search_fields  = ('ch_name','en_name')
    list_display = [ 'id', 'ch_name', 'en_name', 'type', 'brand_name' ]
    autocomplete_fields = ('ch_alias','en_alias','related_drugs')
admin.site.register(Drug,DrugAdmin)

class Disease2drugAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    search_fields  = ('disease__ch_name','drug__ch_name','descr',)
    list_filter  = ('FDA_approve','NMPA_approve','test_required')
    list_display = [ 'id', 'disease', 'drug', 'FDA_approve', 'NMPA_approve', 'test_required','descr' ]
    autocomplete_fields = ('disease','drug')
admin.site.register(Disease2drug,Disease2drugAdmin)

class Drug2Mut_summaryAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    search_fields  = ('ch_name',)
    list_filter  = ('support_ev_level','drug_impact')
    list_display = [ 'id', 'ch_name', 'drug_impact', 'support_ev_level','descr' ]
    autocomplete_fields = ('evidences_new',)
admin.site.register(Drug2Mut_summary,Drug2Mut_summaryAdmin)

class Mut2drugAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    search_fields  = ('disease2drug__disease__ch_name','disease2drug__drug__ch_name','mut__mut_name')
    list_display = [ 'id', 'disease2drug', 'mut', 'drug_impact', 'dose_impact', 'toxi_impact'  ]
    list_filter  = ('drug_impact__support_ev_level','drug_impact__drug_impact')
    autocomplete_fields = ('disease2drug','mut','drug_impact','dose_impact','toxi_impact')
admin.site.register(Mut2drug,Mut2drugAdmin)

class Mut2evidenceAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    search_fields  = ('mut2drug__disease2drug__disease__ch_name','mut2drug__disease2drug__drug__ch_name','evidence__ev_name')
    list_display = [ 'id', 'mut2drug', 'drug_impact', 'support_direct', 'ev_level', 'evidence', ]
    list_filter  = ('drug_impact', 'support_direct', 'ev_level')
    autocomplete_fields = ('mut2drug','evidence')
admin.site.register(Mut2evidence,Mut2evidenceAdmin)

class EvidenceAdmin(ImportExportActionModelAdmin):
    list_max_show_all = 10000
    search_fields  = ('ev_name','pubmed_id','ev_abstract')
    list_display = [ 'id', 'ev_name', 'ev_type', 'pubmed_id', 'pub_version', 'nct_id', 'ev_abstract' ]
    list_filter  = ('ev_type',)
    autocomplete_fields = ('ev_file',)
admin.site.register(Evidence,EvidenceAdmin)
