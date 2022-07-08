from django.contrib import admin
from .models import *

# Register your models here.
class EvidenceAdmin(admin.ModelAdmin):
    search_fields  = ('ev_type','ev_id')
    ordering = ['ev_type','ev_id']
admin.site.register(Evidence,EvidenceAdmin)
