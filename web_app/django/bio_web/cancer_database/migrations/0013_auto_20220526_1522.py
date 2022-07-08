# Generated by Django 3.2.13 on 2022-05-26 07:22

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('cancer_database', '0012_auto_20220526_1421'),
    ]

    operations = [
        migrations.AddField(
            model_name='mutation',
            name='mut_function',
            field=models.CharField(blank=True, choices=[('LOF', 'Loss_of_function'), ('LLOF', 'Likely_Loss_of_function'), ('GOF', 'Gain_of_function'), ('LGOF', 'Likely_Gain_of_function'), ('NEU', 'Neutral'), ('LNEU', 'Likely_Neutral'), ('SWITCH', 'Switch_of_function'), ('LSWITCH', 'Likely_Switch_of_function'), ('NO', 'NO_effect'), ('unknown', 'Unknown'), ('na', '未定义')], max_length=200, verbose_name='变异功能'),
        ),
        migrations.AlterField(
            model_name='mut2evidence',
            name='support_direct',
            field=models.CharField(blank=True, choices=[('support', '支持'), ('not_support', '不支持')], max_length=100, verbose_name='证据方向'),
        ),
    ]
