# -*- coding: utf-8 -*-
# Generated by Django 1.11.2 on 2017-06-18 20:51
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Sequences',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('contig_id', models.CharField(max_length=128)),
                ('gene_description', models.CharField(max_length=256)),
                ('sequence', models.TextField()),
                ('species', models.CharField(max_length=128)),
            ],
            options={
                'db_table': 'sequences',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=128)),
            ],
            options={
                'db_table': 'species',
                'managed': False,
            },
        ),
    ]
