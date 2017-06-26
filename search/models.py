from django.db import models


class Sequences(models.Model):
    contig_id = models.CharField(max_length=128)
    gene_description = models.CharField(max_length=256)
    sequence = models.TextField()
    species = models.CharField(max_length=128)

    def __str__(self):
        return self.contig_id + ' - ' + self.species

    class Meta:
        managed = False
        db_table = 'sequences'


class Species(models.Model):
    name = models.CharField(max_length=128)

    def __str__(self):
        return self.name

    class Meta:
        managed = False
        db_table = 'species'
        unique_together = (('id', 'name'),)
