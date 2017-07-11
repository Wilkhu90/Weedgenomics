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

class herbiscide(models.Model):
    gene_id = models.CharField(max_length=128)
    genus = models.CharField(max_length=128)
    organism = models.TextField()
    genbankId = models.TextField()
    description = models.TextField()
    submissionDate = models.TextField()
    seqTechnology = models.TextField()
    taxonomy = models.TextField()
    sequence = models.TextField()
    sequenceLength = models.TextField()
    authors = models.TextField()
    comment = models.TextField()
    consrtm = models.TextField()
    journal = models.TextField()
    location = models.TextField()
    medline_id = models.TextField()
    pubmed_id = models.TextField()
    title = models.TextField()

    def __str__(self):
        return self.gene_id + ' - ' + self.genus

    class Meta:
        managed = False
        db_table = 'herbiscide'
