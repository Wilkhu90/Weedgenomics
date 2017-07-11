from django.conf.urls import url
from . import views

app_name = 'search'

urlpatterns = [
    url(r'^$', views.main),
    url(r'^main/$', views.main, name='main'),
    url(r'^species/$', views.species, name='species'),
    url(r'^search_keyword/$', views.search_keyword, name='search_keyword'),
    url(r'^search_keyword_render/$', views.search_keyword_render, name='search_keyword_render'),
    url(r'^search_sequence_render/$', views.search_sequence_render, name='search_sequence_render'),
    url(r'^(?P<contig_id>[0-9]+)/$', views.detail, name='detail'),
    url(r'^download_fasta/(?P<seq_id>[0-9]+)/$', views.download_file, name='download_fasta'),
    url(r'^blastn_search/$', views.blastn_search, name='blastn_search'),
    url(r'^blastn_render/$', views.blastn_render, name='blastn_render'),
    url(r'^blastn_ncbi/(?P<seq_id>[0-9]+)/$', views.blast_at_ncbi, name='blast_at_ncbi'),
    url(r'^herbiscide_search/$', views.herbiscide_search, name='herbiscide_search'),
    url(r'^herbiscide_render/$', views.herbiscide_render, name='herbiscide_render'),
]
