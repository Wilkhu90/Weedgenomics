from Bio import pairwise2
from Bio.Seq import Seq
from django.http import Http404
from django.shortcuts import render
from .models import Sequences, herbiscide
from django.http import HttpResponse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import tempfile


def index(request):
    template = "search/index.html"
    context = {}
    return render(request, template, context)


def main(request):
    template = "search/main.html"
    context = {}
    return render(request, template, context)


def species(request):
    template = "search/species.html"
    context = {}
    return render(request, template, context)


# Sequences.objects.exclude(gene_description__icontains=query).filter(species__iexact=species)
def search_keyword(request):
    query = request.GET.get("query")
    species = request.GET.get("species")
    all_sequences = []
    if len(query) > 0:
        all_sequences = Sequences.objects.filter(gene_description__icontains=query).filter(species__iexact=species)
    context = {"results": all_sequences, "searchQuery": { "query": query, "species": species}}
    template = "search/results.html"
    return render(request, template, context)


def detail(request, contig_id):
    try:
        sequence = Sequences.objects.get(pk=contig_id)
    except Sequences.DoesNotExist:
        raise Http404("Sequence Does Not Exist")
    context = {"sequence": sequence}
    template = "search/details.html"
    return render(request, template, context)


def search_sequence(request):
    seq_query = request.GET.get("seq_query")
    threshold = request.GET.get("threshold")
    species = request.GET.get("species")
    keyword = request.GET.get("keyword")
    all_sequences = Sequences.objects.all()
    count = 0
    total = len(all_sequences)

    if threshold == "None":
        threshold = len(seq_query)*2
    hit_list = []
    print("searching")
    for sequence in all_sequences:
        count += 1
        inc = int(total/100)
        if count % inc == 0:
            print(str(int((count/total)*100))+' % completed')
        if sequence.species == species and keyword in sequence.gene_description:
            record = Seq(sequence.sequence)
            alignments = pairwise2.align.localms(record, seq_query, 2, -2, -1, -0.5)
            for alignment in alignments:
                if alignment[2] == threshold:
                    # change to sequence based on the results returned
                    hit_list.append(sequence)
                    break

    context = {"hit_list": hit_list}
    template = "search/sequence_results.html"
    return render(request, template, context)


def search_sequence_render(request):
    return render(request, "search/sequence_search.html", {})


def search_keyword_render(request):
    return render(request, "search/keyword_search.html", {})


def download_file(request, seq_id):
    seq_id = int(seq_id)
    sequence = Sequences.objects.get(id=seq_id)
    temp_sequence = str(sequence.sequence)

    # insert next line for every 80 characters
    length = len(temp_sequence)

    temp_list = []
    begin = 0
    for i in range(0, length):
        if i % 80 == 0:
            temp_list.append(temp_sequence[begin:i]+'\n')
            begin = i
    temp_list.append(temp_sequence[begin:])

    temp_sequence = ''.join(temp_list)

    my_sequence = "> "+sequence.contig_id+" | "+sequence.gene_description
    my_sequence = my_sequence+temp_sequence

    temp_file = tempfile.TemporaryFile()
    my_sequence = my_sequence.encode()

    temp_file.write(my_sequence)
    temp_file.seek(0)

    response = HttpResponse(temp_file, content_type='application/fasta;charset=UTF-8')
    response['Content-Disposition'] = "attachment; filename=%s" % sequence.species+".fasta"

    temp_file.close()
    return response


def blastn_search(request):
    name = request.GET.get("query")
    database_name = request.GET.get("species")
    query = request.GET.get("query")

    # check which database because each database has different names
    if database_name == 'Cyperus_esculentus':
        database = 'search/data/Cyperus_esculentus/YNS_NewNames.fasta'
    if database_name == 'Eleusine_indica':
        database = 'search/data/Eleusine_indica/EleusineIndicaFinal.fasta'
    if database_name == 'Cyperus_rotundus':
        database = 'search/data/Cyperus_rotundus/PNS_TrinityNewNames2.fasta'
    if database_name == 'Poa_annua_infirma':
        database = 'search/data/Poa_annua_infirma/AnnuaInfirmaHomeoFinal.fasta'
    if database_name == 'Poa_annua_supina':
        database = 'search/data/Poa_annua_supina/AnnuaSupinaHomeoFinal.fasta'
    if database_name == 'Poa_infirma':
        database = 'search/data/Poa_infirma/InfirmaFinal.fasta'
    if database_name == 'Poa_supina':
        database = 'search/data/Poa_supina/SupinaFinal.fasta'

    hits = []
    if len(query) > 0:
        # write query in fasta format
        query = '>test query \n'+query
        file = open('search/tempfiles/query.fasta', 'w')
        file.write(query)
        file.close()

        # select files
        query = "search/tempfiles/query.fasta"
        output = "search/tempfiles/search.xml"

        cline_blast = NcbiblastnCommandline(query=query, db=database, evalue=0.00001, outfmt=5, out=output)
        stdout, stderr = cline_blast()

        # read output file
        results_handle = open(output, 'r')
        blast_records = NCBIXML.read(results_handle)
        results_handle.close()

        results = blast_records.alignments

        for result in results:
            for hit in result.hsps:
                # get the id based on contg_id / hit.hit_id
                sequence = Sequences.objects.get(contig_id=result.hit_id)
                hits.append({"Hit_exp": hit.expect, "Hit_query": hit.query,
                             "Hit_match": hit.match, "Hit_sbject": hit.sbjct,
                             "description": sequence.gene_description, "ID": sequence.id, "Hit_id": result.hit_id})

    context = {"hits": hits, "searchQuery": { "species": database_name, "query": name }}
    return render(request, "search/blastn_results.html", context)


def blastn_render(request):
    return render(request, "search/blastn_search.html", {})

def herbiscide_search(request):
    geneId= request.GET.get("geneId")
    genus = request.GET.get("genus")

    herbs = herbiscide.objects.filter(gene_id__iexact=geneId).filter(genus__iexact=genus)
    context = {"results": herbs, "searchQuery": { "geneId": geneId, "genus": genus}}

    return render(request, "search/herbiscide_results.html", context)

def herbiscide_render(request):
    return render(request, "search/herbiscide_search.html", {})
