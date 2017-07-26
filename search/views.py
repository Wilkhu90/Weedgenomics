from django.http import Http404
from django.shortcuts import render
from .models import Sequences, herbiscide
from django.http import HttpResponse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from django.core.mail import send_mail
from .forms import ContactForm
from django.conf import settings
import tempfile
import os

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


def team(request):
    template = "search/team.html"
    context = {}
    return render(request, template, context)

def sponsor(request):
    template = "search/sponsor.html"
    context = {}
    return render(request, template, context)

# Keyword search
def search_keyword(request):
    query = request.GET.get("query")
    species = request.GET.get("species")
    sql_query = 'SELECT * FROM sequences Where species = "'+species+'"'+' and gene_description LIKE "%%'+query+'%%"'
    all_sequences = []
    if len(query) > 0:
        all_sequences = Sequences.objects.raw(sql_query)
    context = {"results": all_sequences, "searchQuery": {"query": query, "species": species}}
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

def herbicideDetail(request, genbankId):
    try:
        sequence = herbiscide.objects.get(genbankId=genbankId)
    except Sequences.DoesNotExist:
        raise Http404("Sequence Does Not Exist")
    context = {"sequence": sequence}
    template = "search/herbicideDetail.html"
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
    # evalue of search given by user
    threshold = request.GET.get("threshold")
    if threshold != 'None':
        threshold = float(threshold)
    else:
        threshold = 0.0001
    database_name = request.GET.get("species")
    query = request.GET.get("query")
    query = validate_and_replace(query)
    print(query)
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

        cline_blast = NcbiblastnCommandline(query=query, db=database, evalue=threshold, outfmt=5, out=output)
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

    context = {"hits": hits, "searchQuery": {"species": database_name, "query": name}}
    template = "search/blastn_results.html"
    return render(request, template, context)


def blastn_at_ncbi(request, seq_id):
    query_sequence = Sequences.objects.get(pk=seq_id)
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence.sequence)
    blast_records = NCBIXML.read(result_handle)
    results = blast_records.alignments
    hits = []
    for result in results:
        for hit in result.hsps:
            hits.append({"Hit_exp": hit.expect, "Hit_query": hit.query,
                         "Hit_match": hit.match, "Hit_sbject": hit.sbjct,
                         "description": query_sequence.gene_description, "ID": query_sequence.id, "Hit_id": result.hit_id})
    context = {"hits": hits}
    template = "search/blastn_results.html"
    return render(request, template, context)


def blastn_render(request):
    return render(request, "search/blastn_search.html", {})


def validate_and_replace(sequence):
    sequence = list(sequence)
    sequence_list = ['A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N', 'Z']
    lower_seq = list(map(lambda a: a.lower(), sequence_list))
    for s in range(0, len(sequence)):
        if sequence[s] not in sequence_list:
            if sequence[s] not in lower_seq:
                sequence[s] = '#'
    # blank replacement
    new_sequence = []
    for s in range(0, len(sequence)):
        if sequence[s] != ' ':new_sequence.append(sequence[s])
    return ''.join(new_sequence)


# contact the admin of the site
def contact(request):
    form = ContactForm(request.POST or None)
    template = 'search/contact.html'
    message = None

    if form.is_valid():
        name = form.cleaned_data['name']
        affiliation = form.cleaned_data['affiliation']
        position = form.cleaned_data['position']
        comment = form.cleaned_data['comment']
        subject = 'Website Email'
        email_from = form.cleaned_data['email']
        email_to = [settings.EMAIL_HOST_USER]
        message = comment + '\n\nby \n\n' + name + '\nemail: '+email_from+'\naffiliation: '+affiliation+'\nposition: '+position

        # password obtained from environment variable
        # settings.EMAIL_HOST_PASSWORD = os.environ['Web_password']
        settings.EMAIL_HOST_PASSWORD = 'Herbs@123'

        send_mail(subject, message, email_from, email_to, fail_silently=True)
        message = 'Your message has been received, I will get back to you soon! Thanks!'
        form = ContactForm(request.GET)

    context = {'form': form, 'message': message}
    return render(request, template, context)


# protein search at ncbi
def blastx_at_ncbi(request, seq_id):
    query_sequence = Sequences.objects.get(pk=seq_id)
    result_handle = NCBIWWW.qblast("blastx", "nr", query_sequence.sequence)
    blast_records = NCBIXML.read(result_handle)
    results = blast_records.alignments
    hits = []
    for result in results:
        for hit in result.hsps:
            hits.append({"Hit_exp": hit.expect, "Hit_query": hit.query,
                         "Hit_match": hit.match, "Hit_sbject": hit.sbjct,
                         "description": query_sequence.gene_description, "ID": query_sequence.id, "Hit_id": result.hit_id})
    context = {"hits": hits}

    # change to new html page if needed later
    template = "search/blastn_results.html"
    return render(request, template, context)


def herbiscide_search(request):
    geneId= request.GET.get("geneId")
    genus = request.GET.get("genus")

    herbs = herbiscide.objects.filter(gene_id__iexact=geneId).filter(genus__iexact=genus)
    context = {"results": herbs, "searchQuery": {"geneId": geneId, "genus": genus}}

    template = "search/herbiscide_results.html"
    return render(request, template, context)


def herbiscide_render(request):
    template = "search/herbiscide_search.html"
    return render(request, template, {})

def download_herb_file(request, genbankId):
    sequence = herbiscide.objects.get(genbankId=genbankId)
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

    my_sequence = "> "+sequence.genbankId+" | "+sequence.description
    my_sequence = my_sequence+temp_sequence

    temp_file = tempfile.TemporaryFile()
    my_sequence = my_sequence.encode()

    temp_file.write(my_sequence)
    temp_file.seek(0)

    response = HttpResponse(temp_file, content_type='application/fasta;charset=UTF-8')
    response['Content-Disposition'] = "attachment; filename=%s" % sequence.organism+".fasta"

    temp_file.close()
    return response
