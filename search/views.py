from Bio import pairwise2
from Bio.Seq import Seq
from django.http import Http404
from django.shortcuts import render
from .models import Sequences
from django.http import HttpResponse
import tempfile


def index(request):
    template = "search/index.html"
    context = {}
    return render(request, template, context)


def search_keyword(request):
    query = request.GET.get("query")
    species = request.GET.get("species")
    all_sequences = Sequences.objects.all()
    results = []
    for sequence in all_sequences:
        description = sequence.gene_description
        if description.find(query) != -1 and sequence.species == species:
            results.append(sequence)
    context = {"results": results}
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
    my_sequence = my_sequence+"\n"+temp_sequence

    temp_file = tempfile.TemporaryFile()
    my_sequence = my_sequence.encode()

    temp_file.write(my_sequence)
    temp_file.seek(0)

    response = HttpResponse(temp_file, content_type='application/fasta;charset=UTF-8')
    response['Content-Disposition'] = "attachment; filename=%s" % sequence.species+".fasta"

    temp_file.close()
    return response
