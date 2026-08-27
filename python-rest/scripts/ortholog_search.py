import requests 
import json
from collections import Counter
import pandas as pd

file_name = "/Users/Dani/OneDrive/GenDataPro/TFs.txt"
with open(file_name) as gene_text:
    gene_file = gene_text.read()
gene_ids = gene_file.split()

for gene in gene_ids:    
    new_gene = gene
    response = requests.get("http://rest.ensembl.org/homology/id/{}?content-type=application/json".format(new_gene))
    
    json_data = response.content.decode("utf-8")
    data = json.loads(json_data)

    homologies = data["data"][0]["homologies"]

    orthologs = [ortho for ortho in homologies if ortho["type"] == "ortholog_one2one" or "ortholog_one2many"]

    for transcript in orthologs:
        if (
            transcript["target"]["species"] == "mus_musculus" or transcript["target"]["species"] == "drosophila_melanogaster" or
            transcript["target"]["species"] == "caenorhabditis_elegans" or transcript["target"]["species"] == "xenopus_tropicalis"
            ):
            #print transcript.get('source', {}).get('id'),
            print transcript.get('target', {}).get('id')
            #print transcript.get('target', {}).get('species')


    print transcript.get('source', {}).get('id')
