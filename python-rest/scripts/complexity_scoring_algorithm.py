import time
try:
    from queue import Queue
except ImportError:
    from Queue import Queue
from threading import Thread
import requests 
import json
from collections import Counter
import pandas as pd



isoforms = {}
errors = {key:[] for key in ["paralog", "isoform", "motif"]}
file_name = "ortho_tf_9174.txt"
with open(file_name) as gene_text:
    gene_file = gene_text.read()
all_gene_ids = gene_file.split()#[:200]
ts = time.time()

hd_out = pd.HDFStore("out.h5")

try:
    frame = hd_out.get("gene_data")
    gene_ids = [gene_id for gene_id in all_gene_ids if gene_id not in hd_out.get("genes").index]

except KeyError:
    gene_ids = all_gene_ids
    columns = [
        "isoforms",
        "paralogues",
        "motifs",
        "complexity",
        "name",
        "species",
    ]

    frame = pd.DataFrame(columns=columns, data=0, index=gene_ids)
    hd_out.put("gene_data", frame, format='table', data_columns=True)

paralogues = Queue()
isoforms = Queue()
motifs = Queue()

class IsoformWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
    
    def run(self):
        while True:
            gene_id = self.queue.get()
            response_isoform = requests.get("http://rest.ensembl.org/lookup/id/{}?content-type=application/json;expand=2".format(gene_id))
            try:
                json_data_isoform = response_isoform.content.decode("utf-8")
                data_isoform = json.loads(json_data_isoform)
                count_isoform = Counter([transcript["biotype"] for transcript in data_isoform["Transcript"]])
                frame.ix[gene_id, "isoforms"] = count_isoform["protein_coding"]
                frame.ix[gene_id, "complexity"] += count_isoform["protein_coding"]
                frame.ix[gene_id, "name"] = data_isoform["display_name"]
                frame.ix[gene_id, "species"] = data_isoform["species"]
                for transcript in data_isoform["Transcript"]:
                    if transcript["biotype"] == "protein_coding":
                        motifs.put((gene_id, transcript["Translation"].get("id")))
            except Exception as e:
                print(e)
                errors["isoform"].append(gene_id)
            self.queue.task_done()
            time.sleep(1)

class ParalogWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
    
    def run(self):
        while True:
            gene_id = self.queue.get()
            response_paralog = requests.get("http://rest.ensembl.org/homology/id/{}?content-type=application/json".format(gene_id))
            json_data_paralog = response_paralog.content.decode("utf-8")
            data_paralog = json.loads(json_data_paralog)
            try:
                json_data_paralog = response_paralog.content.decode("utf-8")
                data_paralog = json.loads(json_data_paralog)
                count_paralog = Counter([homo["type"] for homo in data_paralog["data"][0]["homologies"]])
                frame.ix[gene_id, "paralogues"] = count_paralog["within_species_paralog"] + 1
                frame.ix[gene_id, "complexity"] += count_paralog["within_species_paralog"] + 1
            except Exception as e:
                print(e)
                errors["paralog"].append(gene_id)
            self.queue.task_done()
            time.sleep(1)
            
class MotifWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
    
    def run(self):
        while True:
            gene_id, translation = self.queue.get()
            response_motif = requests.get("http://rest.ensembl.org/overlap/translation/{}?content-type=application/json;type=Prosite_profiles".format(translation))
            try:
                json_data_motif = response_motif.content.decode("utf-8")
                data_motif = json.loads(json_data_motif)
                frame.ix[gene_id, "motifs"] += len(data_motif)
                frame.ix[gene_id, "complexity"] += len(data_motif)
            except Exception as e:
                print(e)
                errors["motif"].append([gene_id, translation])
            self.queue.task_done()
            time.sleep(1)

class HDF5WriterThread(Thread):
    def run(self):
        while not motifs.empty() and not paralogues.empty() and not isoforms.empty():
            hd_out.append("gene_data", frame, format='table', data_columns=True)
            time.sleep(30)
threads = []
for _ in range(10):
    isoform_worker = IsoformWorker(isoforms)
    isoform_worker.daemon = True
    isoform_worker.start()
    paralog_worker = ParalogWorker(paralogues)
    paralog_worker.daemon = True
    paralog_worker.start()
    motif_worker = MotifWorker(motifs)
    motif_worker.daemon = True
    motif_worker.start()


for gene in gene_ids:
    isoforms.put(gene)
    paralogues.put(gene)

writerthread = HDF5WriterThread()
writerthread.daemon = True
writerthread.start()

isoforms.join()
paralogues.join()
motifs.join()

print("Took {}".format(time.time() - ts))
frame
