from Bio import SeqIO

lncrna = {}
list_lncrna = list(SeqIO.parse("./data/gencode.v33.lncRNA_transcripts.fa",
                               format="fasta"))
for x in list_lncrna:
    id = str(x.id).split("|")[0]
    seq = str(x.seq).replace("U", "T")
    if len(seq) > 200:
        lncrna[id] = seq

f = open("data/gencode.v33.lncRNA_transcripts_new.fa","w")
for k,v in lncrna.items():
    f.write(">"+k+"\n")
    f.write(v+"\n")
f.close()