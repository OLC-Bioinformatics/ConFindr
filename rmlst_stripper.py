from Bio import SeqIO

f = open('list.txt')
genes = f.readlines()
f.close()

for i in range(len(genes)):
    genes[i] = genes[i].rstrip()

rmlst_stuff = SeqIO.parse('/mnt/nas/Adam/assemblypipeline/rMLST/2017-03-29/rMLST_combined.fasta', 'fasta')
for item in rmlst_stuff:
    gene = item.id.split('_')[0]
    if gene not in genes:
        print('>' + item.id)
        print(str(item.seq))