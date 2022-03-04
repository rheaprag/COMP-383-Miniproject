#This is my mini-project code

#import necessary libraries
import os
import datetime
from Bio import SeqIO

#write paths
directory = os.popen('pwd'). read().rstrip()
path = (directory + '/Users/rhea/PycharmProjects/COMP383_MiniProject')
os.system('mkdir' + path)
os.chdir(path)

#1. Retrieve the Illumina reads for the resequencing of K-12 project
#get files from sra - https://www.ncbi.nlm.nih.gov/sra/SRX5005282
os.system('wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR818/SRR8185310/SRR8185310.sra')
os.system('fastq-dump -I --split-files SRR8185310.sra')

#2. Use SPAdes to assemble the genome
#http://cab.spbu.ru/software/spades/
os.system('spades -k 55,77,99,127 -t 2 --only-assembler -s SRR8185310_1.fastq -o ' + path)
command = 'spades -k 55,77,99,127 -t 2 --only-assembler -s SRR8185310_1.fastq -o ' + path

#write to log file
log_file = open(path + "miniproject.log", 'w')
for i in range(1):
    log_file.write("Spades:" + command + "\n" + "\n")

#3. Calculate the number of contigs with length > 1000
#"There are # contigs > 1000 in the assembly"
#Only consider contigs > 1000 bp in length
sequences = []
#open and read fasta file, if the length of the sequence is > 1000, add record to list
with open(path + "contigs.fasta", 'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > 1000:
            sequences.append(record)

#create list for sequences with length > 1000
sequences_length = str(sequences)

assembly = ("There are " + sequences_length + "contigs > 1000 in the assembly.")
print(assembly)
#write to fasta file
SeqIO.write(sequences, "sequences.fasta", "fasta")

#write to lof gile
log_file = open(path + "miniproject.log", 'a')
for i in range(1):
    log_file.write(assembly + "\n" + "\n")

#4. Calculate the length of the assembly
#"There are # bp in the assembly"
length = 0
temp = 0
for record in SeqIO.parse("sequences.fasta", "fasta"):
    temp = len(record)
    length += temp
final_length = str(length)
total_length = ("There are " + final_length + "bp in the assembly")

#write to log file
log_file = open(path + "miniproject.log", 'a')
for i in range(1):
    log_file.write(total_length + "\n" + "\n")

#5. Use Prokka to annotate the assembly, use Escherichia genus database
os.system('prokka --force --outdir ' + path + 'Prokka_Output/ --genus Escherichia --locustag ECOL long_sequences.fasta')
prokka = "prokka --outdir " + path + "Prokka_Output/ --genus Escherichia --locustag ECOL long_sequences.fasta"

#write to log file
log_file = open(path + "miniproject.log", 'a')
for i in range(1):
    log_file.write("\n" + prokka + "\n" + "\n")

#6. Write the results of the annotation in the *.txt file to log file in the same format as the *.txt file
current = datetime.datatine.now()
hold = str(current)

month = (current.month)
day = ("%d" % current.day)
year = ("%d" % current.year)

if month >= 10:
    current_date = (month + day + year)
else:
    month = ("" + str(month))
    current_date = (month + day + year)

txt_file = "PROKKA_" + current_date + ".txt"

with open(path + 'Prokka_Output/' + txt_file) as hold:
    with open(path + "miniproject.log", 'a') as only:
        for x in hold:
            only.write(x)

#7. Write to the log file if a discrepency is found
#uhhhhhh
file_data = open(path + "Prokka_Output/" + txt_file)
CDS = 4140
TRNA = 89

for x in txt_file:
    if x.startswith("CDS"):
        cds_line = x[5:]
    if x.startswith("tRNA"):
        trna_line = x[5:]

cds_line = int(cds_line)
cds_count = CDS - cds_line

trna_line = int(trna_line)
trna_count = TRNA - trna_line

if cds_count > 0:
    x = "less"
elif cds_count < 0:
    x = "additional"
    cds_count = abs(cds_count)
if trna_count > 0:
    n = "additional"
    trna_count = abs(trna_count)

final = ("Prokka found " + str(cds_count) + " " + x + " CDS and " + str(trna_count) + " " + n + " tRNA than the RefSeq.")

# Write to log file
append = open(path + "/OptionA.log", "a")
for n in range(1):
    append.write(final)
