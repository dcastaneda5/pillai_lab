#Pillai Lab
#Daniel Castaneda Mogollon

#This code reports the number of snps, indels, and divergence from a reference sequence against query sequences.
#It takes a reference sequence from a .fasta file and gets the difference from every other sequence. The reference
#sequence must have 'reference' in its header. It compares nucleotides A,U,G,C,T,N and indels as '-'. The input file
#must be previously aligned and be in a .fasta format. It only asks the user for the path where the file is and a name
#for the output file.

print("This code has the purpose of printing the number of snps found between a reference vs a query sequence(s).")
print("The user must input an aligned file of said sequences, and label the reference as '>reference'.")

import os

headers=[]
sequences=[]
nucleotides=[]
reference_header=''
reference_sequence=''
#opening the MSA aligned file
file_analyze = input("Please type the path to your aligned file (must be in .txt or .fasta format): ")
with open(file_analyze) as f:
    content = f.readlines()
    for lines in content:
        if lines.startswith(">"):
            headers.append(lines.replace("\n",""))
        else:
            sequences.append(lines.upper().replace("\n",""))
i=0
#getting the reference sequence, assuming it has the word 'reference' in it
for item in headers:
    if item.__contains__("reference"):
        # getting the reference sequence from the header index
        reference_sequence = sequences[i]
        reference_header = headers[i]
    else:
        i=i+1
#making sure the sequences have the same size (considering indels as -)
for item in sequences:
    nucleotides.append(len(item))
first_sequence_length = nucleotides[0]
for items in nucleotides:
    if len(item)!=first_sequence_length:
        print("The sequences are different size. Exiting the program now.")
        exit()
    else:
        continue

print("Sequences are the same size. Analyzing snps . . .")
print(reference_sequence)
#getting an output file as csv
f_out_name = input("Please name your output file: ")
f_out = open(f_out_name+'.csv','a')
snps_list=[]
indels_list=[]
DNA_nucleotides=['A','T','G','C','N','U']
snp_counter = 0
indel_counter = 0
for sequence in sequences:
    for j in range(0,len(reference_sequence),1):
        if (reference_sequence[j] in DNA_nucleotides) and (sequence[j] in DNA_nucleotides):
            if reference_sequence[j]!=sequence[j]:
                snp_counter = snp_counter+1
            else:
                #in this case the nucleotides match each other
                continue
        else:
            if (reference_sequence[j]=='-') and (sequence[j]=='-'):
                continue
            elif (reference_sequence[j]=='-' and sequence[j]!='-'):
                indel_counter= indel_counter+1
            elif (reference_sequence[j]!='-' and sequence[j]=='-'):
                indel_counter= indel_counter+1
    snps_list.append(snp_counter)
    indels_list.append(indel_counter)
    snp_counter = 0
    indel_counter=0
#printing the snps and indels into the output file
f_out.write("Sequence name,SNPs,Indels,Divergence(%),\n")
divergence_list=[]
for m in range(0,len(headers),1):
    divergence_list.append(((snps_list[m]+indels_list[m])/len(sequences[m])*100))
    f_out.write(headers[m]+','+str(snps_list[m])+','+str(indels_list[m])+','+str(divergence_list[m])+',\n')














