#Pillai lab
#Daniel Castaneda Mogollon

#This code reports a summary stats of fastq files, in regards of the length of the sequences, not the
#quality score.

import os
import statistics
path = input("Please input the path to your folder where you wish to analyze your samples: ")
os.chdir(path)
length_list=[]
counter=0
reads=0
out_file = open("Summary_output.csv","a")
out_file.write("Sample_name,Mean,StDev,Median,Reads,\n")
for item in os.listdir(path):
    if item[-6:]!='.fastq':
        continue
    else:
        print("Analyzing item: "+ item)
        with open(item) as f:
            content = f.readlines()
            for lines in content:
                #After the "@", it  counts the length of the lines
                if counter==1:
                    reads=reads+1
                    length_list.append(len(lines))
                    counter=0
                #Serves as a counter for the next .fastq line, the one with the sequence.
                if lines.startswith("@"):
                    counter = counter+1
        mean = statistics.mean(length_list)
        stdev = statistics.stdev(length_list)
        median = statistics.median(length_list)
        print("Mean: "+str(mean))
        print("Standard Deviation: "+str(stdev))
        print("Median: "+str(median))
        print("Reads: "+str(reads))
        out_file.write(item+','+str(mean)+','+str(stdev)+','+str(median)+','+str(reads)+',\n')
        reads=0

        #This makes sure the list gets empty again for the next file
        del(length_list[:])
        #print(length_list)

