#Pipeline for Pseudo analysis. Run this script in a directory containing pairs of fastq files, gzipped or uncompressed.
#Names should be of the format NAME_1.fastq NAME_2.fastq. Please remove any underscore characters except for before the final 1 or 2 in the filename.
trimmomatic=/export/home/chris/Packages/Trimmomatic-0.38/trimmomatic-0.38.jar
PAN_INDEX=PAN

#To unzip fastq.gz files
#gunzip *.fastq.gz

#Loop for each set of fastq files.
for FILE in *.fastq
do

#Isolate the base name of the fastq files.
name=`echo $FILE | cut -d"." -f1`

#Pair end remove adapters, low quality and N bases sliding window and less than 36b long
java -jar $trimmomatic SE -phred33 $FILE $name\_trimmed.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Trimming adapter from either end of reads.
cutadapt -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o $name\_clean.fastq $name\_trimmed.fastq

#Mapping Reads
bowtie2 -N 1 -t -x $PAN_INDEX -U $name\_clean.fastq -S $name.sam

#Remove any reads that map to reverse strand, and remove any that don't map.
samtools view -F 20 $name.sam > $name\_mapped.sam

#Count the reads start with T at the 5'end reflects the stops at the 3' end cDNA
awk '$10 ~ /^T/' $name\_mapped.sam | awk '{print $4}' | sort | uniq -c | sort -n -r > $name.txt



done
