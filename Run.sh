#Cut adaptor option suggested by Demian, worked for CONTROLS
cutadapt -j 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 25,25 -o trimm_1_cleaned.fastq -p trimm_2_cleaned.fastq 1.fastq 2.fastq

#Cut adaptor indicated by RAP-AMT protocol
cutadapt -a CACGACGCTCTTCCGATCT --overlap 3 --minimum-length 100  -o trimmed_T1-CMC_S3.fastq T1-CMC_S3.fastq

#For the reverse adapter removal in the no 2 reads
#cutadapt -A AGATCGGAAGAGCGTCGTG --overlap 3 --minimum-length 100  -o trimmed_T1-CMC_S3.fastq T1-CMC_S3.fastq

#For pair-end reads
cutadapt -a CACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTG -o trimmed_T1-CMC_1.fastq -p trimmed_T1-CMC_2.fastq T1-CMC_S3_L001_R1_001.fastq T1-CMC_S3_L001_R2_001.fastq


#Cutting with trimmomatic for quality control
java -jar /home/chris/bin/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 T1-CMC_S3_L001_R1_001.fastq T1-CMC_S3_L001_R2_001.fastq ofp.fq ofu.fq orp.fq oru.fq ILLUMINACLIP:/home/chris/bin/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

#TruSeq2-PE.fa has the adapter we are using to make the libraries
#Universal Primer in 5'->3': AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#grep 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' ~/Downloads/trimmomatic-0.38/adapters/*.fa
#/home/chris/Downloads/trimmomatic-0.38/adapters/TruSeq2-PE.fa:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#/home/chris/Downloads/trimmomatic-0.38/adapters/TruSeq2-PE.fa:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT


#I tried:long adaptor Universal complete and partial, none of them reduce the peak of the lenght

#FASTQC analysis for Quality Standards
fastqc trimmed_T1-CMC_S3.fastq 

#For mapping from here:
bowtie2 -t -x ../Pseudo-Seq/PAN_noK7 -U trimmed_T1-CMC_S3.fastq -S mapped_trimmed_T1-CMC_S3.sam

#For the controls I used tophat 2.1.1
#I build the index as CONTROL with 
bowtie-build CONTROL.fa CONTROL
tophat2 --no-novel-indels CONTROL SM-16_1_clean.fastq.gz SM-16_2_clean.fastq.gz 
mv tophat_out/ SM-16/
#samtools view -Sb mapped_trimmed_T1-CMC_S3.sam > mapped_trimmed_T1-CMC_S3.bam
#samtools view -Sb SM-16/accepted_hits.bam
samtools view -h SM-16/accepted_hits.bam > SM-16.sam
cat SM-16.sam | awk '{print $4}' > SM-16.tex

#Batch script to align multiple sets of fastq files to a reference with bwa
#Run from within a directory that contains fastq files of the form NAME_1.fastq and NAME_2.fastq
#Remove any extra "_" or "." symbols from NAME
#Additionally, you will need a reference fasta file and supply it to the script as you call it
#Usage: batch_bwa_alignment.sh <name_of_reference.fasta>

#Comment out one of the following two index lines depending on your reference.
#Index for long reference, like a genome. Comment out if your reference is short
#bwa index -a bwtsw -p index $1
#Index for short reference, like a single locus. Comment out if your reference is long
bwa index -a is -p index $1

for FILE in *_1.fastq
do
name=`echo $FILE | cut -d"_" -f1`
bwa mem -t 5 index $name\_1.fastq $name\_2.fastq > $name.sam
done


#For Sam to Bam
samtools view -Sb mapped_trimmed_T1-CMC_S3.sam > mapped_trimmed_T1-CMC_S3.bam

#Create a file with only column contains information POS (Position Number Start)
cat mapped_trimmed_T1-CMC_S3.sam | awk '{print $4}' > temporary02.tex 

#Extract and count only Ts
grep '^T' temporary02.tex | wc -l
#count all lines
wc -l temporary02.tex

##############################################################
Col	Field	Type	Brief Description
1	QNAME	String	Query template NAME
2	FLAG	Int	bitwise FLAG
3	RNAME	String	References sequence NAME
4	POS	Int	1- based leftmost mapping POSition
5	MAPQ	Int	MAPping Quality
6	CIGAR	String	CIGAR String
7	RNEXT	String	Ref. name of the mate/next read
8	PNEXT	Int	Position of the mate/next read
9	TLEN	Int	observed Template LENgth
10	SEQ	String	segment SEQuence
11	QUAL	String	ASCII of Phred-scaled base QUALity+33
##############################################################

for FILE in *_results.txt
do
name=`echo $FILE | cut -d"_" -f1`
	for (( i=1; i<1078; i++ ))
	do
		found=0
		while read line
		do
			position=`echo "$line" | cut -d" " -f1`
			if [[ $position = $i ]]
			then
				found=1
				break
			fi
		done < $FILE
		if [[ $found = 0 ]]
		then
			depth_one=0
			depth_two=0
			while read depth_line
			do
				depth_position=`echo $depth_line | cut -d" " -f2`
				if [[ $i = $depth_position ]]
				then
					depth_one=1
					break
				fi
			done < $name\_CMC_depth.txt
			while read depth_line
			do
				depth_position=`echo $depth_line | cut -d" " -f2`
				if [[ $i = $depth_position ]]
				then
					depth_two=1
					break
				fi
			done < $name\_CTL_depth.txt
			if [[ $depth_one = 1 ]] && [[ $depth_two = 1 ]]
			then
				echo
				echo "$i 0" >> $name\_results.txt
			else
				echo "$i N/A" >> $name\_results.txt
			fi
		fi
	done
done
