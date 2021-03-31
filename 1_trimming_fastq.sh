#Batch script for cleaning reads in fastq files
#Run from within a directory that contains fastq files of the form NAME_1.fastq and NAME_2.fastq
#Remove any extra "_" or "." symbols from NAME

for FILE in *_1.fastq.gz
do
file_prefix=`echo $FILE | cut -d"_" -f1`
cutadapt -j 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 25,25 -o $file_prefix\_1_cleaned.fastq -p $file_prefix\_2_cleaned.fastq $file_prefix\_1.fastq $file_prefix\_2.fastq
mv $file_prefix\_1_cleaned.fastq $file_prefix\_1.fastq
mv $file_prefix\_2_cleaned.fastq $file_prefix\_2.fastq
done

#Adding mapping with tophat2
 tophat2 --no-novel-indels INDEX SM-16_1_clean.fastq.gz SM-16_2_clean.fastq.gz
