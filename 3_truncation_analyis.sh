TRUNCATION_END=T

create_counts_file () {

# Extract mapped reads from sam file
samtools view -h -F 4 $1\_$2.sam > $1\_$2_mapped.sam
# Convert sam to bam
samtools view -b $1\_$2_mapped.sam > $1\_$2.bam
# Sort bam file, required for downstream analyses
samtools sort $1\_$2.bam > $1\_$2_sorted.bam
# Convert bam to bed
bedtools bamtobed -i $1\_$2_sorted.bam > $1\_$2.bed
# Generate text file showing the coverage depth at each base of reference
samtools depth -b $1\_$2.bed $1\_$2_sorted.bam > $1\_$2_depth.txt
# Remove headers from mapped sam file
samtools view -F 4 -F 16 $1\_$2.sam > $1\_$2_mapped.sam
# Record positions of reads that end in X at each position in sam file
while read line
do
	sequence=`echo "$line" | awk '{print $10}'`
	if echo "$sequence" | grep -q "^$TRUNCATION_END"
	then
		echo "$line" | awk '{print $4}' >> $1\_$2_counts.txt
	fi
done < $1\_$2_mapped.sam
# Count these positions and record the copies of each unique position
sort $1\_$2_counts.txt | uniq -c > tmp.counts
mv tmp.counts $1\_$2_counts.txt
}

# Preprocess experimental sam file
for FILE in *_CMC.sam
do
# Isolate prefix of file before "_" character
file_prefix=`echo $FILE | cut -d'_' -f1`
# Create file with contains counts of reads that end in X at each position.
create_counts_file $file_prefix CMC
create_counts_file $file_prefix CTL

# Extract mapped reads from sam file
#samtools view -h -F 4 $file_prefix\_CMC.sam > $file_prefix\_CMC_mapped.sam
# Convert sam to bam
#samtools view -b $file_prefix\_CMC_mapped.sam > $file_prefix\_CMC.bam
# Sort bam file, required for downstream analyses
#samtools sort $file_prefix\_CMC.bam > $file_prefix\_CMC_sorted.bam
# Convert bam to bed
#bedtools bamtobed -i $file_prefix\_CMC_sorted.bam > $file_prefix\_CMC.bed
# Generate text file showing the coverage depth at each base of reference
#samtools depth -b $file_prefix\_CMC.bed $file_prefix\_CMC_sorted.bam > $file_prefix\_CMC_depth.txt
# Count the number of reads that end in X at each position
#samtools view -F 4 $file_prefix\_CMC.sam > $file_prefix\_CMC_mapped.sam
#while read line
#do
#sequence=`echo "$line" | awk '{print $10}'`
#if echo "$sequence" | grep -q "^$TRUNCATION_END"
#then
#echo "$line" | awk '{print $4}' >> $file_prefix\_CMC_counts.txt
#fi
#done < $file_prefix\_CMC_mapped.sam
#sort $file_prefix\_CMC_counts.txt | uniq -c > tmp.counts
#mv tmp.counts $file_prefix\_CMC_counts.txt

# Repeat these steps for the control
#samtools view -h -F 4 $file_prefix\_CTL.sam > $file_prefix\_CTL_mapped.sam
#samtools view -b $file_prefix\_CTL_mapped.sam > $file_prefix\_CTL.bam
#samtools sort $file_prefix\_CTL.bam > $file_prefix\_CTL_sorted.bam
#bedtools bamtobed -i $file_prefix\_CTL_sorted.bam > $file_prefix\_CTL.bed
#samtools depth -b $file_prefix\_CTL.bed $file_prefix\_CTL_sorted.bam > $file_prefix\_CTL_depth.txt
#samtools view -F 4 $file_prefix\_CTL.sam > $file_prefix\_CTL_mapped.sam
#while read line
#do
#sequence=`echo "$line" | awk '{print $10}'`
#if echo "$sequence" | grep -q "^$TRUNCATION_END"
#then
#echo "$line" | awk '{print $4}' >> $file_prefix\_CTL_counts.txt
#fi
#done < $file_prefix\_CTL_mapped.sam
#sort $file_prefix\_CTL_counts.txt | uniq -c > tmp.counts
#mv tmp.counts $file_prefix\_CTL_counts.txt

# Run the following commands on each line of the experimental T counts file
while read line
do
# Isolate number of reads ending at X for current site
CMC_X_sites=`echo "$line" | awk '{print $1}'`
# Isolate current position on reference
CMC_position=`echo "$line" | awk '{print $2}'`
# Find the current position in the file containing coverage depth and pull out the depth at that position
CMC_depth=`awk -v pos="$CMC_position" '$2 == pos {print $3}' $file_prefix\_CMC_depth.txt`
# Calculate number of reads ending in X base divided by total coverage at the current position
CMC_mean=$(bc <<< "scale=5 ; $CMC_X_sites / $CMC_depth")
echo "$CMC_position $CMC_mean" >> CMC_ratios.txt
done < $file_prefix\_CMC_counts.txt

# Repeat these steps for the control
while read line
do
CTL_X_sites=`echo "$line" | awk '{print $1}'`
CTL_position=`echo "$line" | awk '{print $2}'`
CTL_depth=`awk -v pos="$CTL_position" '$2 == pos {print $3}' $file_prefix\_CTL_depth.txt`
CTL_mean=$(bc <<< "scale=5 ; $CTL_X_sites / $CTL_depth")
echo "$CTL_position $CTL_mean" >> CTL_ratios.txt
done < $file_prefix\_CTL_counts.txt

decimal_regex='[0-9]*[.][0-9]+'
while read line
do
position=`echo "$line" | awk '{print $1}'`
CMC_mean=`echo "$line" | awk '{print $2}'`
CTL_mean=`awk -v pos="$position" '$1 == pos {print $2}' CTL_ratios.txt`
if ! [[ $CTL_mean =~ $decimal_regex ]]
then
echo "$position CTL_0" >> $file_prefix\_results.txt
else
ratio=$(bc <<< "scale=5 ; $CMC_mean / $CTL_mean")
log=`echo "l($ratio)/l(2)" | bc -l`
echo "$position $log" >> $file_prefix\_results.txt
fi
done < CMC_ratios.txt

while read line
do
position=`echo "$line" | awk '{print $1}'`
CTL_mean=`echo "$line" | awk '{print $2}'`
CMC_mean=`awk -v pos="$position" '$1 == pos {print $2}' CMC_ratios.txt`
if ! [[ $CMC_mean =~ $decimal_regex ]]
then
echo "$position CMC_0" >> $file_prefix\_results.txt
else
ratio=$(bc <<< "scale=5 ; $CMC_mean / $CTL_mean")
log=`echo "l($ratio)/l(2)" | bc -l`
echo "$position $log" >> $file_prefix\_results.txt
fi
done < CTL_ratios.txt

sort -n $file_prefix\_results.txt | uniq > tmp.txt
mv tmp.txt $file_prefix\_results.txt

rm CTL_ratios.txt
rm CMC_ratios.txt

done

create_counts_file () {

# Extract mapped reads from sam file
samtools view -h -F 4 $file_prefix\_$1.sam > $file_prefix\_$1_mapped.sam
# Convert sam to bam
samtools view -b $file_prefix\_$1_mapped.sam > $file_prefix\_$1.bam
# Sort bam file, required for downstream analyses
samtools sort $file_prefix\_$1.bam > $file_prefix\_$1_sorted.bam
# Convert bam to bed
bedtools bamtobed -i $file_prefix\_$1_sorted.bam > $file_prefix\_$1.bed
# Generate text file showing the coverage depth at each base of reference
samtools depth -b $file_prefix\_$1.bed $file_prefix\_$1_sorted.bam > $file_prefix\_$1_depth.txt
# Count the number of reads that end in X at each position
samtools view -F 4 $file_prefix\_$1.sam > $file_prefix\_$1_mapped.sam
while read line
do
sequence=`echo "$line" | awk '{print $10}'`
if echo "$sequence" | grep -q "^$TRUNCATION_END"
then
echo "$line" | awk '{print $4}' >> $file_prefix\_$1_counts.txt
fi
done < $file_prefix\_$1_mapped.sam
sort $file_prefix\_$1_counts.txt | uniq -c > tmp.counts
mv tmp.counts $file_prefix\_$1_counts.txt
}
