# *************
# Generate Reads in Consensus Peaks for DESeq2 
# Yuyan Cheng 
# 01-12-20
# *************

START=$SECONDS

INDIR=$1
OUTDIR=$2
BAMDIR=$3
N=$4

fc_dir=/u/home/y/ycheng41/apps/subread-2.0.0-Linux-x86_64/bin

mergedPeak=$OUTDIR'/merged.narrowPeak'
consensusPeak=$OUTDIR'/consensus.narrowPeak'

#bedtools merge get ALL peaks (overlapped at least 1 bp; or non-overlapped) to one consensus file; 

cat $INDIR/*.narrowPeak_rmBlk | sort -k1,1 -k2,2n | bedtools merge -i - -d -1 > $mergedPeak


bedtools intersect -a $mergedPeak -b $INDIR/*.narrowPeak_rmBlk \
-sorted -wa -wb | cut -f 1-3 | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' - > $consensusPeak'.all'

# remove the peaks that appear in less than 20% of all samples (2 in this experiment)
awk -v t=$N '($4 >= t)' $consensusPeak'.all' > $consensusPeak

# Make SAF file (+1 because SAF is 1-based, BED/narrowPeak is 0-based)
awk 'OFS="\t" {print $1"."$2+1"."$3, $1, $2+1, $3, "."}' $consensusPeak > $consensusPeak'.saf'

## Then run fc: 
consensusPeakMatrix=$OUTDIR'/peaks_countMatrix.txt'
$fc_dir/featureCounts -a $consensusPeak'.saf' -F SAF -T 4 -p -o $consensusPeakMatrix $BAMDIR/*.unique.sorted.rmDup.bam  
#-T = number of threads; -p indicates paired-end data

if [[ -e "$consensusPeakMatrix" ]]; then
    rm $consensusPeak'.saf' $consensusPeak $mergedPeak
fi

sleep 10
# hack to ensure job lasts 10 minutes to ensure no throttling
END=$SECONDS
ELAPSED=$((END-START))
echo $ELAPSED
if [ $ELAPSED -lt 600 ]; then
  TOSLEEP=$((600 - ELAPSED))
  sleep $TOSLEEP
fi