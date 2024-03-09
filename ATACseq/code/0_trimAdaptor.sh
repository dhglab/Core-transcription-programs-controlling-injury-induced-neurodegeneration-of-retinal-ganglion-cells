#Trim adaptor and filter low-quality reads with TrimGalore
#Yuyan Chenng
#10-03-19

START=$SECONDS

INDIR=$1
OUTDIR=$2
trim_galore=~/.local/bin/trim_galore

cd $OUTDIR/FASTQ
for fastq1 in $INDIR/*R1.fastq.gz;do
    fastq2=${fastq1/R1.fastq.gz/R2.fastq.gz}
    name_prefix=$(basename "$fastq1")
    name_prefix=${name_prefix/.fastq.gz/}
    trim_fastq1=$OUTDIR/FASTQ/$name_prefix'_val_1.fq.gz'
    
    if [[ ! -e "$trim_fastq1" ]]; then
    $trim_galore --stringency 5 --paired --nextera \
    --fastqc_args "--outdir $OUTDIR/fastQC" $fastq1 $fastq2
    fi
done



sleep 10
# hack to ensure job lasts 10 minutes to ensure no throttling
END=$SECONDS
ELAPSED=$((END-START))
echo $ELAPSED
if [ $ELAPSED -lt 600 ]; then
  TOSLEEP=$((600 - ELAPSED))
  sleep $TOSLEEP
fi
