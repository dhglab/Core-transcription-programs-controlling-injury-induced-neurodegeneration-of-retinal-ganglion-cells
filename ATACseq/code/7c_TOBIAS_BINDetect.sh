#ATAC correct Tobias - subsample and then perform tn5 shift
#load python/3.7.1
#load samtools
#Yuyan 11-09-19


START=$SECONDS  

MOTIF=$1
INDIR=$2
GENOME=$3
PEAK=$4
OUTDIR=$5
CMP=$6

if [[ $CMP == nonTarget ]]; then
TOBIAS BINDetect --motifs $MOTIF \
                 --signals $INDIR/nonTarget_day{1,0,3}.unique.sorted.rmDup.merge.footprint.bw \
                 --genome $GENOME \
                 --peaks $PEAK \
                 --outdir $OUTDIR --time-series --cores 4
fi

if [[ $CMP == CTCFvsNon_day1 ]]; then
TOBIAS BINDetect --motifs $MOTIF \
                 --signals $INDIR/{nonTarget,CTCF}_day1.unique.sorted.rmDup.merge.footprint.bw \
                 --genome $GENOME \
                 --peaks $PEAK \
                 --outdir $OUTDIR --time-series --cores 4
fi

if [[ $CMP == CTCFvsNon_day3 ]]; then
TOBIAS BINDetect --motifs $MOTIF \
                 --signals $INDIR/{nonTarget,CTCF}_day3.unique.sorted.rmDup.merge.footprint.bw \
                 --genome $GENOME \
                 --peaks $PEAK \
                 --outdir $OUTDIR --time-series --cores 4
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

