#***********
# MACS2 for ATAC-seq calling narrow and broad peaks
# Yuyan 12-29-19
# shifted tagAlign have more peaks identifed than non-shifted BAM
#***********

START=$SECONDS

INDIR=$1
OUTDIR=$2



MACS2_path=~/.local/bin/macs2

for tagAlign in $INDIR/*.unique.sorted.rmDup.bam; do
    name_prefix=$(basename "$tagAlign")
    name_prefix=${name_prefix/.unique.sorted.rmDup.bam/}
    
    #**************
    # derive the narrowPeaks
    #**************
    
    narrow_peak=$OUTDIR'/'$name_prefix'.macs2_peaks.narrowPeak'
    if [[ ! -f $narrow_peak ]]; then
	   MACS2_cmd='macs2 callpeak -t '$tagAlign' -f BAMPE -g mm --nolambda --call-summits --keep-dup all --min-length 100 --q 0.05 -n '$name_prefix'.macs2 --outdir' $OUTDIR
	   # min-length was recommended for atac-seq to save false-negative as the defalt will be fragment length (usually too big)
	   # -SPMR MACS will SAVE signal per million reads
	   eval $MACS2_cmd

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