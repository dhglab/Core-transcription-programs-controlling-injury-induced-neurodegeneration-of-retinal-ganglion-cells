#!/bin/bash
# Shell script to map atac-seq fastq files onto mm10 genome using Bowtie2 and mark duplicate using picard tools
# Yuyan Cheng
# 10-22-2021

START=$SECONDS

INDIR=$1
OUTDIR=$2
INDEX=$3
GENOME=$4

for fastq in $INDIR/*R1.fastq.gz;do
	fastq2=${fastq/R1.fastq.gz/R2.fastq.gz}
	
	name_prefix=${fastq/_R1.fastq.gz/}
	name_prefix=${name_prefix//$INDIR/$OUTDIR}
	
	sam=$name_prefix'.sam'
	bam=$name_prefix'.bam'
	sortedname=$name_prefix'.sorted'
	if [ ! -e "$sortedname.markDup.bam" ]; then
		if [ ! -e "$sortedname.bam" ]; then
			if [ ! -e "$sam" ]; then
				bowtie2 -p 6 --very-sensitive -X 2000 -x $INDEX -1 $fastq -2 $fastq2 > $sam
			fi
			if [ ! -e "$bam" ]; then
				samtools view -bS $sam > $bam # convert SAM to BAM
			fi
			samtools sort -O bam -o $sortedname.bam -T temp $bam # Sort BAM file
			samtools index $sortedname.bam # Index sorted BAM file
			rm $sam
			rm $bam
		fi
		## Picard: MarkDuplicates, CollectAlignmentSummaryMetrics, CollectGcBiasMetrics	
		java -Xmx4G -jar $PICARD CollectAlignmentSummaryMetrics R=$GENOME I=$sortedname.bam O=$sortedname.metricsAlign
		java -Xmx4G -jar $PICARD CollectGcBiasMetrics R=$GENOME I=$sortedname.bam O=$sortedname.metricsGCFull S=$sortedname.metricsGCSum CHART=$sortedname.metricsGC.pdf
		java -Xmx4G -jar $PICARD MarkDuplicates I=$sortedname.bam O=$sortedname.markDup.bam M=$sortedname.metricsDup
		samtools index $sortedname.markDup.bam
		if [ -e "$sortedname.markDup.bam" ]; then
			rm $sortedname.bam
			rm $sortedname.bam.bai
		fi
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