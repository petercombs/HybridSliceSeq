DEXSEQ_PY=/home/pcombs/R/x86_64-pc-linux-gnu-library/3.3/DEXSeq/python_scripts


%.dexseq_counts: %.bam.bai
	./qsubber \
		--job-name $(notdir $(@D))_$(*F) \
		--threads 1 \
	python2 $(DEXSEQ_PY)/dexseq_count.py \
		--paired=yes \
		--stranded=no \
		--format=bam \
		--order=pos \
		/home/pcombs/HybridSliceSeq/Reference/dexseq_mel.gtf \
		$*.bam \
		$@

