$(ANALYSIS_DIR)/%_realign_targets.list: $(ANALYSIS_DIR)/%_gdna_dedup.bam $(REFDIR)/d%_prepend.fasta.bwt $(REFDIR)/d%_prepend.dict
	gatk -T RealignerTargetCreator \
		-R $(REFDIR)/d$*_prepend.fasta \
		-fixMisencodedQuals \
		-I $< \
		-o $@

$(ANALYSIS_DIR)/%_realigned_reads.bam : $(ANALYSIS_DIR)/%_realign_targets.list $(ANALYSIS_DIR)/%_gdna_dedup.bam $(REFDIR)/d%_prepend.dict
	gatk -T IndelRealigner
		-R $(REFDIR)/d$*_prepend.fasta \
		-I $(ANALYSIS_DIR)/$*_gdna_dedup.bam \
		-fixMisencodedQuals \
		-targetIntervals $(ANALYSIS_DIR)/%_realign_targets.list \
		-o $@

$(ANALYSIS_DIR)/on_sim/%_bwa.bam: $(%) $(SIMFASTA2).bwt | $(ANALYSIS_DIR)
	bwa mem -M $(SIMFASTA2) \
		-R "@RG	ID:$*	SM:$*	PL:illumina	LB:lib1	PU:unit1" \
		-t 8 \
		$(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
		$(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) \
		| samtools view -bS -o $(basename $@)_unsorted.bam -
	samtools sort $(basename $@)_unsorted.bam $(basename $@)


$(ANALYSIS_DIR)/sec_gdna_mapped.bam: sequence/sec_gdna_reads $(SECFASTA2).bwt | $(ANALYSIS_DIR)
	bwa mem -M -t 10 $(SECFASTA2) \
		-R "@RG	ID:group2	SM:sample1	PL:illumina	LB:lib1	PU:unit1" \
		$<_1 $<_2 \
		| samtools view -bS -o $(basename $@)_unsorted.bam -;
	samtools sort -@ 4 $(basename $@)_unsorted.bam $(basename $@)
