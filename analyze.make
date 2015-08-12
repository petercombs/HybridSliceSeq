include gmsl

comma := ,
space := $(null) #

wittkopp_reads:  sequence/SRR869573 sequence/SRR869575 sequence/SRR869590 sequence/SRR869596 sequence/SRR869597 sequence/SRR869598 sequence/SRR869599 sequence/SRR869600 sequence/SRR869601 sequence/SRR869602 sequence/SRR869844 sequence/SRR869603 sequence/SRR869905 sequence/SRR869592 sequence/SRR869579 sequence/SRR869580 sequence/SRR869587
	echo "All Downloaded"

sim_gdna = sequence/SRR869579 sequence/SRR869580
sec_gdna = sequence/SRR869587
simXsec = sequence/SRR869602 sequence/SRR869844
secXsim = sequence/SRR869603 sequence/SRR869905

sequence/%:
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_1.fastq.gz
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_2.fastq.gz
	touch $@


$(ANALYSIS_DIR)/on_sim/%_on_sim_bowtie.bam: $(%) $(basename $(SIMFASTA2)).1.ebwt | $(ANALYSIS_DIR)
	bowtie \
		--try-hard \
		-p 8 \
		--sam \
		--sam-RG "ID:$*" \
		--sam-RG "SM:$*" \
		--sam-RG "PL:illumina" \
		--sam-RG "LB:lib1"\
		--sam-RG "PU:unit1" \
		$(basename $(SIMFASTA2)) \
		-1 $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
		-2 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) \
		| samtools view -bS -o $(basename $@)_unsorted.bam -
	samtools sort $(basename $@)_unsorted.bam $(basename $@)

$(ANALYSIS_DIR)/on_sim/%_on_sim_bowtie2.bam: $(%) $(basename $(SIMFASTA2)).1.bt2 | $(ANALYSIS_DIR)
	bowtie2 \
		--very-sensitive \
		-p 8 \
		--rg-id $* \
		--rg "SM:$*" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x $(basename $(SIMFASTA2)) \
		-1 $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
		-2 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) \
		| samtools view -bS -o $(basename $@)_unsorted.bam -
	samtools sort $(basename $@)_unsorted.bam $(basename $@)

$(ANALYSIS_DIR)/on_sec/%_on_sec_bowtie2.bam: $(%) $(basename $(SECFASTA2)).1.bt2 | $(ANALYSIS_DIR)
	bowtie2 \
		--very-sensitive \
		-p 8 \
		--rg-id $* \
		--rg "SM:$*" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x $(basename $(SECFASTA2)) \
		-1 $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
		-2 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) \
		\| samtools view -bS -o $(basename $@)_unsorted.bam -
	samtools sort -@ 4 $(basename $@)_unsorted.bam $(basename $@)
	rm $(basename $@)_unsorted.bam


$(ANALYSIS_DIR)/%_dedup.bam: $(ANALYSIS_DIR)/%.bam
	picard MarkDuplicates INPUT=$< OUTPUT=$@ METRICS_FILE=$(basename $@)_metrics.txt
	picard BuildBamIndex INPUT=$@

$(ANALYSIS_DIR)/on_sim/%_raw_variants_uncalibrated.vcf: $(ANALYSIS_DIR)/on_sim/%_dedup.bam
	gatk -T HaplotypeCaller \
		-R $(SIMFASTA2) \
		-I $< \
		--genotyping_mode DISCOVERY \
		-fixMisencodedQuals \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sec/%_raw_variants_uncalibrated.vcf: $(ANALYSIS_DIR)/on_sec/%_dedup.bam
	gatk -T HaplotypeCaller \
		-R $(SECFASTA2) \
		-I $< \
		--genotyping_mode DISCOVERY \
		-fixMisencodedQuals \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sim/%_on_sim.bam: $(%) $(REFDIR)/Dsim_unspliced/Genome | $(ANALYSIS_DIR)/on_sim
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(REFDIR)/Dsim_unspliced \
		--outFileNamePrefix $(dir $@)$* \
		--outTmpDir $(basename $@)_tmp/ \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) 
	samtools view -bS -o $@ $(@D)/$*Aligned.out.sam
	rm $(@D)/$*Aligned.out.sam

$(ANALYSIS_DIR)/on_sec/%_on_sec.bam: $(%) $(REFDIR)/Dsec_unspliced/Genome | $(ANALYSIS_DIR)/on_sec
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(REFDIR)/Dsec_unspliced \
		--outFileNamePrefix $(dir $@)$* \
		--outTmpDir $(basename $@)_tmp \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) 
	samtools view -bS -o $@ $(@D)/$*Aligned.out.sam
	rm $(@D)/$*Aligned.out.sam

$(ANALYSIS_DIR)/on_sim: | $(ANALYSIS_DIR)
	mkdir $@

$(ANALYSIS_DIR)/on_sec: | $(ANALYSIS_DIR)
	mkdir $@
