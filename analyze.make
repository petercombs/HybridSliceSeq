include gmsl

comma := ,
space := $(null) #

wittkopp_reads:  sequence/SRR869573 sequence/SRR869575 sequence/SRR869590 sequence/SRR869596 sequence/SRR869597 sequence/SRR869598 sequence/SRR869599 sequence/SRR869600 sequence/SRR869601 sequence/SRR869602 sequence/SRR869844 sequence/SRR869603 sequence/SRR869905 sequence/SRR869592
	echo "All Downloaded"

wittkopp_sim_gdna = sequence/SRR869579 sequence/SRR869580
wittkopp_sec_gdna = sequence/SRR869587
simXsec = sequence/SRR869602 sequence/SRR869844
secXsim = sequence/SRR869603 sequence/SRR869905

sequence/%:
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_1.fastq.gz
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_2.fastq.gz
	touch $@

$(ANALYSIS_DIR)/on_sim/%_on_sim.bam: $(%) $(REFDIR)/Dsim_unspliced/Genome | $(ANALYSIS_DIR)/on_sim
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(REFDIR)/Dsim_unspliced \
		--outFileNamePrefix $(dir $@) \
		--outTmpDir $(basename $@)_tmp \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) 
	samtools view -bS -o $(@D)/accepted_hits.bam $(@D)/Aligned.out.sam
	rm $(@D)/Aligned.out.sam

$(ANALYSIS_DIR)/on_sec/%_on_sec.bam: $(%) $(REFDIR)/Dsec_unspliced/Genome | $(ANALYSIS_DIR)/on_sec
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(REFDIR)/Dsec_unspliced \
		--outFileNamePrefix $(dir $@) \
		--outTmpDir $(basename $@)_tmp \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) 
	samtools view -bS -o $(@D)/accepted_hits.bam $(@D)/Aligned.out.sam
	rm $(@D)/Aligned.out.sam

$(ANALYSIS_DIR)/on_sim: | $(ANALYSIS_DIR)
	mkdir $@

$(ANALYSIS_DIR)/on_sec: | $(ANALYSIS_DIR)
	mkdir $@
