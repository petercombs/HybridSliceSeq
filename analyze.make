include gmsl

comma := ,
space := $(null) #

wittkopp_reads:  sequence/SRR869573 sequence/SRR869575 sequence/SRR869590 sequence/SRR869596 sequence/SRR869597 sequence/SRR869598 sequence/SRR869599 sequence/SRR869600 sequence/SRR869601 sequence/SRR869602 sequence/SRR869844 sequence/SRR869603 sequence/SRR869905 sequence/SRR869592 sequence/SRR869579 sequence/SRR869580 sequence/SRR869587
	echo "All Downloaded"

sim_gdna = sequence/SRR869579 sequence/SRR869580
sec_gdna = sequence/SRR869587
simXsec = sequence/SRR869602 sequence/SRR869844
secXsim = sequence/SRR869603 sequence/SRR869905
sim_rna = sequence/SRR869573 sequence/SRR869575
sec_rna = sequence/SRR869601

mel_gdna = sequence/SRR492062 sequence/SRR492061
simXmel = sequence/SRR869592
melXsim = sequence/SRR869590
simPLUSmel = sequence/SRR869596	sequence/SRR869597 sequence/SRR869598 sequence/SRR869599 sequence/SRR869600

QSUBBER_ARGS = --keep-temporary tmp --log-base $(basename $@)  --job-name $(basename $(@F))

SIMSEC_REFSIMFASTA2 = analysis/on_sim/simsec_variants.fasta
SIMSEC_REFSECFASTA2 = analysis/on_sec/simsec_variants.fasta

sequence/%:
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_1.fastq.gz
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_2.fastq.gz
	touch $@


$(ANALYSIS_DIR)/on_sim/%_bowtie.bam: $(%) $(basename $(SIMFASTA2)).1.ebwt | $(ANALYSIS_DIR)
	echo $<
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

$(ANALYSIS_DIR)/on_sim/%_bowtie2.bam: $(%) $(basename $(SIMFASTA2)).1.bt2 | $(ANALYSIS_DIR)
	./qsubber --job-name $*_bowtie2 --queue batch --keep-temporary tmp -t 8 \
		-l mem=2gb -l pmem=2gb --log-base $(basename $@) \
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
		\| samtools view -bS -o $(basename $@)_unsorted.bam -
	samtools sort $(basename $@)_unsorted.bam $(basename $@)
	rm $(basename $@)_unsorted.bam

$(ANALYSIS_DIR)/on_sec/%_bowtie2.bam: $(%) $(basename $(SECFASTA2)).1.bt2 | $(ANALYSIS_DIR)
	./qsubber --job-name $*_bowtie2 --keep-temporary tmp \
		--queue batch -l mem=2gb -l pmem=2gb -t 8 --log-base $(basename $@) \
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
	picard MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=$< OUTPUT=$@ METRICS_FILE=$(basename $@)_metrics.txt
	picard BuildBamIndex INPUT=$@

$(ANALYSIS_DIR)/on_sim/%_raw_variants_uncalibrated.vcf: $(ANALYSIS_DIR)/on_sim/%_dedup.bam $(SIMFASTA2).fai $(basename $(SIMFASTA2)).dict
	./qsubber $(QSUBBER_ARGS) -t 1  \
	gatk -T HaplotypeCaller \
		-R $(SIMFASTA2) \
		-I $< \
		--genotyping_mode DISCOVERY \
		-fixMisencodedQuals \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sim/%_raw_variants_uncalibrated.g.vcf: $(ANALYSIS_DIR)/on_sim/%_dedup.bam $(SIMFASTA2).fai $(basename $(SIMFASTA2)).dict
	./qsubber $(QSUBBER_ARGS) -t 1  \
	gatk -T HaplotypeCaller \
		-R $(SIMFASTA2) \
		-I $< \
		--genotyping_mode DISCOVERY \
		-fixMisencodedQuals \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sim/%_raw_variants_uncalibrated.p.g.vcf: $(ANALYSIS_DIR)/on_sim/%_dedup.bam $(SIMFASTA2).fai $(basename $(SIMFASTA2)).dict
	./qsubber $(QSUBBER_ARGS) -t 8  \
	gatk -T HaplotypeCaller \
		-R $(SIMFASTA2) \
		-I $< \
		-nct 8  \
		--genotyping_mode DISCOVERY \
		-fixMisencodedQuals \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sec/%_raw_variants_uncalibrated.vcf: $(ANALYSIS_DIR)/on_sec/%_dedup.bam $(SECFASTA2).fai $(basename $(SECFASTA2)).dict
	./qsubber $(QSUBBER_ARGS) -t 1  \
	gatk -T HaplotypeCaller \
		-R $(SECFASTA2) \
		-I $< \
		--genotyping_mode DISCOVERY \
		-fixMisencodedQuals \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sec/%_raw_variants_uncalibrated.g.vcf: $(ANALYSIS_DIR)/on_sec/%_dedup.bam $(SECFASTA2).fai $(basename $(SECFASTA2)).dict
	./qsubber $(QSUBBER_ARGS) -t 1  \
	gatk -T HaplotypeCaller \
		-R $(SECFASTA2) \
		-I $< \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-fixMisencodedQuals \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_sec/%_raw_variants_uncalibrated.p.g.vcf: $(ANALYSIS_DIR)/on_sec/%_dedup.bam $(SECFASTA2).fai $(basename $(SECFASTA2)).dict
	./qsubber $(QSUBBER_ARGS) -t 8  \
	gatk -T HaplotypeCaller \
		-R $(SECFASTA2) \
		-I $< \
		-nct 8 \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-fixMisencodedQuals \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@

$(ANALYSIS_DIR)/on_%_STAR_RNASEQ.bam: 
	./qsubber  $(QSUBBER_ARGS) --resource "mem=2gb" -t 8 \
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(dir $@)masked \
		--outFileNamePrefix $(@D)/$(notdir $*) \
		--outFilterMultimapNmax 1 \
		--outSAMattributes MD NH --clip5pNbases 6 \
		--outTmpDir $(basename $@)_tmp/ \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($(notdir $*))))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($(notdir $*))))) 
	./qsubber $(QSUBBER_ARGS) -t 4 \
	samtools view -bS $(@D)/$(notdir $*)Aligned.out.sam \| samtools sort -@ 4 - $(basename $@)
	rm $(@D)/$(notdir $*)Aligned.out.sam
	@echo pass


$(ANALYSIS_DIR)/on_sim/%_STAR.bam: $(%) $(REFDIR)/Dsim_unspliced/Genome | $(ANALYSIS_DIR)/on_sim
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(REFDIR)/Dsim_unspliced \
		--outFileNamePrefix $(dir $@)$* \
		--outTmpDir $(basename $@)_tmp/ \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*)))) 
	samtools view -bS -o $@ $(@D)/$*Aligned.out.sam
	rm $(@D)/$*Aligned.out.sam

$(ANALYSIS_DIR)/on_sec/%_STAR.bam: $(%) $(REFDIR)/Dsec_unspliced/Genome | $(ANALYSIS_DIR)/on_sec
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


$(ANALYSIS_DIR)/on_%/simsec_variants.tsv: $($(call uc,%)FASTA2)   $(ANALYSIS_DIR)/on_%/sim_gdna_bowtie2_raw_variants_uncalibrated.p.g.vcf $(ANALYSIS_DIR)/on_%/sec_gdna_bowtie2_raw_variants_uncalibrated.p.g.vcf
	echo $^
	gatk -T GenotypeGVCFs \
		-R $($(call uc,$*)FASTA2) \
		-V $(@D)/sim_gdna_bowtie2_raw_variants_uncalibrated.p.g.vcf \
		-V $(@D)/sec_gdna_bowtie2_raw_variants_uncalibrated.p.g.vcf \
		-o $(@D)/simsec_variants_on_$*.gvcf
	gatk -T VariantsToTable \
		-R $($(call uc,$*)FASTA2) \
		-V $(@D)/simsec_variants_on_$*.gvcf \
		-F CHROM -F POS -F REF -F ALT -F QUAL \
		-F HET -F HOM-REF -F HOM-VAR -F NCALLED \
		-o $@

$(ANALYSIS_DIR)/on_%/simsec_masked.fasta $(ANALYSIS_DIR)/on_%/simsec_variant.bed: $(ANALYSIS_DIR)/on_%/simsec_variants.tsv
	python MaskReferenceFromGATKTable.py \
		--threads 4 \
		--emit-bed $(@D)/simsec_variant.bed \
		$($(call uc,$*)FASTA2) \
		$< \
		$(@D)/simsec_masked.fasta

$(ANALYSIS_DIR)/on_%/masked/Genome: $(ANALYSIS_DIR)/on_%/simsec_masked.fasta | $(ANALYSIS_DIR)/on_%/masked
	STAR --runMode genomeGenerate --genomeDir $(dir $@) \
		--outTmpDir $(@D)/_tmp/ \
		--genomeFastaFiles $< \
		--sjdbGTFfile $($(call uc,$*)GTF)

%_wasp_dup_removed.bam : %_STAR_RNASEQ.bam
	./qsubber $(QSUBBER_ARGS) -t 4 \
	python2 rmdup_pe.py $< - \
		\| samtools sort -n -@ 3 - $(basename $@)

%_SNP_COUNTS.txt : %_wasp_dup_removed.bam 
	mkdir -p $*_countsnpase_tmp
	./qsubber $(QSUBBER_ARGS) \
	python2 CountSNPASE.py \
		--mode single \
		--reads $< \
		--snps $(@D)/simsec_variant.bed \
		--prefix $*_countsnpase_tmp/
	mv $*_countsnpase_tmp/_SNP_COUNTS.txt $@
	
%_gene_ase.tsv : %_SNP_COUNTS.txt
	python2 GetGeneASE.py \
		--snpcounts $< \
		--phasedsnps $(@D)/simsec_variant.bed \
		--gff $($(call uc,$(call substr,$(notdir $(@D)),4,6))GTF) \
		-o $@ \
		--writephasedsnps

$(ANALYSIS_DIR)/on_%/masked:
	mkdir $@
