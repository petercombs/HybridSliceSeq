include gmsl

VARIANTS=analysis_godot/on_melr5/melsim_variant.bed

.SECONDEXPANSION:

comma := ,
space := $(null) #


sim_gdna = sequence/SRR520334
mel_gdna = sequence/SRR835939

mel_dnase = sequence/SRR060797_se sequence/SRR060797_se sequence/SRR060798_se sequence/SRR060799_se

oemale = sequence/150812_BRISCOE_0250_AC7K8CACXX_L2_OEmale
oefemale = sequence/150812_BRISCOE_0250_AC7K8CACXX_L2_OEfemale
fbfemale = sequence/150812_BRISCOE_0250_AC7K8CACXX_L2_FBfemale
fbmale = sequence/150812_BRISCOE_0250_AC7K8CACXX_L2_FBmale

oemale_r2   = sequence/150814_BRISCOE_0251_BC7LJDACXX_L1_ACTTGATG_pf
oefemale_r2 = sequence/150814_BRISCOE_0251_BC7LJDACXX_L1_TGACAGAC_pf
fbfemale_r2 = sequence/150814_BRISCOE_0251_BC7LJDACXX_L1_CGATGTTT_pf
fbmale_r2   = sequence/150814_BRISCOE_0251_BC7LJDACXX_L1_TAGAACAC_pf

QSUBBER_ARGS =  $(QSUBBER_LOCAL) --sleep-random 10 --keep-temporary tmp --log-base $(basename $@)  --job-name $(basename $(@F))

sequence/%:
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_1.fastq.gz
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_2.fastq.gz
	touch $@

sequence/%_se:
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*.fastq.gz
	touch $@

$(ANALYSIS_DIR)/on_sim/%_bowtie.bam: $(%) $(basename $(SIMFASTA2)).1.ebwt | $(ANALYSIS_DIR)
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

$(ANALYSIS_DIR)/on_%_bowtie2_se.bam: $$(basename $$($$(call uc,$$(subst /,,$$(dir $$*)))FASTA2)).1.bt2 $$($$(notdir $$*)) | $(ANALYSIS_DIR)
	mkdir -p $(@D)
	./qsubber --job-name $*_bowtie2 --queue batch --keep-temporary tmp -t 4 \
		-l mem=8gb -l pmem=8gb --log-base $(basename $@) \
	bowtie2 \
		--very-sensitive \
		-p 8 \
		--rg-id $(notdir $*) \
		--rg "SM:$(notdir $*)" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x $(basename $(basename $<)) \
		-U $(subst $(space),$(comma),$(strip $(patsubst %_se, %.fastq.gz, $($(notdir $*))))) \
		\| samtools view -bS - \
		\| samtools sort - $(basename $@)

$(ANALYSIS_DIR)/on_%_bowtie2.bam: $$(basename $$($$(call uc,$$(subst /,,$$(dir $$*)))FASTA2)).1.bt2 $$($$(notdir $$*)) | $(ANALYSIS_DIR)
	mkdir -p $(@D)
	./qsubber --job-name $*_bowtie2 --queue batch --keep-temporary tmp -t 4 \
		-l mem=8gb -l pmem=8gb --log-base $(basename $@) \
	bowtie2 \
		--very-sensitive-local \
		-p 8 \
		--rg-id $(notdir $*) \
		--rg "SM:$(notdir $*)" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x $(basename $(basename $<)) \
		-1 $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($(notdir $*))))) \
		-2 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($(notdir $*))))) \
		\| samtools view -bS - \
		\| samtools sort - $(basename $@)




$(ANALYSIS_DIR)/%_dedup.bam: $(ANALYSIS_DIR)/%.bam
	picard MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		READ_NAME_REGEX=null \
		INPUT=$< OUTPUT=$@ METRICS_FILE=$(basename $@)_metrics.txt
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

$(ANALYSIS_DIR)/on_%_raw_variants_uncalibrated.p.g.vcf: $$($$(subst /,,$$(call uc,$$(dir $$*)))FASTA2).fai $$(basename $$($$(subst /,,$$(call uc,$$(dir $$*)))FASTA2)).dict $(ANALYSIS_DIR)/on_%_bowtie2_dedup.bam
	./qsubber $(QSUBBER_ARGS) -t 4  \
	gatk -T HaplotypeCaller \
		-R $(basename $<) \
		-I $(ANALYSIS_DIR)/on_$*_bowtie2_dedup.bam \
		-nct 8 \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o $@


$(ANALYSIS_DIR)/on_%_STAR_RNASEQ.bam: $$(@D)/masked/Genome $$($$(notdir $$*))
	rm -rf $(basename $@)_tmp/
	./qsubber  $(QSUBBER_ARGS) --resource "mem=8gb" -t 4 \
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


$(ANALYSIS_DIR)/on_sim/%_STAR.bam: $(%) $(REFDIR)/Dsim_unspliced/Genome | $(ANALYSIS_DIR)/on_sim
	module load STAR && \
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(REFDIR)/Dsim_unspliced \
		--outFileNamePrefix $(dir $@)$* \
		--outTmpDir $(basename $@)_tmp/ \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($*)))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($*))))
	samtools view -bS -o $@ $(@D)/$*Aligned.out.sam
	rm $(@D)/$*Aligned.out.sam

$(ANALYSIS_DIR)/on_%_STAR.bam: $(REFDIR)/D$$(subst /,,$$(dir $$*))_unspliced/Genome
	mkdir -p $(@D)
	rm -rf $(basename $@)_tmp
	./qsubber $(QSUBBER_ARGS) -t 4 \
	STAR \
		--parametersFiles $(STARCONFIG) \
		--genomeDir $(patsubst %/,%,$(dir $<)) \
		--outFileNamePrefix $(dir $@)$(notdir $*) \
		--outTmpDir $(basename $@)_tmp \
		--readFilesIn $(subst $(space),$(comma),$(strip $(patsubst %, %_1.fastq.gz, $($(notdir $*))))) \
				 $(subst $(space),$(comma),$(strip $(patsubst %, %_2.fastq.gz, $($(notdir $*)))))
	samtools view -bS -o $@ $(@D)/$(notdir $*)Aligned.out.sam
	rm $(@D)/$(notdir $*)Aligned.out.sam

$(ANALYSIS_DIR)/on_sim: | $(ANALYSIS_DIR)
	mkdir $@

$(ANALYSIS_DIR)/on_%/ : | $(ANALYSIS_DIR)
	mkdir $@


$(ANALYSIS_DIR)/on_%_variants.tsv: $$($$(call uc,$$(subst /,,$$(dir $$*)))FASTA2)   \
	$(ANALYSIS_DIR)/on_$$(dir $$*)$$(call substr,$$(notdir $$*),1,3)_gdna_raw_variants_uncalibrated.p.g.vcf \
	$(ANALYSIS_DIR)/on_$$(dir $$*)$$(call substr,$$(notdir $$*),4,6)_gdna_raw_variants_uncalibrated.p.g.vcf
	gatk -T GenotypeGVCFs \
		-R $($(call uc,$(patsubst %/,%,$(dir $*)))FASTA2) \
		-V $(@D)/$(call substr,$(notdir $*),1,3)_gdna_raw_variants_uncalibrated.p.g.vcf \
		-V $(@D)/$(call substr,$(notdir $*),4,6)_gdna_raw_variants_uncalibrated.p.g.vcf \
		-o $(@D)/$(notdir $*)_variants_on_$(patsubst %/,%,$(dir $*)).gvcf
	gatk -T VariantsToTable \
		-R $($(call uc,$(patsubst %/,%,$(dir $*)))FASTA2) \
		-V $(@D)/$(notdir $*)_variants_on_$(patsubst %/,%,$(dir $*)).gvcf \
		-F CHROM -F POS -F REF -F ALT -F QUAL \
		-F HET -F HOM-REF -F HOM-VAR -F NCALLED \
		-GF GT \
		-o $@

$(ANALYSIS_DIR)/on_%_masked.fasta $(ANALYSIS_DIR)/on_%_variant.bed: $(ANALYSIS_DIR)/on_%_variants.tsv
	python MaskReferenceFromGATKTable.py \
		--target-species $(patsubst %/,%,$(dir $*)) \
		--emit-bed $(@D)/$(notdir $*)_variant.bed \
		--outfasta $(@D)/$(notdir $*)_masked.fasta \
		$($(call uc,$(patsubst %/,%,$(dir $*)))FASTA2) \
		$<

$(ANALYSIS_DIR)/on_mel/masked/Genome: $(ANALYSIS_DIR)/on_mel/melsim_masked.fasta | $(ANALYSIS_DIR)/on_%/masked
	rm -rf $(@D)/_tmp
	module load STAR && \
	STAR --runMode genomeGenerate --genomeDir $(dir $@) \
		--outTmpDir $(@D)/_tmp/ \
		--genomeFastaFiles $< \
		--sjdbGTFfile $(MELGTF)

$(ANALYSIS_DIR)/on_%/masked/Genome: $(ANALYSIS_DIR)/on_%/melsim_masked.fasta | $(ANALYSIS_DIR)/on_%/masked
	rm -rf $(@D)/_tmp
	module load STAR && \
	STAR --runMode genomeGenerate --genomeDir $(dir $@) \
		--outTmpDir $(@D)/_tmp/ \
		--genomeFastaFiles $< \
		--sjdbGTFfile $($(call uc,$*)GTF)

Reference/melsim/Genome: $(ANALYSIS_DIR)/on_melr5/melsim_masked.fasta
	rm -rf $(@D)/_tmp
	module load STAR && \
	STAR --runMode genomeGenerate --genomeDir $(dir $@) \
		--outTmpDir $(@D)/_tmp/ \
		--genomeFastaFiles $< \
		--sjdbGTFfile $(MELGTFR5)

$(ANALYSIS_DIR)/on_%/masked/transcriptome: $(ANALYSIS_DIR)/on_%/melsim_masked.1.bt2 | $(ANALYSIS_DIR)/on_%/masked
	mkdir -p $(@D)
	tophat2 --transcriptome-index $@ \
		--GTF $($(call uc,$*)GTF) \
		$(basename $(basename $<))
	touch $@

$(ANALYSIS_DIR)/on_mel/masked/transcriptome: $(ANALYSIS_DIR)/on_mel/melsim_masked.1.bt2 | $(ANALYSIS_DIR)/on_%/masked
	tophat2 --transcriptome-index $@ \
		--GTF $(MELGTF) \
		$(basename $(basename $<))
	touch $@

$(ANALYSIS_DIR)/on_%/masked/index_kal: $(ANALYSIS_DIR)/on_%/masked/transcriptome
	kallisto index \
		-i $@ \
		$<.fa

$(ANALYSIS_DIR)/on_%/abundance.tsv: $(ANALYSIS_DIR)/on_$$(firstword $$(call split,/,$$*))/masked/index_kal
	mkdir -p $(@D)
	./qsubber $(QSUBBER_ARGS) -t 1 \
	kallisto quant \
		--index $<\
		--plaintext \
		--output-dir $(@D) \
		$(foreach F,$($(notdir $*)), $F_1.fastq.gz $F_2.fastq.gz)


$(ANALYSIS_DIR)/on_%/abundance.tsv: $(ANALYSIS_DIR)/on_$$(firstword $$(call split,/,$$*))/masked/index_kal
	mkdir -p $(@D)
	./qsubber $(QSUBBER_ARGS) -t 1 \
	kallisto quant \
		--index $<\
		--plaintext \
		--output-dir $(@D) \
		$(foreach F,$($(notdir $*)), $F_1.fastq.gz $F_2.fastq.gz)

%_wasp_dedup.bam: %.bam $(ANALYSIS_DIR)/deduplicate
	./qsubber $(QSUBBER_ARGS) -t 6 --resource mem=32000m --load-module picard/2.8.1 \
	picard MarkDuplicates \
		I=$< O=$@ \
		METRICS_FILE=$@.metrics \
		REMOVE_DUPLICATES=true \
		DUPLICATE_SCORING_STRATEGY=RANDOM

%/melsim_SNP_COUNTS.txt : %/assigned_dmelR_wasp_dedup.bam $(VARIANTS)
	mkdir -p $*/melsim_countsnpase_tmp
	./qsubber --resource mem=4000m $(QSUBBER_ARGS) \
	python2 CountSNPASE.py \
		--mode single \
		--reads $< \
		--snps $(VARIANTS) \
		$(if $(wildcard $*/is_single), -t) \
		--prefix $*/melsim_countsnpase_tmp/
	mv $*/melsim_countsnpase_tmp/_SNP_COUNTS.txt $@

%/melsim_wasp_SNP_COUNTS.txt : %/accepted_hits.remap.keep.sorted_wasp_dedup.bam
	mkdir -p $*/melsim_countsnpase_tmp
	./qsubber $(QSUBBER_ARGS) \
	python2 CountSNPASE.py \
		--mode single \
		--reads $< \
		--snps $(VARIANTS) \
		$(if $(wildcard $*/is_single), -t) \
		--prefix $*/melsim_countsnpase_tmp/
	mv $*/melsim_countsnpase_tmp/_SNP_COUNTS.txt $@


%/simmel_SNP_COUNTS.txt : %/on_simaccepted_hits_remapped_sorted_wasp_dedup.bam $(ANALYSIS_DIR)/on_sim/melsim_variant.bed
	mkdir -p $*/simmel_countsnpase_tmp
	./qsubber $(QSUBBER_ARGS) \
	python2 CountSNPASE.py \
		--mode single \
		--reads $< \
		--snps $(ANALYSIS_DIR)/on_sim/melsim_variant.bed \
		$(if $(wildcard $*/is_single), -t) \
		--prefix $*/simmel_countsnpase_tmp/
	mv $*/simmel_countsnpase_tmp/_SNP_COUNTS.txt $@

%_gene_ase.tsv : %_SNP_COUNTS.txt GetGeneASE.py $(ANALYSIS_DIR)/recalc_ase $(ANALYSIS_DIR)/on_$$(call substr,$$(notdir $$@),1,3)/true_hets.tsv
	./qsubber $(QSUBBER_ARGS) -t 1 \
	python2 GetGeneASE.py \
		--snpcounts $< \
		--phasedsnps $(ANALYSIS_DIR)/on_$(call substr,$(notdir $@),1,3)/melsim_variant.bed \
		--allele-min 0 \
		--min 5 \
		--gff $($(call uc,$(call substr,$(notdir $@),1,3))GTF) \
		-o $@ \
		--preference-index \
		--true-hets $(ANALYSIS_DIR)/on_$(call substr,$(notdir $@),1,3)/true_hets.tsv \
		--writephasedsnps

%/simmel_gene_ase_by_read.tsv:  %/on_simaccepted_hits_remapped_sorted_wasp_dedup.sorted.bam $(ANALYSIS_DIR)/on_sim/melsim_variant.bed $(ANALYSIS_DIR)/recalc_ase
	./qsubber $(QSUBBER_ARGS) -t 1 \
	python ~/ASEr/bin/GetGeneASEbyReads.py \
		--outfile $@ \
		--id-name gene_name \
		--ase-function pref_index \
		$(ANALYSIS_DIR)/on_sim/melsim_variant.bed \
		$(SIMGTF) \
		$<

%/melsim_gene_ase_by_read.tsv : %/assigned_dmelR_wasp_dedup.sorted.bam %/assigned_dmelR_wasp_dedup.sorted.bam.bai $(VARIANTS) $(ANALYSIS_DIR)/recalc_ase
	./qsubber $(QSUBBER_ARGS) -t 1 \
	python ~/ASEr/bin/GetGeneASEbyReads.py \
		--outfile $@ \
		--id-name gene_name \
		--ase-function pref_index \
		$(VARIANTS) \
		$(MELGTFR5) \
		$<

%/melsim_gene_ase_by_read_with_wasp.tsv : %/assigned_dmelR_wasp_dedup.remap.kept_sorted.bam %/assigned_dmelR_wasp_dedup.remap.kept.sorted.bam.bai $(VARIANTS) $(ANALYSIS_DIR)/recalc_ase
	./qsubber $(QSUBBER_ARGS) -t 1 \
	python ~/ASEr/bin/GetGeneASEbyReads.py \
		--outfile $@ \
		--ase-function pref_index \
		$(VARIANTS) \
		$(MELGTFR5) \
		$<

%/wasp_genes.fpkm_tracking: %/assigned_dmelR_wasp_dedup.remap.kept_sorted.bam $(MELGTFR5) $(MELFASTA2) $(MELBADGTF)
	./qsubber $(QSUBBER_ARGS)_$(*F) --load-module cufflinks -t 4 \
	cufflinks \
		--num-threads 8 \
		--output-dir $(@D)/tmp \
		--multi-read-correct \
		--frag-bias-correct $(MELR5FASTA2) \
		--GTF $(MELGTFR5) \
		--mask-file $(MELBADGTFR5) \
		$<
	mv $(@D)/tmp/genes.fpkm_tracking $@


%_cds_ase.tsv : %_SNP_COUNTS.txt GetGeneASE.py $(ANALYSIS_DIR)/recalc_ase $(ANALYSIS_DIR)/on_mel/true_hets.tsv
	./qsubber $(QSUBBER_ARGS) -t 1 \
	python2 GetGeneASE.py \
		--snpcounts $< \
		--phasedsnps $(VARIANTS) \
		--allele-min 0 \
		--min 5 \
		--type CDS \
		--gff $($(call uc,$(call substr,$(notdir $@),1,3))GTF) \
		-o $@ \
		--preference-index \
		--true-hets $(ANALYSIS_DIR)/on_mel/true_hets.tsv \
		--writephasedsnps

$(ANALYSIS_DIR)/recalc_ase:
	touch $@

%_gene_ase.tsv : $$(firstword $$(subst _counts_, ,$$@))_SNP_COUNTS.txt GetGeneASE.py analyze.make
	echo $(firstword $(subst _counts_, ,$*))_SNP_COUNTS.txt
	echo $(lastword $(subst _counts_, ,$*))
	python2 GetGeneASE.py \
		--snpcounts $< \
		--phasedsnps $(VARIANTS) \
		--allele-min 0 \
		--min $(lastword $(subst _counts_, ,$*)) \
		--gff $($(call uc,$(call substr,$(notdir $@),1,3))GTF) \
		-o $@ \
		--writephasedsnps
		#--true-hets analysis/on_sim/true_hets.tsv \

$(ANALYSIS_DIR)/on_%/masked:
	mkdir -p $@

%/genes.fpkm_tracking: %_STAR_RNASEQ.bam
	mkdir -p $(@D)
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 4 \
	cufflinks \
		--GTF $($(call uc,$(call substr,$(notdir $(*D)),4,6))GTF) \
		--num-threads 8 \
		--frag-bias-correct $(@D)/../melsim_masked.fasta \
		--multi-read-correct \
		--output-dir $(@D) \
		--quiet \
		$<

%accepted_hits_remapped.bam: #$$(dir $$*)/accepted_hits.bam #$(ANALYSISDIR)/$$(notdir $$*)/masked/Genome
	./qsubber $(QSUBBER_ARGS)_$(*F) --load-module STAR -t 4 \
	python RemapOnOtherGenome.py $(@D)/accepted_hits.bam $(ANALYSIS_DIR)/$(notdir $*)/masked/Genome $@

WASPMAP = $(HOME)/FWASP/mapping

%.remap.fq1.gz %.remap.fq2.gz %.keep.bam %.to.remap.bam : %.bam $(WASPMAP)/find_intersecting_snps_2.py
	./qsubber $(QSUBBER_ARGS) -t 1 python $(WASPMAP)/find_intersecting_snps_2.py --phased -p $< analysis_godot/on_melr5 \
		\> $*_find_snps.log

%.remap.bam: %.remap.fq1.gz %.remap.fq2.gz
	./qsubber $(QSUBBER_ARGS) --resource "mem=8gb" -t 4 --load-module STAR \
		STAR --parametersFiles $(STARCONFIG) \
		--genomeDir Reference/melsim --outFileNamePrefix $(@D)/ \
		--outSAMattributes MD NH --clip5pNbases 6 \
		--readFilesIn $^
	./qsubber $(QSUBBER_ARGS) --resource "mem=6gb" -t 4 --load-module samtools \
	samtools view -bS -o $@ $(@D)/Aligned.out.sam
	rm $(@D)/Aligned.out.sam

%.remap.kept.bam: %.remap.bam %.to.remap.num.gz
	./qsubber $(QSUBBER_ARGS) --resource "mem=8gb" -t 1 \
	python2 $(WASPMAP)/filter_remapped_reads.py \
		-p \
		$*.to.remap.bam \
		$*.remap.bam \
		$*.remap.kept.bam \
		$*.to.remap.num.gz

%.to.remap.num.gz: %.remap.fq1.gz
	python $(WASPMAP)/make_num_from_fq.py $<

%.keep.merged.bam: %.keep.bam %.remap.keep.bam
	samtools merge -f $@ $^

%.sorted.bam: %.bam
	./qsubber $(QSUBBER_ARGS) --resource "mem=8gb" -t 4 \
	samtools sort $(SORT_OPTS) -@ 4 $< $*.sorted
	samtools index $@

%/mel_only.bam %/sim_only.bam: %/assigned_dmelR_wasp_dedup.sorted.bam
	./qsubber $(QSUBBER_ARGS) --resource "mem=8gb" -t 4 \
	 python PartitionReads.py \
		 --reference-species mel --alt-species sim \
		 --paired \
		 $(VARIANTS) \
		 $<
	mv $(@D)/assigned_dmelR_wasp_dedup.sorted_sim.bam $(@D)/sim_only.bam
	mv $(@D)/assigned_dmelR_wasp_dedup.sorted_mel.bam $(@D)/mel_only.bam
	mv $(@D)/assigned_dmelR_wasp_dedup.sorted_ambig.bam $(@D)/ambig_only.bam


%_only/genes.fpkm_tracking: %_only_sorted.bam
	mkdir -p $(@D)
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 4 --load-module cufflinks \
	cufflinks \
		--GTF-guide $(MELGTF) \
		--num-threads 8 \
		--frag-bias-correct $(ANALYSIS_DIR)/on_mel/melsim_masked.fasta \
		--multi-read-correct \
		--output-dir $(@D) \
		--quiet \
		$<


%/sign_test.tsv: %/melsim_gene_ase_by_read.tsv $(ANALYSIS_DIR)/redo_sign_test
	./qsubber $(QSUBBER_ARGS)_$(*F) \
	python sign_test.py \
		--drop-genes analysis/results/strict_maternal_gene_names.txt \
		--print-header \
		--print-counts \
		--min-genes 10 \
		$< \
		Reference/mel_goterms.gmt \
		$@ \
		200

$(ANALYSIS_DIR)/redo_sign_test:
	touch $@



$(ANALYSIS_DIR)/sign_test_summary.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,sign_test.tsv,$$(FPKMS)) $(ANALYSIS_DIR)/map_stats.tsv
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename sign_test.tsv \
			--column "pval" \
			--strip-low-reads 3000000 \
			--strip-on-unique \
			--strip-as-nan \
			--map-stats $(ANALYSIS_DIR)/map_stats.tsv \
			--key "category" \
			--out-basename sign_test_summary \
			--float-format "%0.18g" \
			$(ANALYSIS_DIR)


$(ANALYSIS_DIR)/sign_test_or_summary.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,sign_test.tsv,$$(FPKMS)) $(ANALYSIS_DIR)/map_stats.tsv
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename sign_test.tsv \
			--column "oddsratio" \
			--strip-low-reads 3000000 \
			--strip-on-unique \
			--strip-as-nan \
			--map-stats $(ANALYSIS_DIR)/map_stats.tsv \
			--key "category" \
			--out-basename sign_test_or_summary \
			--float-format "%0.18g" \
			$(ANALYSIS_DIR)
