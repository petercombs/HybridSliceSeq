# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis_godot

# Reference FASTA and GFF files from FlyBase
MEL5RELEASE= r5.57_FB2014_03
MELRELEASE = r6.09_FB2016_01
SIMRELEASE = r2.01_FB2016_01
SECRELEASE = r1.3_FB2016_01

MELMAJORVERSION = $(word 1, $(subst ., , $(MELRELEASE)))
MELVERSION = $(word 1, $(subst _FB, ,$(MELRELEASE)))
MELDATE = $(word 2, $(subst _FB, ,$(MELRELEASE)))

SIMMAJORVERSION = $(word 1, $(subst ., , $(SIMRELEASE)))
SIMVERSION = $(word 1, $(subst _FB, ,$(SIMRELEASE)))
SIMDATE = $(word 2, $(subst _FB, ,$(SIMRELEASE)))


SECMAJORVERSION = $(word 1, $(subst ., , $(SECRELEASE)))
SECVERSION = $(word 1, $(subst _FB, ,$(SECRELEASE)))
SECDATE = $(word 2, $(subst _FB, ,$(SECRELEASE)))

PREREQDIR = prereqs
MELFASTA = $(PREREQDIR)/dmel-all-chromosome-$(MELVERSION).fasta
MELR5FASTA = $(PREREQDIR)/dmel-all-chromosome-r5.57.fasta
SIMFASTA = $(PREREQDIR)/dsim-all-chromosome-$(SIMVERSION).fasta

REFDIR = Reference

MELFASTA2= $(REFDIR)/dmel_prepend.fasta
MELR5FASTA2= $(REFDIR)/dmelr5_prepend.fasta
SIMFASTA2= $(REFDIR)/dsim_prepend.fasta

ORTHOLOGS = $(PREREQDIR)/gene_orthologs_fb_$(MELDATE).tsv

MELGFF   = $(PREREQDIR)/dmel-all-$(MELVERSION).gff
MELGFFR5 = $(PREREQDIR)/dmel-all-r5.57.gff
MELGTF   = $(REFDIR)/mel_good.gtf
MELGTFR5 = $(REFDIR)/mel_r5_good.gtf
MELALLGTF   = $(REFDIR)/mel_all.gtf
MELALLGTFR5 = $(REFDIR)/mel_r5_all.gtf
MELBADGTF   = $(REFDIR)/mel_bad.gtf
MELBADGTFR5 = $(REFDIR)/mel_r5_bad.gtf

SIMGFF   = $(PREREQDIR)/dsim-all-$(SIMVERSION).gff
SIMGTF   = $(REFDIR)/sim_good.gtf
SIMALLGTF= $(REFDIR)/sim_all.gtf
SIMBADGTF= $(REFDIR)/sim_bad.gtf

GENEMAPTABLE = $(PREREQDIR)/gene_map_table_fb_$(MELDATE).tsv


.SECONDARY:


all :  $(REFDIR)/mel_$(MELMAJORVERSION) $(REFDIR)/mel_$(MELVERSION) $(FPKMS) $(ANALYSIS_DIR)/ase_summary.tsv

genomes: Reference/Dmel/Genome $(SIMFASTA2) $(MELFASTA2) $(SECFASTA2)
	echo "Genomes Made"


# Read the per-project make-file
include config.make
include hybrids.make
include processing.make

$(ANALYSIS_DIR)/retabulate:
	touch $@

$(ANALYSIS_DIR)/summary.tsv : $(ANALYSIS_DIR)/retabulate MakeSummaryTable.py $(FPKMS) $(RUNCONFIG) $(ANALYSIS_DIR)/map_stats.tsv| $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
       --params $(RUNCONFIG) \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 63 \
	   --map-stats $(ANALYSIS_DIR)/map_stats.tsv \
	   --filename $(QUANT_FNAME) \
	   --key $(QUANT_KEY) \
	   --column $(QUANT_COL) \
		$(ANALYSIS_DIR) \
		\| tee $(ANALYSIS_DIR)/mst.log


$(ANALYSIS_DIR)/summary_fb.tsv : $(ANALYSIS_DIR)/retabulate MakeSummaryTable.py $(FPKMS) $(RUNCONFIG) $(ANALYSIS_DIR)/map_stats.tsv | $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
       --params $(RUNCONFIG) \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 65 \
	   --map-stats $(ANALYSIS_DIR)/map_stats.tsv \
	   --out-basename summary_fb \
	   --filename $(QUANT_FNAME) \
	   --key tracking_id \
	   --column $(QUANT_COL) \
		$(ANALYSIS_DIR) \
		\| tee $(ANALYSIS_DIR)/mst_fb.log

$(ANALYSIS_DIR)/summary_wasp.tsv : $(ANALYSIS_DIR)/retabulate MakeSummaryTable.py $$(subst genes,wasp_genes,$$(FPKMS)) $(RUNCONFIG) $(ANALYSIS_DIR)/map_stats.tsv | $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
       --params $(RUNCONFIG) \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 65 \
	   --map-stats $(ANALYSIS_DIR)/map_stats.tsv \
	   --out-basename summary_wasp \
	   --filename wasp_$(QUANT_FNAME) \
	   --key gene_short_name \
	   --column $(QUANT_COL) \
		$(ANALYSIS_DIR) \
		\| tee $(ANALYSIS_DIR)/mst_fb.log


$(ANALYSIS_DIR)/ase_summary_by_read.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,melsim_gene_ase_by_read.tsv,$$(FPKMS)) $(ANALYSIS_DIR)/map_stats.tsv
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename melsim_gene_ase_by_read.tsv \
			--column "ase_value" \
			--strip-low-reads 3000000 \
			--strip-on-unique \
			--strip-as-nan \
			--map-stats $(ANALYSIS_DIR)/map_stats.tsv \
			--key "gene" \
			--out-basename ase_summary_by_read \
			--float-format "%0.6f" \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/ase_summary_by_read_with_wasp.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,melsim_gene_ase_by_read_with_wasp.tsv,$$(FPKMS)) $(ANALYSIS_DIR)/map_stats.tsv
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename melsim_gene_ase_by_read_with_wasp.tsv \
			--column "ase_value" \
			--strip-low-reads 3000000 \
			--strip-on-unique \
			--strip-as-nan \
			--map-stats $(ANALYSIS_DIR)/map_stats.tsv \
			--key "gene" \
			--out-basename ase_summary_by_read_with_wasp \
			--float-format "%0.6f" \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/ase_summary_by_read_in_%subset.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,$$(*)subset/melsim_gene_ase_by_read.tsv,$$(FPKMS)) $(ANALYSIS_DIR)/map_stats.tsv
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename melsim_gene_ase_by_read.tsv \
			--in-subdirectory $*subset \
			--column "ase_value" \
			--strip-low-reads 90000 \
			--strip-on-unique \
			--strip-as-nan \
			--map-stats $(ANALYSIS_DIR)/map_stats.tsv \
			--key "gene" \
			--out-basename ase_summary_by_read \
			--float-format "%0.6f" \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/ase_summary.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,melsim_gene_ase.tsv,$$(FPKMS)) $(ANALYSIS_DIR)/map_stats.tsv
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename melsim_gene_ase.tsv \
			--strip-low-reads 3000000 \
			--strip-on-unique \
			--strip-as-nan \
			--map-stats $(ANALYSIS_DIR)/map_stats.tsv \
			--column "REF-ALT_RATIO" \
			--key "FEATURE" \
			--out-basename ase_summary \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/ase_summary_on_sim.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,simmel_gene_ase.tsv,$$(FPKMS))
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename simmel_gene_ase.tsv \
			--column "REF-ALT_RATIO" \
			--key "FEATURE" \
			--out-basename ase_summary_on_sim \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/ase_summary_on_sim_by_read.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,simmel_gene_ase_by_read.tsv,$$(FPKMS))
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename simmel_gene_ase.tsv \
			--column "REF-ALT_RATIO" \
			--key "FEATURE" \
			--out-basename ase_summary_on_sim_by_read \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/cds_ase_summary.tsv: $(ANALYSIS_DIR)/retabulate $$(subst genes.fpkm_tracking,melsim_cds_ase.tsv,$$(FPKMS))
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 6 \
	python MakeSummaryTable.py \
			--params Parameters/RunConfig.cfg \
			--filename melsim_cds_ase.tsv \
			--column "REF-ALT_RATIO" \
			--key "FEATURE" \
			--out-basename cds_ase_summary \
			$(ANALYSIS_DIR)

$(ANALYSIS_DIR)/map_stats.tsv: $$(subst genes.fpkm_tracking,assigned_dmelR.mapstats,$$(FPKMS))
	python GetMapStats.py \
		--params Parameters/RunConfig.cfg \
		-u \
		--count-all \
		-T \
		$(ANALYSIS_DIR)

%.mapstats: %.bam
	./qsubber $(QSUBBER_ARGS) -t 1 \
	python GetSingleMapStats.py $<


%/genes.fpkm_tracking : %/assigned_dmelR.bam $(MELGTF) $(MELFASTA2) $(MELBADGTF)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	touch $@
	./qsubber $(QSUBBER_ARGS)_$(*F) --load-module cufflinks -t 4 \
	cufflinks \
		--num-threads 8 \
		--output-dir $(@D) \
		--multi-read-correct \
		--frag-bias-correct $(MELFASTA2) \
		--GTF $(MELGTF) \
		--mask-file $(MELBADGTF) \
		$<

%/assigned_dmelR.bam : %/accepted_hits.bam AssignReads2.py
	samtools view -H $< \
		| grep -Pv 'SN:(?!dmel)' \
		> $(@D)/mel_only.header.sam
	samtools view -H $< \
		| grep -oP 'SN:....' \
		| cut -c 4- \
		| sort -u \
		> $(@D)/species_present
	ns=`wc -l $(@D)/species_present | cut -f 1`
	./qsubber $(QSUBBER_ARGS)_$(*F) --load-module samtools -t 3 \
	if [ `wc -l $(@D)/species_present | cut -d ' ' -f 1` -eq "1" ]\; then \
		samtools sort -@ 3 -o $@ $< \; \
	else \
		python AssignReads2.py $(@D)/accepted_hits.bam\; \
		samtools sort -o $@ $(@D)/assigned_dmel.bam \; \
		samtools view $(@D)/assigned_dmel_sorted.bam \
			\| cat $(@D)/mel_only.header.sam - \
			\| samtools view -bS -o $@ -\; \
		rm $(@D)/assigned_dmel_sorted.bam\; \
	fi
	samtools index $@

$(MELALLGTF): $(MELGFF) | $(REFDIR)
	gffread $< -C -E -T -o- | \
		awk '{print "dmel_"$$0}' > \
		$@

$(MELALLGTFR5): $(MELGFFR5) | $(REFDIR)
	gffread $< -C -E -T -o- | \
		awk '{print "dmel_"$$0}' > \
		$@

$(MELGTF): $(MELALLGTF) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		| grep 'gene_id' \
		> $@

$(MELGTFR5): $(MELALLGTFR5) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		| grep 'gene_id' \
		> $@

$(SIMGTF): $(SIMALLGTF) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		| grep -vP 'dsim_Scf_NODE_(108665)' \
		| grep 'gene_id' \
		> $@

$(SECGTF): $(SECALLGTF) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		| grep 'gene_id' \
		> $@

$(MELBADGTF): $(MELALLGTF) | $(REFDIR)
	cat $< \
		| grep -P '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> $@

$(MELBADGTFR5): $(MELALLGTFR5) | $(REFDIR)
	cat $< \
		| grep -P '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> $@

$(SIMALLGTF): $(SIMGFF) | $(REFDIR)
	gffread $< -C -E -T -o- | \
		awk '{print "dsim_"$$0}' > \
		$@

$(SIMBADGTF): $(SIMALLGTF) | $(REFDIR)
	cat $< \
		| grep -v 'gene_id' \
		> $@

$(MELFASTA): $(REFDIR)/mel_$(MELMAJORVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/fasta/dmel-all-chromosome-$(MELVERSION).fasta.gz
	gunzip --force $@.gz

$(SIMFASTA): $(REFDIR)/sim_$(SIMMAJORVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_$(SIMRELEASE)/fasta/dsim-all-chromosome-$(SIMVERSION).fasta.gz
	gunzip --force $@.gz

$(MELTRANSCRIPTS) : $(REFDIR)/mel_$(MELVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/fasta/dmel-all-transcript-$(MELVERSION).fasta.gz
	gunzip --force $@.gz


$(MELGFF): $(REFDIR)/mel_$(MELVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/gff/dmel-all-$(MELVERSION).gff.gz
	gunzip --force $@.gz

$(SIMGFF): $(REFDIR)/sim_$(SIMVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_$(SIMRELEASE)/gff/dsim-all-$(SIMVERSION).gff.gz
	gunzip --force $@.gz

$(MELFASTA2): $(MELFASTA) $(REFDIR)/mel_$(MELMAJORVERSION) | $(REFDIR)
	perl -pe 's/>/>dmel_/' $(MELFASTA) > $@

$(MELR5FASTA2): $(MELR5FASTA) | $(REFDIR)
	perl -pe 's/>/>dmel_/' $< > $@

$(SIMFASTA2): $(SIMFASTA) $(REFDIR)/sim_$(SIMMAJORVERSION) | $(REFDIR)
	perl -pe 's/>/>dsim_/' $(SIMFASTA) > $@

$(REFDIR)/Dmel/transcriptome : $(MELGTF) |  $(REFDIR)/Dmel
	tophat --GTF $(MELGTF) \
		--transcriptome-index $@ \
		$(REFDIR)/Dmel
	touch $@


$(ORTHOLOGS) : | $(PREREQDIR)
	wget -O $@.gz -i ftp.flybase.org/releases/FB$(MELDATE)/precomputed_files/genes/gene_orthologs_fb_$(MELDATE).tsv.gz
	gunzip --force $@.gz

$(REFDIR) :
	mkdir $@

$(PREREQDIR):
	mkdir $@

$(ANALYSIS_DIR):
	mkdir $@

##### MEL GENOMES ####
$(REFDIR)/Dmel/Genome : $(REFDIR)/mel_$(MELMAJORVERSION) | $(MELGTF)  $(REFDIR)/Dmel $(MELFASTA2) $(REFDIR)
	rm -rf $(@D)/_tmp
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dmel \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(MELFASTA2) \
		--sjdbGTFfile $(MELGTF)

$(REFDIR)/Dmel_unspliced/Genome : $(REFDIR)/mel_$(MELMAJORVERSION) | $(REFDIR)/Dmel_unspliced $(MELFASTA2) $(REFDIR)
	rm -rf $(@D)/_tmp
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dmel_unspliced \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(MELFASTA2) \

$(REFDIR)/Dmel : | $(REFDIR)
	mkdir $@

$(REFDIR)/Dmel_unspliced : | $(REFDIR)
	mkdir $@

##### SIM GENOMES ####
$(REFDIR)/Dsim/Genome : $(REFDIR)/sim_$(SIMMAJORVERSION) | $(SIMGTF)  $(REFDIR)/Dsim $(SIMFASTA2) $(REFDIR)
	rm -rf $(@D)/_tmp
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dsim \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(SIMFASTA2) \
		--sjdbGTFfile $(SIMGTF)

$(REFDIR)/Dsim_unspliced/Genome : $(REFDIR)/sim_$(SIMMAJORVERSION) | $(REFDIR)/Dsim_unspliced $(SIMFASTA2) $(REFDIR)
	rm -rf $(@D)/_tmp
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dsim_unspliced \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(SIMFASTA2)

$(REFDIR)/Dsim : | $(REFDIR)
	mkdir $@

$(REFDIR)/Dsim_unspliced : | $(REFDIR)
	mkdir $@


%/:
	mkdir $@

$(GENEMAPTABLE):
	wget ftp://ftp.flybase.net/releases/FB$(MELDATE)/precomputed_files/genes/$(notdir $(GENEMAPTABLE)).gz \
		-O $(GENEMAPTABLE).gz
	gunzip --force $(GENEMAPTABLE).gz

%_sorted.bam: %.bam
	./qsubber $(QSUBBER_ARGS)_$(*F) -t 1 \
	samtools sort $< $*_sorted
	samtools index $@

%.bam.bai: %.bam
	samtools index $<

$(REFDIR)/mel_$(MELVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/mel_$(MELMAJORVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/sim_$(SIMVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/sim_$(SIMMAJORVERSION): | $(REFDIR)
	touch $@


$(REFDIR)/d%_masked/done: $(REFDIR)/d%_masked.fasta
	python faSplitter.py $< $(@D)
	touch $@

%.bwt : %
	bwa index $*
	samtools faidx $*

%.fai : %
	samtools faidx $<

%.1.ebwt: $(basename $(basename %)).fasta
	bowtie-build --offrate 3 $< $(basename $(basename $@))

%.1.bt2: $(basename $(basename %)).fasta
	bowtie2-build --offrate 3 $< $(basename $(basename $@))

%.dict : %.fasta
	picard CreateSequenceDictionary R=$< O=$@

%.dict : %.fa
	picard CreateSequenceDictionary R=$< O=$@

%_transcriptome: %.1.bt2
	echo $(call uc,$(patsubst Reference/d%_prepend,%,$*))GTF
	tophat2 --transcriptome-index $@ \
			--GTF $($(call uc,$(patsubst Reference/d%_prepend,%,$*))GTF) \
			$*


$(REFDIR)/d%_masked.fasta: $(REFDIR)/d%_prepend.fasta
	trfBig $< $@

%.bw : %.bam
	python ChromSizes.py $<
	bamToBed -i $< -bed12 | bed12ToBed6 -i stdin | genomeCoverageBed -bga -i stdin -g $<.chromsizes > $(basename $@).bed
	bedSort $(basename $@).bed $(basename $@)_sorted.bed
	bedGraphToBigWig $(basename $@)_sorted.bed $<.chromsizes $@
	rm $(basename $@).bed $(basename $@)_sorted.bed

%.cxb : %.bam
	./qsubber --job-name $(@F)_cuffquant --queue batch --keep-temporary tmp -t 2 \
	cuffquant \
		--output-dir $(@D) \
		--num-thread 2 \
		--multi-read-correct \
		--mask-file $(MELBADGTF) \
		--frag-bias-correct $(MELFASTA2) \
		$(MELGTF) \
		$<
	mv $(@D)/abundances.cxb $@
	#cp $(*D)/$(*F)/abundances.cxb $@

Reference/tss: $(MELGTF)
	cat $< \
		| python FindTSSs.py \
	    | rev \
	    | uniq -f 1 \
	    | rev \
	    > $@

Reference/dmel_R5.gtf: prereqs/dmel-all-r5.57.gff
	gffread $< -C -E -T -o-  \
		> $@

Reference/tss_r5: Reference/dmel_R5.gtf
	cat $< \
		| python FindTSSs.py \
	    | rev \
	    | uniq -f 1 \
	    | rev \
	    > $@

%/unmapped_sample.fasta: %/Unmapped.out.mate1
	python SampleUnmapped.py $< $@

%/blasted.tsv : %/unmapped_sample.fasta
	./qsubber --job-name $(basename $(@D))_blast --load-module blast --keep-temporary tmp -t 4 --log-base $(basename $@) \
	blastn -db nt \
		-query $< \
		-outfmt \"6 qseqid staxids sskingdoms sblastnames sscinames\" \
		-num_threads 8  -max_target_seqs 1 \
		-out $@

%subset/accepted_hits.bam: $$(dir $$(*))/assigned_dmelR_wasp_dedup.sorted.bam
	mkdir -p $(dir $*)/subset$(notdir $*)M/
	ln -f -s subset$(notdir $*)M/ $(@D)
	./qsubber --job-name $(notdir $(*D))_subset$(notdir $*) --keep-temporary tmp -t 1 --log-base $(basename $@) \
	python ./SubSample.py $(dir $*)/assigned_dmelR_wasp_dedup.sorted.bam $(notdir $*)
