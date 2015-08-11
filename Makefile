# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis

# Reference FASTA and GFF files from FlyBase
MELRELEASE = r6.06_FB2015_03
SIMRELEASE = r2.01_FB2015_01
SECRELEASE = r1.3_FB2015_01

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
SIMFASTA = $(PREREQDIR)/dsim-all-chromosome-$(SIMVERSION).fasta
SECFASTA = $(PREREQDIR)/dsec-all-chromosome-$(SECVERSION).fasta

REFDIR = Reference

MELFASTA2= $(REFDIR)/dmel_prepend.fasta
SIMFASTA2= $(REFDIR)/dsim_prepend.fasta
SECFASTA2= $(REFDIR)/dsec_prepend.fasta

ORTHOLOGS = $(PREREQDIR)/gene_orthologs_fb_$(MELDATE).tsv

MELGFF   = $(PREREQDIR)/dmel-all-$(MELVERSION).gff
MELGTF   = $(REFDIR)/mel_good.gtf
MELALLGTF   = $(REFDIR)/mel_all.gtf
MELBADGTF   = $(REFDIR)/mel_bad.gtf

SIMGFF   = $(PREREQDIR)/dsim-all-$(SIMVERSION).gff
SIMGTF   = $(REFDIR)/sim_good.gtf
SIMALLGTF= $(REFDIR)/sim_all.gtf
SIMBADGTF= $(REFDIR)/sim_bad.gtf

SECGFF   = $(PREREQDIR)/dsec-all-$(SECVERSION).gff
SECGTF   = $(REFDIR)/sec_good.gtf
SECALLGTF= $(REFDIR)/sec_all.gtf
SECBADGTF= $(REFDIR)/sec_bad.gtf

GENEMAPTABLE = $(PREREQDIR)/gene_map_table_fb_$(MELDATE).tsv


all :  $(REFDIR)/mel_$(MELMAJORVERSION) $(REFDIR)/mel_$(MELVERSION)

genomes: Reference/Dmel/Genome $(SIMFASTA2) $(MELFASTA2) $(SECFASTA2)
	echo "Genomes Made"


# Read the per-project make-file
include config.make
include analyze.make

$(MELALLGTF): $(MELGFF) | $(REFDIR)
	gffread $< -E -T -o- | \
		awk '{print "dmel_"$$0}' > \
		$@

$(MELGTF): $(MELALLGTF) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> $@

$(SIMGTF): $(SIMALLGTF) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> $@

$(SECGTF): $(SECALLGTF) | $(REFDIR)
	cat $< \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> $@

$(MELBADGTF): $(MELALLGTF) | $(REFDIR)
	cat $< \
		| grep -P '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> $@
$(SIMALLGTF): $(SIMGFF) | $(REFDIR)
	gffread $< -E -T -o- | \
		awk '{print "dsim_"$$0}' > \
		$@

$(SECALLGTF): $(SECGFF) | $(REFDIR)
	gffread $< -E -T -o- | \
		awk '{print "dsec_"$$0}' > \
		$@

$(MELFASTA): $(REFDIR)/mel_$(MELMAJORVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/fasta/dmel-all-chromosome-$(MELVERSION).fasta.gz
	gunzip --force $@.gz

$(SIMFASTA): $(REFDIR)/sim_$(SIMMAJORVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_$(SIMRELEASE)/fasta/dsim-all-chromosome-$(SIMVERSION).fasta.gz
	gunzip --force $@.gz

$(SECFASTA): $(REFDIR)/sec_$(SECMAJORVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_sechellia/dsec_$(SECRELEASE)/fasta/dsec-all-chromosome-$(SECVERSION).fasta.gz
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

$(SECGFF): $(REFDIR)/sec_$(SECVERSION) | $(REFDIR) $(PREREQDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_sechellia/dsec_$(SECRELEASE)/gff/dsec-all-$(SECVERSION).gff.gz
	gunzip --force $@.gz

$(MELFASTA2): $(MELFASTA) $(REFDIR)/mel_$(MELMAJORVERSION) | $(REFDIR)
	perl -pe 's/>/>dmel_/' $(MELFASTA) > $@

$(SIMFASTA2): $(SIMFASTA) $(REFDIR)/sim_$(SIMMAJORVERSION) | $(REFDIR)
	perl -pe 's/>/>dsim_/' $(SIMFASTA) > $@

$(SECFASTA2): $(SECFASTA) $(REFDIR)/sec_$(SECMAJORVERSION) | $(REFDIR)
	perl -pe 's/>/>dsec_/' $(SECFASTA) > $@
	
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
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dmel \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(MELFASTA2) \
		--sjdbGTFfile $(MELGTF)

$(REFDIR)/Dmel_unspliced/Genome : $(REFDIR)/mel_$(MELMAJORVERSION) | $(REFDIR)/Dmel_unspliced $(MELFASTA2) $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dmel_unspliced \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(MELFASTA2) \

$(REFDIR)/Dmel : | $(REFDIR)
	mkdir $@

$(REFDIR)/Dmel_unspliced : | $(REFDIR)
	mkdir $@

##### SIM GENOMES ####
$(REFDIR)/Dsim/Genome : $(REFDIR)/sim_$(SIMMAJORVERSION) | $(SIMGTF)  $(REFDIR)/Dsim $(SIMFASTA2) $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dsim \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(SIMFASTA2) \
		--sjdbGTFfile $(SIMGTF)

$(REFDIR)/Dsim_unspliced/Genome : $(REFDIR)/sim_$(SIMMAJORVERSION) | $(REFDIR)/Dsim_unspliced $(SIMFASTA2) $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dsim_unspliced \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(SIMFASTA2) 

$(REFDIR)/Dsim : | $(REFDIR)
	mkdir $@

$(REFDIR)/Dsim_unspliced : | $(REFDIR)
	mkdir $@

##### SEC GENOMES ####
$(REFDIR)/Dsec/Genome : $(REFDIR)/sec_$(SECMAJORVERSION) | $(SECGTF)  $(REFDIR)/Dsec $(SECFASTA2) $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dsec \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(SECFASTA2) \
		--sjdbGTFfile $(SECGTF)

$(REFDIR)/Dsec_unspliced/Genome : $(REFDIR)/sec_$(SECMAJORVERSION) | $(REFDIR)/Dsec_unspliced $(SECFASTA2) $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dsec_unspliced \
		--outTmpDir $(@D)/_tmp \
		--genomeFastaFiles $(SECFASTA2) 

$(REFDIR)/Dsec : | $(REFDIR)
	mkdir $@

$(REFDIR)/Dsec_unspliced : | $(REFDIR)
	mkdir $@

$(GENEMAPTABLE):
	wget ftp://ftp.flybase.net/releases/FB$(MELDATE)/precomputed_files/genes/$(notdir $(GENEMAPTABLE)).gz \
		-O $(GENEMAPTABLE).gz
	gunzip --force $(GENEMAPTABLE).gz

%_sorted.bam: %.bam
	touch $@
	samtools sort $< $*_sorted 
	samtools index $@

$(REFDIR)/mel_$(MELVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/mel_$(MELMAJORVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/sim_$(SIMVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/sim_$(SIMMAJORVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/sec_$(SECVERSION): | $(REFDIR)
	touch $@

$(REFDIR)/sec_$(SECMAJORVERSION): | $(REFDIR)
	touch $@
