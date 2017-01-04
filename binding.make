TFS=OTF0070.1 hb  kni gt D FBgn0001325_4 tll hkb  FBgn0000251_3 FBgn0003448_3 twi #Med dl"
TF_NAMES=bcd  hb  kni gt D Kr            tll hkb cad            sna           twi #med dl"
MATCH_DIM=0.7


#TFS=OTF0070.1 hkb FBgn0003448_3
#TF_NAMES=bcd hkb  sna
#MATCH_DIM=0.7

.SECONDEXPANSION:

.SECONDARY:


analysis/targets/%/aligned.needleall.svg: analysis/targets/%/aligned.needleall analysis/targets/%/mel/fimo.txt analysis/targets/%/sim/fimo.txt analysis/targets/%/melsim.fasta PatserAlignToSVG.py
	python PatserAlignToSVG.py -v -X 100 \
		--show-alignment \
		-x .75 --y-sep 60 --y-scale 3.5 \
		--bar-width 10 \
		--comp1 mel --comp2 sim \
		--match-dim $(MATCH_DIM) \
		--meme-suite \
		--tf $(TFS) \
		--tf-names $(TF_NAMES) \
		--sequence --fasta $(@D)/melsim.fasta \
		--needleall $(@D)/aligned.needleall \
		$(@D)

analysis/targets/%/melyak.needleall.svg: analysis/targets/%/melyak.needleall analysis/targets/%/mel/fimo.txt analysis/targets/%/yak/fimo.txt analysis/targets/%/melyak.fasta PatserAlignToSVG.py
	python PatserAlignToSVG.py -v -X 100 \
		--show-alignment \
		-x .75 --y-sep 60 --y-scale 3.5 \
		--bar-width 10 \
		--comp1 mel --comp2 yak \
		--match-dim $(MATCH_DIM) \
		--meme-suite \
		--tf $(TFS) \
		--tf-names $(TF_NAMES) \
		--sequence --fasta $(@D)/melyak.fasta \
		--needleall $< \
		$(@D)

analysis/targets/%/fimo.txt: analysis/targets/%.fasta
	mkdir -p $(@D)
	fimo -oc $(@D) --thresh 1e-3  Reference/all_meme_filtered.meme $< 2> $@.log



analysis/targets/%/melsim.fasta: analysis/targets/%/mel.fasta analysis/targets/%/sim.fasta
	cat $^ > $@

analysis/targets/%/melyak.fasta: analysis/targets/%/mel.fasta analysis/targets/%/yak.fasta
	cat $^ > $@


analysis/targets/%/mel.fasta : analysis/targets/%/mel.bed
	module load bedtools && \
	bedtools getfasta -s -name -fi Reference/d$(basename $(@F))_prepend.fasta -fo $@ -bed $<

analysis/targets/%/sim.fasta: analysis/targets/%/sim.bed
	module load bedtools && \
	bedtools getfasta -s -name -fi Reference/d$(basename $(@F))_prepend.fasta -fo $@ -bed $<

analysis/targets/%/yak.fasta: analysis/targets/%/yak.bed
	module load bedtools && \
	bedtools getfasta -s -name -fi Reference/dyak_prepend.fa -fo $@ -bed $<

analysis/targets/%/sim.bed: analysis/targets/%/mel.fasta
	#                   1      2     3      4      5     6       7     8    9    10     11
	module load blast bioawk && \
		blastn -db Reference/dsim_prepend.fasta \
		-outfmt "6 sseqid sstart send qseqid evalue sstrand length qlen slen qstart qend" \
		-gapextend 0\
		-query $< \
		| bioawk -t '$$10 > 1 && $$6 ~ /minus/ {$$2 += $$10 + 1}; \
		             $$10 > 1 && $$6 ~ /plus/ {$$2 -= $$10 + 1}; \
			     $$11 < $$8 && $$6 ~ /minus/ {$$3 -= ($$8 - $$11) + 1}; \
			     $$2 > $$3 {gsub("mel", "sim", $$4); print $$1,$$3,$$2+1,$$4,$$7/($$8+1),"-", $$3, $$2 }; \
		             $$2 < $$3 {gsub("mel", "sim", $$4); print $$1,$$2,$$3+1,$$4,$$7/($$8+1),"+", $$2, $$3 }; '\
		> $@

analysis/targets/%/yak.bed: analysis/targets/%/mel.fasta
	module load blast bioawk && \
		blastn -db Reference/yakpse \
		-outfmt "6 sseqid sstart send qseqid evalue sstrand length qlen slen qstart qend" \
		-gapextend 0\
		-max_target_seqs 1\
		-query $< \
		| awk '!_[$$4]++' \
		| bioawk -t '$$10 > 1 && $$6 ~ /minus/ {$$2 += $$10 + 1}; \
		             $$10 > 1 && $$6 ~ /plus/ {$$2 -= $$10 + 1}; \
			     $$11 < $$8 && $$6 ~ /minus/ {$$3 -= ($$8 - $$11) + 1}; \
			     $$2 > $$3 {gsub("mel", "yak", $$4); print $$1,$$3,$$2+1,$$4,$$7/($$8+1),"-", $$3, $$2 }; \
		             $$2 < $$3 {gsub("mel", "yak", $$4); print $$1,$$2,$$3+1,$$4,$$7/($$8+1),"+", $$2, $$3 }; '\
		> $@

analysis/targets/%/aligned.needleall: analysis/targets/%/mel.fasta analysis/targets/%/sim.fasta
	module load EMBOSS && \
	needleall -aseq analysis/targets/$*/mel.fasta \
		-bseq analysis/targets/$*/sim.fasta \
		-aformat3 srspair \
		-gapopen 10.0 -gapextend 0.5 \
		-outfile $@

analysis/targets/%/melyak.needleall: analysis/targets/%/mel.fasta analysis/targets/%/yak.fasta
	module load EMBOSS && \
	needleall -aseq analysis/targets/$*/mel.fasta \
		-bseq analysis/targets/$*/yak.fasta \
		-aformat3 srspair \
		-gapopen 10.0 -gapextend 0.5 \
		-outfile $@

analysis/targets/%/mel.bed: | analysis/targets/%/
	module load bedtools && \
		grep '"$*"' Reference/mel_r5_good.gtf \
		| bedtools sort \
		| bedtools merge  \
		| bedtools window -b - -a Reference/dnase_peaks_prepend.bed \
		-u -w 5000 \
		| bioawk -t '{print $$0, "mel_$*_" NR, "0", ".", $$2, $$3}' \
		> $@


analysis/targets/%/:
	mkdir $@


