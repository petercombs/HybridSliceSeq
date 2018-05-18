$(REFDIR)/dgrp_%.vcf: prereqs/dgrp2.vcf
	grep -P '^(#|$*)' $< > $@

$(REFDIR)/melsim_snps.delta : $(MELFASTA2) $(SIMFASTA2)
	nucmer --prefix $(basename $@) $^


BDTNP_BIND_DIR=prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables

Reference/binding/bcd_r5_peaks: $(BDTNP_BIND_DIR)
	for TF in bcd cad D da dl gt hb hkb kni \
				kr mad med prd run shn slp1 sna tll twi z; do \
		cat $(BDTNP_BIND_DIR)/$${TF}_[1-3]_* \
			| awk 'NR > 1 {print $$2":"$$3".."$$4}' \
			> /tmp/$${TF}_peaks; \
		echo 'OldPeak	NewPeak' > Reference/binding/$${TF}_peaks; \
		curl -L  'http://flybase.org/cgi-bin/coord_converter.html' \
			-H 'Referer: http://flybase.org/static_pages/downloads/COORD.html' \
			-F species=dmel \
			-F inr=4\
			-F outr=5\
			-F saveas=File\
			-F ids="" \
			-F idfile="@/tmp/$${TF}_peaks" \
			-F .submit=Go \
			| grep -v '?' \
			>> Reference/binding/$${TF}_r5_peaks; \
		sleep 5; \
	done

Reference/binding/bcd_peaks: $(BDTNP_BIND_DIR)
	for TF in bcd cad D da dl gt hb hkb kni \
				kr mad med run shn slp1 sna tll twi z; do \
		cat $(BDTNP_BIND_DIR)/$${TF}_[1-3]_* \
			| awk 'NR > 1 {print $$2":"$$3".."$$4}' \
			> /tmp/$${TF}_peaks; \
		echo 'OldPeak	NewPeak' > Reference/binding/$${TF}_peaks; \
		curl -L  'http://flybase.org/cgi-bin/coord_converter.html' \
			-H 'Referer: http://flybase.org/static_pages/downloads/COORD.html' \
			-F species=dmel \
			-F inr=4\
			-F outr=6\
			-F saveas=File\
			-F ids="" \
			-F idfile="@/tmp/$${TF}_peaks" \
			-F .submit=Go \
			| grep -v '?' \
			>> Reference/binding/$${TF}_peaks; \
		sleep 5; \
	done

Reference/binding/bcd_%_peaks Reference/binding/hb_%_peaks Reference/binding/kr_%_peaks Reference/binding/prd_%_peaks: $(BDTNP_BIND_DIR)
	for TF in bcd cad D da dl gt hb hkb kni \
				kr mad med prd run shn slp1 sna tll twi z; do \
		cat `ls $(BDTNP_BIND_DIR)/$${TF}_{1,2,3,BQ,FQ}_*` \
			| sort -rnk  5 \
			| head -n $* \
			| awk 'NR > 1 {print $$2":"$$3".."$$4}' \
			> /tmp/$${TF}_peaks; \
		echo 'OldPeak	NewPeak' > Reference/binding/$${TF}_peaks; \
		curl -L  'http://flybase.org/cgi-bin/coord_converter.html' \
			-H 'Referer: http://flybase.org/static_pages/downloads/COORD.html' \
			-F species=dmel \
			-F inr=4\
			-F outr=6\
			-F saveas=File\
			-F ids="" \
			-F idfile="@/tmp/$${TF}_peaks" \
			-F .submit=Go \
			| grep -v '?' \
			>> Reference/binding/$${TF}_$*_peaks; \
		sleep 5; \
	done

Reference/binding/%_peaks.bed : Reference/binding/%_peaks
	cat $< \
		| cut -f 2 \
		| perl -pe 's/(:|\.\.)/\t/g' \
		| awk 'NR > 1' \
		> $@

%peaks_no_exons.bed: %peaks.bed $(REFDIR)/exons.bed
	module load bedtools
	bedtools subtract -a $< -b $(REFDIR)/exons.bed \
		| sort -k 1 -nk2 -nk 3 -u \
		> $@

%peaks_no_CDS.bed: %peaks.bed $(REFDIR)/CDS.bed
	module load bedtools
	bedtools subtract -a $< -b $(REFDIR)/CDS.bed \
			> $@

%_prepend.bed : %.bed
	cat $< | perl -pe 's/^/dmel_/' > $@

analysis/results/%bp.txt: Reference/$$(shell echo $$* | perl -pe 's/peak_genes.*/peaks_prepend/' ).bed
	bedtools window -w $(lastword $(subst _, ,$*)) -u -a $(MELGTF) -b $< \
		| awk '{print $$NF}' \
		| perl -pe 's/"//g' \
		| perl -pe 's/; *$$//' \
		| uniq | sort -u \
		> $@

analysis/results/%bp_tss.txt: Reference/$$(shell echo $$* | perl -pe 's/peak_genes.*/peaks_prepend/' ).bed
	bedtools window -w $(lastword $(subst _, ,$*)) -u -a Reference/tss.bed -b $< \
		| awk '{print $$4}' \
		| perl -pe 's/"//g' \
		| perl -pe 's/; *$$//' \
		| uniq | sort -u \
		> $@

$(REFDIR)/CDS.bed: $(REFDIR)/mel_good.gtf
	grep CDS $< \
			| perl -pe 's/FBtr[^"]*//' \
			| cut -c 6- | cut -f 1,4,5,9 \
			| sort -u -k1 -nk2 -nk3 \
			> $@


analysis/results/has_svase_%.bed : $(REFDIR)/%_peaks_no_exons_prepend.bed analysis/results/has_svase_prepend.bed
	module load bedtools
	bedtools window -w 5000 -b analysis/results/has_svase_prepend.bed -a $< -u \
		> $@

analysis/results/no_svase_%.bed: $(REFDIR)/%_peaks_no_exons_prepend.bed analysis/results/has_svase_%.bed analysis/results/no_svase_prepend.bed
	module load bedtools
	bedtools window -w 5000 -b analysis/results/no_svase_prepend.bed -a $< -u \
		| bedtools window -v -w 5000 -a - -b analysis/results/has_svase_$*.bed \
		> $@


%_snpcounts.bed : %.bed analysis/on_mel/melsim_variant.bed
	bedtools intersect -a $< -b analysis/on_mel/melsim_variant.bed -c > $@

%_snpcounts_merged.bed : %_snpcounts.bed
	cat $< | awk -v OFS="\t" '{print $$0,$$4 / ($$3 -  $$2)}' \
		| bedtools sort -i - \
		| bedtools merge -d 5000 -i - -c 5 -o max \
		> $@

%_ucsc.bed: %.bed
	echo track name='$(@F)' | cat - $< | perl -pe 's/^(dmel_)?/chr/' \
		| perl -pe 's/chrtrack/track/' \
		| grep -P "(chr(X|Y|2L|2R|3L|3R|4)\t|track)" > $@


analysis/targets/BDTNP_binds:
	mkdir -p $@

analysis/targets/BDTNP_binds/all: analysis/targets/BDTNP_binds/bcd.cm analysis/targets/BDTNP_binds/hb.cm analysis/targets/BDTNP_binds/cad.cm analysis/targets/BDTNP_binds/D.cm analysis/targets/BDTNP_binds/dl.cm analysis/targets/BDTNP_binds/ftz.cm analysis/targets/BDTNP_binds/gt.cm analysis/targets/BDTNP_binds/kr.cm analysis/targets/BDTNP_binds/shn.cm analysis/targets/BDTNP_binds/slp1.cm analysis/targets/BDTNP_binds/z.cm analysis/targets/BDTNP_binds/tll.cm analysis/targets/BDTNP_binds/sna.cm



analysis/targets/BDTNP_binds/%.cm: prereqs/BDTNP_pwms/%.cm
	cat $< \
			| awk 'BEGIN {maxnf=0};  \
				   {print $$0; maxnf=((maxnf > NF)?maxnf:NF)}; \
				   END {printf "N|"; for (i = 1; i < maxnf; i++) printf "\t0"; printf "\n"}' \
			> $@
