$(REFDIR)/dgrp_%.vcf: prereqs/dgrp2.vcf
	grep -P '^(#|$*)' $< > $@

$(REFDIR)/melsim_snps.delta : $(MELFASTA2) $(SIMFASTA2)
	nucmer --prefix $(basename $@) $^
