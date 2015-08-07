include gmsl

wittkopp_reads:  sequence/SRR869573 sequence/SRR869575 sequence/SRR869590 sequence/SRR869596 sequence/SRR869597 sequence/SRR869598 sequence/SRR869599 sequence/SRR869600 sequence/SRR869601 sequence/SRR869602 sequence/SRR869844 sequence/SRR869603 sequence/SRR869905 sequence/SRR869592
	echo "All Downloaded"

sequence/%:
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_1.fastq.gz
	wget -c -P sequence ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$(call substr,$*,1,6)/$*/$*_2.fastq.gz
	touch $@


