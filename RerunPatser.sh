#!/usr/bin/bash 

source /etc/profile.d/modules.sh
module load bedtools

echo "Getting FASTAs"
bedtools getfasta -s -name -fi Reference/dmel_prepend.fasta -fo Reference/mel_hb_regregions.fasta -bed Reference/hb_with_imgs.bed
bedtools getfasta -s -name -fi Reference/dsim_prepend.fasta -fo Reference/sim_hb_regregions.fasta -bed Reference/sim_hb_regregions.bed

cat Reference/{mel,sim}_hb_regregions.fasta > Reference/melsim_hb_regregions.fasta

#nucmer -g 400 -b 400 -p Reference/simsec_both_promoters Reference/sim_eloF_promoters.fasta Reference/sec_eloF_promoters.fasta
#nucmer -g 400 -p Reference/melsec_both_promoters Reference/mel_eloF_promoters.fasta Reference/sec_eloF_promoters.fasta
#needleall -aseq Reference/mel_hb_regregions.fasta \
#	-bseq Reference/sim_hb_regregions.fasta \
#	-aformat3 srspair \
#	-gapopen 10.0 -gapextend 0.5 \
#	-outfile Reference/melsim_hb_regregions.needleall


echo "Getting PSEQs"
python FastaToPseq.py Reference/mel_hb_regregions.fasta
python FastaToPseq.py Reference/sim_hb_regregions.fasta

echo "Mapping binding sites"
#parallel 'patser-v3e -c -s -lp -6.3 -w -m {} -v -f ~/HybridSliceSeq/Reference/mel_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/mel_{/.}.txt' ::: analysis/targets/OnTheFly/*
#parallel 'patser-v3e -M -6 -c -s -lp -6.3 -m {} -f ~/HybridSliceSeq/Reference/mel_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/mel_{/.}.txt' ::: analysis/targets/BDTNP_binds/*
#parallel 'patser-v3e -c -s -lp -6.3 -w -m {} -v -f ~/HybridSliceSeq/Reference/sim_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/sim_{/.}.txt' ::: analysis/targets/OnTheFly/*
#parallel  'patser-v3e -M -6 -c -s -lp -6.3 -m {} -f ~/HybridSliceSeq/Reference/sim_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/sim_{/.}.txt' ::: analysis/targets/BDTNP_binds/*

#fimo -oc analysis/targets/hb/mel --thresh 1e-3  prereqs/fly_factor_survey.meme Reference/mel_hb_regregions.fasta 2> /dev/null
#fimo -oc analysis/targets/hb/sim --thresh 1e-3  prereqs/fly_factor_survey.meme Reference/sim_hb_regregions.fasta 2> /dev/null

TFS="OTF0070.1 hb  kni gt D FBgn0001325_4 tll hkb  FBgn0000251_3 FBgn0003448_3 twi" #Med dl"
TF_NAMES="bcd  hb  kni gt D Kr            tll hkb cad            sna           twi" #med dl"

#
echo python PatserAlignToSVG.py -X 100 --show-alignment -x .75 --y-sep 60 --y-scale 3.5 --bar-width 10 --comp1 mel --comp2 sim \
	--tf ${TFS} \
	--tf-names ${TF_NAMES} \
	--sequence --fasta Reference/melsim_hb_regregions.fasta \
	--needleall Reference/melsim_hb_regregions.needleall \
	analysis/targets/hb

python PatserAlignToSVG.py -X 100 --show-alignment -x .75 --y-sep 60 --y-scale 3.5 --bar-width 10 --comp1 mel --comp2 sim \
        --meme-suite \
        --tf ${TFS} \
        --tf-names ${TF_NAMES} \
        --sequence --fasta Reference/melsim_hb_regregions.fasta \
        --needleall Reference/melsim_hb_regregions.needleall \
        analysis/targets/hb
