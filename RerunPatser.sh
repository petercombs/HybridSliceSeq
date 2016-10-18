#!/usr/bin/bash

#source /etc/profile.d/modules.sh
PATH=$PATH:$HOME/meme/bin
module load bedtools

echo "Getting FASTAs"
bedtools getfasta -s -name -fi Reference/dmel_prepend.fasta -fo Reference/mel_hb_regregions.fasta -bed Reference/hb_with_imgs.bed
bedtools getfasta -s -name -fi Reference/dsim_prepend.fasta -fo Reference/sim_hb_regregions.fasta -bed Reference/sim_hb_regregions.bed
bedtools getfasta -s -name -fi Reference/dmel_prepend.fasta -fo Reference/mel_Kr_regregions.fasta -bed Reference/mel_Kr_regregions.bed
bedtools getfasta -s -name -fi Reference/dsim_prepend.fasta -fo Reference/sim_Kr_regregions.fasta -bed Reference/sim_Kr_regregions.bed
bedtools getfasta -s -name -fi Reference/dmel_prepend.fasta -fo Reference/mel_pinocchio_regregions.fasta -bed Reference/mel_pinocchio.bed
bedtools getfasta -s -name -fi Reference/dsim_prepend.fasta -fo Reference/sim_pinocchio_regregions.fasta -bed Reference/sim_pinocchio.bed

cat Reference/{mel,sim}_hb_regregions.fasta Reference/sim_hkb_2016_10_12{,F,_manstitch}.fa > Reference/melsim_hb_regregions.fasta
cat Reference/{mel,sim}_Kr_regregions.fasta > Reference/melsim_Kr_regregions.fasta
cat Reference/{mel,sim}_pinocchio_regregions.fasta > Reference/melsim_pinocchio_regregions.fasta

cat Reference/sim_hb_regregions.fasta Reference/sim_hkb_2016_10_12{,F,_manstitch}.fa  > Reference/sim_hb_regregions_plussanger.fasta

#nucmer -g 400 -b 400 -p Reference/simsec_both_promoters Reference/sim_eloF_promoters.fasta Reference/sec_eloF_promoters.fasta
#nucmer -g 400 -p Reference/melsec_both_promoters Reference/mel_eloF_promoters.fasta Reference/sec_eloF_promoters.fasta
needleall -aseq Reference/mel_hb_regregions.fasta \
	-bseq Reference/sim_hb_regregions_plussanger.fasta \
	-aformat3 srspair \
	-gapopen 10.0 -gapextend 0.5 \
	-outfile Reference/melsim_hb_regregions.needleall

needleall -aseq Reference/mel_Kr_regregions.fasta \
	-bseq Reference/sim_Kr_regregions.fasta \
	-aformat3 srspair \
	-gapopen 10.0 -gapextend 0.5 \
	-outfile Reference/melsim_Kr_regregions.needleall

needleall -aseq Reference/mel_pinocchio_regregions.fasta \
	-bseq Reference/sim_pinocchio_regregions.fasta \
	-aformat3 srspair \
	-gapopen 10.0 -gapextend 0.5 \
	-outfile Reference/melsim_pinocchio_regregions.needleall

echo "Getting PSEQs"
python FastaToPseq.py Reference/mel_hb_regregions.fasta
python FastaToPseq.py Reference/sim_hb_regregions.fasta
python FastaToPseq.py Reference/mel_Kr_regregions.fasta
python FastaToPseq.py Reference/sim_Kr_regregions.fasta

echo "Mapping binding sites"
#parallel 'patser-v3e -c -s -lp -6.3 -w -m {} -v -f ~/HybridSliceSeq/Reference/mel_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/mel_{/.}.txt' ::: analysis/targets/OnTheFly/*
#parallel 'patser-v3e -M -6 -c -s -lp -6.3 -m {} -f ~/HybridSliceSeq/Reference/mel_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/mel_{/.}.txt' ::: analysis/targets/BDTNP_binds/*
#parallel 'patser-v3e -c -s -lp -6.3 -w -m {} -v -f ~/HybridSliceSeq/Reference/sim_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/sim_{/.}.txt' ::: analysis/targets/OnTheFly/*
#parallel  'patser-v3e -M -6 -c -s -lp -6.3 -m {} -f ~/HybridSliceSeq/Reference/sim_hb_regregions.pseq -a ~/oenocytes/Reference/alphabet_mel > analysis/targets/hb/sim_{/.}.txt' ::: analysis/targets/BDTNP_binds/*

fimo -oc analysis/targets/hb/mel --thresh 1e-3  Reference/all_meme_filtered.meme Reference/mel_hb_regregions.fasta 2> /dev/null
fimo -oc analysis/targets/hb/sim --thresh 1e-3  Reference/all_meme_filtered.meme Reference/sim_hb_regregions_plussanger.fasta 2> /dev/null
fimo -oc analysis/targets/Kr/mel --thresh 1e-3  Reference/all_meme_filtered.meme Reference/mel_Kr_regregions.fasta 2> /dev/null
fimo -oc analysis/targets/Kr/sim --thresh 1e-3  Reference/all_meme_filtered.meme Reference/sim_Kr_regregions.fasta 2> /dev/null
#fimo -oc analysis/targets/pinocchio/mel --thresh 1e-3  Reference/all_meme_filtered.meme Reference/mel_pinocchio_regregions.fasta 2> /dev/null
#fimo -oc analysis/targets/pinocchio/sim --thresh 1e-3  Reference/all_meme_filtered.meme Reference/sim_pinocchio_regregions.fasta 2> /dev/null

TFS="OTF0070.1 hb  kni gt D FBgn0001325_4 tll hkb  FBgn0000251_3 FBgn0003448_3 twi" #Med dl"
TF_NAMES="bcd  hb  kni gt D Kr            tll hkb cad            sna           twi" #med dl"
MATCH_DIM=0.7
TFS="OTF0070.1 hkb FBgn0003448_3"
TF_NAMES="bcd hkb  sna"
MATCH_DIM=0.0

#
python PatserAlignToSVG.py -v -X 100 --show-alignment -x .75 --y-sep 60 --y-scale 3.5 --bar-width 10 --comp1 mel --comp2 sim \
        --match-dim ${MATCH_DIM} \
        --meme-suite \
        --tf ${TFS} \
        --tf-names ${TF_NAMES} \
        --sequence --fasta Reference/melsim_hb_regregions.fasta \
        --needleall Reference/melsim_hb_regregions.needleall \
        analysis/targets/hb


python PatserAlignToSVG.py -v -X 100 --show-alignment -x .75 --y-sep 60 --y-scale 3.5 --bar-width 10 --comp1 mel --comp2 sim \
        --match-dim ${MATCH_DIM} \
        --meme-suite \
        --tf ${TFS} \
        --tf-names ${TF_NAMES} \
        --sequence --fasta Reference/melsim_Kr_regregions.fasta \
        --needleall Reference/melsim_Kr_regregions.needleall \
        analysis/targets/Kr

python PatserAlignToSVG.py -v -X 100 --show-alignment -x .75 --y-sep 60 --y-scale 3.5 --bar-width 10 --comp1 mel --comp2 sim \
        --match-dim ${MATCH_DIM} \
        --meme-suite \
        --tf ${TFS} \
        --tf-names ${TF_NAMES} \
        --sequence --fasta Reference/melsim_pinocchio_regregions.fasta \
        --needleall Reference/melsim_pinocchio_regregions.needleall \
        analysis/targets/pinocchio
