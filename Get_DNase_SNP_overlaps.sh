#!/usr/bin/env bash

module load bedtools

# Filter out no-svASE genes near svASE genes
bedtools window -w 5000 -b analysis/results/has_svase.bed -a analysis/results/no_svase.bed -u > analysis/results/has_near_no_svase.bed 
bedtools intersect -v -a analysis/results/no_svase.bed -b analysis/results/has_near_no_svase.bed > analysis/results/no_svase_clean.bed

# Find DNase peaks near svASE and nosvASE genes
bedtools window -w 5000 -b analysis/results/no_svase_clean.bed -a Reference/DNase_R6.bed -u > analysis/results/no_svase_dnase.bed
bedtools window -w 5000 -b analysis/results/has_svase.bed -a Reference/DNase_R6.bed -u > analysis/results/has_svase_dnase.bed

# Zelda Peaks
bedtools window -w 5000 -b analysis/results/no_svase_clean.bed -a Reference/li13_r6.bed -u \
	| bedtools subtract -a - -b Reference/exons.bed \
	> analysis/results/no_svase_zld.bed
bedtools window -w 5000 -b analysis/results/has_svase.bed -a Reference/li13_r6.bed -u \
	| bedtools subtract -a - -b Reference/exons.bed \
	> analysis/results/has_svase_zld.bed

# ORegAnno features
bedtools window -w 5000 -b analysis/results/no_svase_clean.bed -a Reference/oreganno_tfbs.bed -u \
	| bedtools subtract -a - -b Reference/exons.bed \
	> analysis/results/no_svase_tfbs.bed
bedtools window -w 5000 -b analysis/results/has_svase.bed -a Reference/oreganno_tfbs.bed -u \
	| bedtools subtract -a - -b Reference/exons.bed \
	> analysis/results/has_svase_tfbs.bed
bedtools window -w 5000 -b analysis/results/no_svase_clean.bed -a Reference/oreganno_reg.bed -u \
	| bedtools subtract -a - -b Reference/exons.bed \
	> analysis/results/no_svase_reg.bed
bedtools window -w 5000 -b analysis/results/has_svase.bed -a Reference/oreganno_reg.bed -u \
	| bedtools subtract -a - -b Reference/exons.bed \
	> analysis/results/has_svase_reg.bed


# Remove DNase peaks that are near genes with and without svASE
bedtools intersect -a analysis/results/no_svase_dnase.bed -b analysis/results/has_svase_dnase.bed -v \
	| bedtools subtract -a - -b Reference/exons.bed \
	>  analysis/results/no_svase_dnase_clean.bed
bedtools intersect -a analysis/results/no_svase_zld.bed -b analysis/results/has_svase_zld.bed -v \
	| bedtools subtract -a - -b Reference/exons.bed \
	>  analysis/results/no_svase_zld_clean.bed


# Count number of SNPs overlapping each DNase peak
bedtools intersect -a analysis/results/has_svase_dnase.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/has_svase_dnase_snps.bed
bedtools intersect -a analysis/results/no_svase_dnase_clean.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/no_svase_dnase_snps.bed

bedtools intersect -a analysis/results/has_svase_zld.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/has_svase_zld_snps.bed
bedtools intersect -a analysis/results/no_svase_zld_clean.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/no_svase_zld_snps.bed

bedtools intersect -a analysis/results/has_svase_tfbs.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/has_svase_tfbs_snps.bed
bedtools intersect -a analysis/results/no_svase_tfbs.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/no_svase_tfbs_snps.bed

bedtools intersect -a analysis/results/has_svase_reg.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/has_svase_reg_snps.bed
bedtools intersect -a analysis/results/no_svase_reg.bed -b analysis/on_mel/melsim_variant_noprepend.bed -c > analysis/results/no_svase_reg_snps.bed

