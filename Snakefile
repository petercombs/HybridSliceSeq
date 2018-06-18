from os import path

mel_release = "r5.57_FB2014_03"
sim_release = "r2.01_FB2016_01"

mel_version, mel_date = mel_release.split('_', 1)
sim_version, sim_date = sim_release.split('_', 1)

num_mel_windows = 17
dates = {'mel': mel_date, 'sim': sim_date}
versions = {'mel': mel_version, 'sim': sim_version}

analysis_dir = 'analysis_godot'
hybrid_dir = 'analysis_godot/on_mel'
mel_gtf = 'Reference/mel_good.gtf'
mel_bad_gtf = 'Reference/mel_bad.gtf'
mel_fasta = 'Reference/dmel_prepend.fasta'

variants = path.join(hybrid_dir, 'melsim_variant.bed')

localrules: all, makedir, all_files_per_sample

configfile: "Parameters/config.json"

mxs_r1 = ['melXsim_cyc14C_rep1_sl{:02d}'.format(i) for i in range(4, 30)]
mxs_r2 = ['melXsim_cyc14C_rep2_sl{:02d}'.format(i) for i in range(7, 33)]
mxs_r3 = ['melXsim_cyc14C_rep3_sl{:02d}'.format(i) for i in range(5, 32)]
sxm_r1 = ['simXmel_cyc14C_rep1_sl{:02d}'.format(i) for i in range(37, 62)]
sxm_r2 = ['simXmel_cyc14C_rep2_sl{:02d}'.format(i) for i in range(3, 30)]
mxm_r1 = ['melXmel_cyc14C_sl{:02d}'.format(i) for i in range(7, 34)]
sxs_r1 = ['simXsim_cyc14C_sl{:02d}'.format(i) for i in range(7, 34)]

#samples = mxs_r1 + mxs_r2 + mxs_r3 + sxm_r1 + sxm_r2 + mxm_r1 + sxs_r1
samples = [sample for sample in config['samples'] if 'gdna' not in sample]

module = '''module () {
        eval `$LMOD_CMD bash "$@"`
        }'''

def getreads(readnum):
    if readnum == 0:
        formatstr = 'sequence/{}.fastq.gz'
    else:
        formatstr = 'sequence/{{}}_{}.fastq.gz'.format(readnum)
    def retfun(wildcards):
        return [formatstr.format(srr)
                for srr in config['samples'][path.basename(wildcards.sample)]]
    return retfun

def getreads_nowc(readnum):
    if readnum == 0:
        formatstr = 'sequence/{}.fastq.gz'
    else:
        formatstr = 'sequence/{{}}_{}.fastq.gz'.format(readnum)
    def retfun(sample):
        return [formatstr.format(srr)
                for srr in config['samples'][path.basename(sample)]]
    return retfun

def getreadscomma(readnum):
    if readnum == 0:
        formatstr = 'sequence/{}.fastq.gz'
    else:
        formatstr = 'sequence/{{}}_{}.fastq.gz'.format(readnum)
    def retfun(wildcards):
        return ",".join([formatstr.format(srr)
                        for srr in config['samples'][path.basename(wildcards.sample)]])
    return retfun


rule all:
    input:
        path.join(analysis_dir, "summary.tsv"),
        path.join(analysis_dir, "ase_summary_by_read.tsv")

rule all_reads:
    input:
        [reads for sample in samples for reads in getreads_nowc(1)(sample)]


rule print_config:
    run:
        print(config)

tfs = "OTF0070.1 hb  kni gt D FBgn0001325_4 tll hkb  FBgn0000251_3 FBgn0003448_3 twi OTF0532.1 OTF0478.1 OTF0181.1 FBgn0013753 FBgn0001204"
tf_names = "bcd  hb  kni gt D Kr            tll hkb cad            sna           twi zen       Doc2      rho/pnt   run/Bgb     hkb_FFS"

tf_dict = {
    'bcd': 'OTF0070.1',
    'hb': 'hb',
    'kni': 'kni',
    'gt': 'gt',
    'D':  'D',
    'Kr': 'FBgn0001325_4',
    'tll': 'tll',
    'hkb': 'hkb',
    'cad': 'FBgn0000251_3',
    'sna': 'FBgn0003448_3',
    'twi': 'twi',
    'zen': 'OTF0532.1',
    'Doc2': 'OTF0478.1',
    'rho/pnt': 'OTF0181.1',
    'rho': 'OTF0181.1',
    'pnt': 'OTF0181.1',
    'run/Bgb': 'FBgn0013753',
    'run': 'FBgn0013753',
    'hkbFFS': 'hkb_FBgn0001204',
    'cic': 'FBgn0028386',
    'retn': 'FBgn0004795',

}

rule get_all_memes:
    output: "prereqs/all_meme.meme"
    shell: """{module}; #module load wget
    wget -O prereqs/motif_databases.tgz \
         http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.17.tgz
    tar -xvf prereqs/motif_databases.tgz -C prereqs
    cat prereqs/motif_databases/FLY/* > {output}
    """
rule condense_memes:
    input:
        "prereqs/all_meme.meme",
    output:
        "Reference/all_meme_filtered.meme"
    shell: "python CondenseMemes.py"


rule alignment_figure:
    input:
            aln="analysis/targets/{gene}/{A}_{B}.needleall",
            fa="analysis/targets/{gene}/{A}_{B}.fasta",
            bed="analysis/targets/{gene}/{A}.bed",
            fimoA="analysis/targets/{gene}/{A}/fimo.txt",
            fimoB="analysis/targets/{gene}/{B}/fimo.txt",
            code="PatserAlignToSVG.py",
    output: "analysis/targets/{gene}/{A}_{B}.needleall.svg"
    log:    "analysis/targets/{gene}/{A}_{B}.needleall.log"
    wildcard_constraints:
        A="[a-z][a-z][a-z]+",
        B="[a-z][a-z][a-z]+",
    shell: """python {input.code} \
            -v --x-ticks 100 --draw-alignment \
            --x-scale 2 --y-sep 60 --y-scale 3.5 \
            --bar-width 10 \
            --comp1 {wildcards.A} --comp2 {wildcards.B} \
            --match-dim 0.7 \
            --meme-suite \
            --tf {tfs} --tf-names {tf_names} \
            --sequence --fasta {input.fa} \
            --needleall {input.aln} \
            --coordinates-bed {input.bed} \
            --bed-track Reference/bdtnp_dnase_2_prepend.bed \
            analysis/targets/{wildcards.gene}/
            """

rule alignment_figure_nondefaulttfs:
    input:
            aln="analysis/targets/{gene}/{A}_{B}.needleall",
            fa="analysis/targets/{gene}/{A}_{B}.fasta",
            bed="analysis/targets/{gene}/{A}.bed",
            fimoA="analysis/targets/{gene}/{A}/fimo.txt",
            fimoB="analysis/targets/{gene}/{B}/fimo.txt",
            code="PatserAlignToSVG.py",
    output: "analysis/targets/{gene}/{A}_{B}.{tfs}.needleall.svg"
    log:    "analysis/targets/{gene}/{A}_{B}.{tfs}.needleall.log"
    wildcard_constraints:
        A="[a-z][a-z][a-z]+",
        B="[a-z][a-z][a-z]+",
    run:
        local_tfnames = wildcards.tfs.split('_')
        local_tfs = [tf_dict[tf] for tf in local_tfnames]
        shell("""python {input.code} \\
                -v --x-ticks 100 --draw-alignment \\
                --x-scale 2 --y-sep 60 --y-scale 3.5 \\
                --bar-width 10 \\
                --comp1 {wildcards.A} --comp2 {wildcards.B} \\
                --match-dim 0.7 \\
                --meme-suite \\
                --tf {local_tfs} --tf-names {local_tfnames} \\
                --sequence --fasta {input.fa} \\
                --needleall {input.aln} \\
                --coordinates-bed {input.bed} \\
                --bed-track Reference/bdtnp_dnase_2_prepend.bed \\
                --outfile {output}\\
                analysis/targets/{wildcards.gene}/
                """)

rule alignment:
    input:
        A="analysis/targets/{gene}/{A}.fasta",
        B="analysis/targets/{gene}/{B}.fasta",
    wildcard_constraints:
        A="[a-z][a-z][a-z]+",
        B="[a-z][a-z][a-z]+",
    output:
        "analysis/targets/{gene}/{A}_{B}.needleall"
    log:
        "analysis/targets/{gene}/{A}_{B}.needleall.log"
    shell: """{module}; module load  EMBOSS
    needleall -aseq {input.A} -bseq {input.B} \
            -aformat3 srspair -gapopen 10.0 -gapextend 0.5 \
            -outfile {output}"""

rule combined_fasta:
    input:
        A="analysis/targets/{gene}/{A}.fasta",
        B="analysis/targets/{gene}/{B}.fasta",
    wildcard_constraints:
        A="[a-z][a-z][a-z]+",
        B="[a-z][a-z][a-z]+",
    output:
        "analysis/targets/{gene}/{A}_{B}.fasta"
    log:
        "analysis/targets/{gene}/{A}_{B}.fasta.log"
    shell:
        "cat {input.A} {input.B} > {output}"

rule fimo:
    input:
        fa="analysis/targets/{gene}/{species}.fasta",
        meme="Reference/all_meme_filtered.meme"
    output:
        "analysis/targets/{gene}/{species}/fimo.txt"
    log: "analysis/targets/{gene}/{species}/fimo.txt.log"
    shell:"""
    mkdir -p `dirname {output}`
    fimo -oc `dirname {output}` --thresh 1e-3 {input.meme} {input.fa} 2> {log}
    """

rule non_mel_bed:
    input:
        fa="analysis/targets/{gene}/mel.fasta",
        blastdb="Reference/d{species}_prepend.fasta.nhr"
    output: "analysis/targets/{gene}/{species}.bed"
    #wildcard_constraints: species="^(?!mel$).*$"
    shell: """{module}; module load blast bioawk
    blastn -db Reference/d{wildcards.species}_prepend.fasta \
            -outfmt "6 sseqid sstart send qseqid evalue sstrand length qlen slen qstart qend" \
            -gapextend 0\
            -query {input.fa} \
        | awk '!_[$4]++' \
        | bioawk -t '$10 > 1 && $6 ~ /minus/ {{$2 += $10 + 1}}; \
                    $10 > 1 && $6 ~ /plus/ {{$2 -= $10 + 1}}; \
                    $11 < $8 && $6 ~ /minus/ {{$3 -= ($8 - $11) + 1}}; \
                    $2 > $3 {{gsub("mel", "{wildcards.species}", $4); print $1,$3,$2+1,$4,$7/($8+1),"-", $3, $2 }}; \
                    $2 < $3 {{gsub("mel", "{wildcards.species}", $4); print $1,$2,$3+1,$4,$7/($8+1),"+", $2, $3 }}; '\
        > {output}
        """

rule mel_bed:
    input:
        gtf="Reference/mel_good.gtf",
        melsize="Reference/dmel.chr.sizes",
        oreganno="Reference/oreganno.prepend.bed",
        dnase="Reference/binding/dnase_peaks_prepend.bed",
    output:
        "analysis/targets/{gene}/mel.bed"
    shell: """{module}; module load bioawk bedtools
        mkdir -p `dirname output`
		grep '"{wildcards.gene}"' {input.gtf} \
		| bedtools sort \
		| bedtools merge  \
		| bedtools window -w 10000 -b - -a {input.dnase} \
        | bioawk -t '{{print $1, $2, $3, "mel_{wildcards.gene}_" NR, "0", "+", $2, $3}}' \
        | uniq --skip-fields 4 \
		> {output}

		grep '"{wildcards.gene}"' {input.gtf} \
		| bedtools sort \
		| bedtools merge  \
		| bedtools window -w 10000 -b - -a {input.dnase} \
        | bioawk -t '{{print $1, $2, $3, "mel_oreganno_{wildcards.gene}_" NR, "0", "+", $2, $3}}' \
        | uniq --skip-fields 4 \
        >> {output}
    """

rule nonmel_fasta_from_bed:
    input:
        bed='analysis/targets/{gene}/{species}.bed',
        full_fasta='Reference/d{species}_prepend.fasta',
    output:
        'analysis/targets/{gene}/{species}.fasta'
    shell: """ {module}; module load bedtools
    bedtools getfasta -fi {input.full_fasta} -bed {input.bed} \
            -fo {output} -s -name
    """
rule make_blastdb:
    input: "Reference/{file}.fasta"
    output: "Reference/{file}.fasta.nhr"
    shell: """{module}; module load blast; makeblastdb -dbtype nucl -in {input}"""

rule mel_chr_sizes:
    input: "Reference/dmel_prepend.fasta"
    output: "Reference/dmel.chr.sizes"
    run:
        out = open(output[0], 'w')
        for line in open(input[0]):
            if not line.startswith('>'):
                continue
            d=line.split()
            chrom = d[0][1:]
            length = next(el for el in d if el.startswith('length')).replace('length=', '').replace(';', '')
            print(chrom, length, sep='\t', file=out)


rule melfasta:
    input:
        bed="{prefix}/mel.bed",
        ref="Reference/dmel_prepend.fasta"
    output:
        "{prefix}/mel.fasta"
    shell: """{module}; module load bedtools
        bedtools getfasta -s -name -fi {input.ref} -fo {output} -bed {input.bed}
        """

rule getfasta:
    input:
        bed="{prefix}/{species}.bed",
        ref="Reference/d{species}_prepend.fasta"
    output:
        "{prefix}/{species}.fasta"
    wildcard_constraints: species="^...$"
    shell: """{module}; module load bedtools
        bedtools getfasta -s -name -fi {input.ref} -fo {output} -bed {input.bed}
        """

rule exons:
    input:
        bed="analysis/targets/{gene}/{species}.bed",
        gtf="Reference/{species}_good.gtf",
    output:
        "analysis/targets/{gene}/{species}_exons.gtf"
    shell: """{module}; module load bedtools
    bedtools intersect -b {input.bed} -a {input.gtf} -u \
        | grep exon \
        > {output}
        """

ruleorder: mel_bed > non_mel_bed

rule make_snpdir:
    input:
        vcf="analysis_godot/on_{target}/melsim_variants_on{target}.gvcf"
    output:
        dir="analysis_godot/on_{target}/snpdir",
        file="analysis_godot/on_{target}/snpdir/all.txt.gz",
    shell:"""
    mkdir -p {output.dir}
    cat {input.vcf}             \
            | grep -v "^#"             \
            | awk 'BEGIN {{OFS="\t"}}; length($4) == 1 && length($5) == 1 {{print $1,$2,$4,$5}};' \
            | gzip -c  \
            > {output.file}
"""

rule wasp_find_snps:
    input:
        bam="{sample}/{prefix}_dedup.bam",
        bai="{sample}/{prefix}_dedup.bam.bai",
        snpdir="analysis_godot/on_mel/snpdir",
        snpfile="analysis_godot/on_mel/snpdir/all.txt.gz"
    output:
        temp("{sample}/{prefix}_dedup.remap.fq1.gz"),
        temp("{sample}/{prefix}_dedup.remap.fq2.gz"),
        temp("{sample}/{prefix}_dedup.keep.bam"),
        temp("{sample}/{prefix}_dedup.to.remap.bam"),

    shell:
        """python ~/FWASP/mapping/find_intersecting_snps.py \
            --progressbar \
            --phased --paired_end \
            {input.bam} {input.snpdir}
        """


rule wasp_remap:
    input:
        R1="{sample}/{prefix}.remap.fq1.gz",
        R2="{sample}/{prefix}.remap.fq2.gz",
        genome="Reference/dmel_prepend/Genome",
        genomedir="Reference/dmel_prepend/"
    output:
        temp("{sample}/{prefix}.remap.bam")
    threads: 16
    shell: """{module}; module load STAR;
    rm -rf {wildcards.sample}/STARtmp
    STAR \
            --genomeDir {input.genomedir} \
            --outFileNamePrefix {wildcards.sample}/remap \
            --outSAMattributes MD NH --clip5pNbases 6 \
            --outSAMtype BAM Unsorted \
            --outTmpDir {wildcards.sample}/STARtmp \
            --limitBAMsortRAM 20000000000 \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn {input.R1} {input.R2}
    mv {wildcards.sample}/remapAligned.out.bam {output}
            """

rule wasp_keep:
    input:
        toremap="{file}.to.remap.bam",
        remapped="{file}.remap.bam",
    output:
        temp("{file}.remap.kept.bam"),
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python ~/FWASP/mapping/filter_remapped_reads.py \
            -p \
            {input.toremap} {input.remapped} \
            {output} """

rule wasp_merge:
    input:
        "{file}.remap.kept.bam",
        "{file}.keep.bam",
    output:
        temp("{file}.keep.merged.bam")
    shell:
        "{module}; module load samtools; samtools merge {output} {input}"

def samples_to_files(fname):
    # Returns a function that is callable by an input rule
    def retfun():
        return [path.join(analysis_dir, sample, fname) for sample in samples]
    return retfun

rule all_files_per_sample:
    output: touch("tmp/all_{fname}")
    input: lambda wildcards: samples_to_files(wildcards.fname)()

rule sample_expr:
    input:
        bam="{sample}/assigned_dmelR.bam",
        bai="{sample}/assigned_dmelR.bam.bai",
        gtf=mel_gtf,
        fasta=mel_fasta,
        mask_gtf=mel_bad_gtf,
        sentinel=path.join(analysis_dir, 'recufflinks')
    threads: 4
    output:
        "{sample}/genes.fpkm_tracking"
    log:
        "{sample}/genes.log"
    shell: """{module}; module load cufflinks
	cufflinks \
		--num-threads 8 \
		--output-dir {wildcards.sample}/ \
		--multi-read-correct \
		--frag-bias-correct {input.fasta} \
		--GTF {input.gtf} \
		--mask-file {input.mask_gtf} \
		{input.bam}
        """


rule sample_gene_ase:
    input:
        bam="{sample}/assigned_dmelR_dedup.bam",
        bai="{sample}/assigned_dmelR_dedup.bam.bai",
        variants=variants,
        hets=path.join(analysis_dir, "on_mel", "melsim_true_hets.tsv"),
        gtf=mel_gtf,
        sentinel=path.join(analysis_dir, 'recalc_ase')
    threads: 1
    output:
        "{sample}/melsim_gene_ase_by_read.tsv"
    log:
        "{sample}/melsim_gene_ase_by_read.log"
    shell: """
    source activate peter
    export PYTHONPATH=$HOME/ASEr/;
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --id-name gene_name \
        --ase-function pref_index \
        --min-reads-per-allele 0 \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """

rule wasp_gene_ase:
    input:
        bam="{sample}/orig_dedup.keep.merged.sorted.bam",
        bai="{sample}/orig_dedup.keep.merged.sorted.bam.bai",
        variants=variants,
        hets=path.join(analysis_dir, "on_mel", "melsim_true_hets.tsv"),
        gtf=mel_gtf,
        sentinel=path.join(analysis_dir, 'recalc_ase')
    threads: 1
    output:
        "{sample}/wasp_gene_ase_by_read.tsv"
    shell: """
    source activate peter
    export PYTHONPATH=$HOME/ASEr/;
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --id-name gene_name \
        --assign-all-reads \
        --ase-function pref_index \
        --min-reads-per-allele 0 \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """

rule sample_cds_ase:
    input:
        bam="{sample}/assigned_dmelR_dedup.bam",
        bai="{sample}/assigned_dmelR_dedup.bam.bai",
        variants=variants,
        hets=path.join(analysis_dir, "on_mel", "melsim_true_hets.tsv"),
        gtf=mel_gtf,
        sentinel=path.join(analysis_dir, 'recalc_ase')
    threads: 1
    output:
        "{sample}/melsim_cds_ase_by_read.tsv"
    log:
        "{sample}/melsim_cds_ase_by_read.log"
    shell: """
    source activate peter
    export PYTHONPATH=$HOME/ASEr/;
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --id-name gene_name \
        --ase-function pref_index \
        --min-reads-per-allele 0 \
        --feature-type CDS \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """

rule sample_exon_ase:
    input:
        bam="{sample}/assigned_dmelR_dedup.bam",
        bai="{sample}/assigned_dmelR_dedup.bam.bai",
        variants=variants,
        hets=path.join(analysis_dir, "on_mel", "melsim_true_hets.tsv"),
        gtf='Reference/mel_renamed_exons.gtf',
        sentinel=path.join(analysis_dir, 'recalc_ase')
    threads: 1
    output:
        "{sample}/melsim_exon_ase_by_read.tsv"
    log:
        "{sample}/melsim_exon_ase_by_read.log"
    shell: """
    source activate peter
    export PYTHONPATH=$HOME/ASEr/;
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --id-name gene_id \
        --ase-function pref_index \
        --min-reads-per-allele 0 \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """


rule exons_gtf:
    input:
        "Reference/{species}_good.gtf"
    output:
        "Reference/{species}_good_exons.gtf"
    shell:"""
    python2 $HOME/R/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
        --aggregate=yes {input} {output}
    """

rule exons_renamed:
    input:
        "Reference/{species}_good_exons.gtf"
    output:
        "Reference/{species}_renamed_exons.gtf"
    run:
        import Utils as ut
        with open(output[0], 'w') as out:
            for line in open(input[0]):
                data = line.split('\t')
                if data[2] != 'exonic_part': continue
                annots = ut.parse_annotation(data[-1])
                data[-1] = 'gene_id "{}_{}"'.format(annots['gene_id'], annots['exonic_part_number'])
                data[2] = 'exon'
                print(*data, file=out, sep='\t')




rule sample_psi:
    input:
        bam="{sample}/assigned_dmelR_dedup.bam",
        bai="{sample}/assigned_dmelR_dedup.bam.bai",
        gtf="Reference/mel_good_exons.gtf",
        sentinel=path.join(analysis_dir, 'recalc_psi')
    output:
        "{sample}/psi.tsv"
    shell:""" {module}; module load fraserconda
    python CalculatePSI.py \
        --outfile {output} \
        {input.bam} {input.gtf}
        """

rule sample_juncs:
    input:
        bam="{analysis_dir}/{sample}/assigned_dmelR_dedup.bam",
        bai="{analysis_dir}/{sample}/assigned_dmelR_dedup.bam.bai",
    output:
        "{analysis_dir}/velvetant/{sample}.junc"
    shell:"""
    ~/leafcutter/scripts/bam2junc.sh {input.bam} {output}
    """

rule all_sample_juncs:
    input: lambda wildcards: [path.join(analysis_dir, 'velvetant', sample +'.junc') for sample in samples]
    output: 'analysis/velvetant/juncfiles.txt'
    run:
        with open(output[0], 'w') as out:
            print(*input, sep='\n', file=out)


rule leafcutter_cluster:
    input:
        juncs='analysis/velvetant/juncfiles.txt',
        firstbam=samples_to_files('assigned_dmelR_dedup.bam')()[0],
    output:
        'analysis/velvetant/clusters_perind.counts.gz'
    shell: """
    {module}; module load fraserconda
    python ~/leafcutter-official/clustering/leafcutter_cluster.py \
        -j {input.juncs}  \
        -o clusters \
        --example-bam {input.firstbam} \
        --rundir analysis/velvetant/
    """

rule sample_velvetant:
    input:
        bam="{sample}/assigned_dmelR_dedup.bam",
        bai="{sample}/assigned_dmelR_dedup.bam.bai",
        juncs="analysis/velvetant/clusters_perind.counts.gz",
        snps="analysis_godot/on_mel/melsim_variant.bed",
    output:
        "{sample}/velvetant.tsv"
    shell: """
    {module}; module load fraserconda
    python VelvetAnt.py -x --snps-bed {input.snps} --splicing-clusters {input.juncs} -o {output} {input.bam}
    """

rule get_all_map_stats:
    input: *samples_to_files('assigned_dmelR.mapstats')(),
    output:
        path.join(analysis_dir, 'map_stats.tsv'),
    log:
        path.join(analysis_dir, 'map_stats.log'),
    shell:"""
    {module}; module load fraserconda;
	python GetMapStats.py \
		--params Parameters/RunConfig.cfg \
		--count-unique \
		--count-all \
		--translate-labels \
		{analysis_dir}
        """

rule get_sample_mapstats:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        bam="{sample}/{fname}.bam",
        bai="{sample}/{fname}.bam.bai",
    output: "{sample}/{fname}.mapstats"
    log: "{sample}/mapstats.log"
    shell: "{module}; module load fraserconda; python GetSingleMapStats.py {input.bam}"

rule snp_counts:
    input:
        bam="{sample}/assigned_dmelR_dedup.bam",
        variants=path.join(analysis_dir,  "on_{parent}", "{parent}{other}_variant.bed"),
    output:
        "{sample}/{parent}{other}_SNP_COUNTS.txt"
    wildcard_constraints:
        parent='[a-z][a-z][a-z]',
        other='[a-z][a-z][a-z]',
    shell:"""
        {module}; module load samtools
        mkdir -p {wildcards.sample}/melsim_countsnpase_tmp
        python2 CountSNPASE.py \
            --mode single \
            --reads {input.bam} \
            --snps {input.variants} \
            --prefix {wildcards.sample}/melsim_countsnpase_tmp/
        mv {wildcards.sample}/melsim_countsnpase_tmp/_SNP_COUNTS.txt {output}
        rm -rf {wildcards.sample}/melsim_countsnpase_tmp

    """

rule true_hets:
    input:
        *samples_to_files('{parent}{other}_SNP_COUNTS.txt')(),
    output:
        path.join(analysis_dir, "on_{parent}", "{parent}{other}_true_hets.tsv")
    shell: """
    {module}; module load fraserconda
    python GetTrueHets.py \
        --min-counts 10 \
        --outfile {output} \
        {input}
    cp {output} `dirname {output}`/true_hets.tsv
    """

rule kallisto_summary:
    input:
        *samples_to_files('abundance.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'summary_kallisto.tsv')
    log:
        path.join(analysis_dir, 'mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename abundance.tsv \
	   --key target_id \
       --out-basename summary_kallisto \
	   --column tpm \
		{analysis_dir} \
		| tee {log}
        """

rule expr_summary:
    input:
        *samples_to_files('genes.fpkm_tracking')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'summary.tsv')
    log:
        path.join(analysis_dir, 'mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename genes.fpkm_tracking \
	   --key gene_short_name \
	   --column FPKM \
		{analysis_dir} \
		| tee {log}
        """

rule ase_summary:
    input:
        *samples_to_files('melsim_gene_ase_by_read.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'ase_summary_by_read.tsv')
    log:
        path.join(analysis_dir, 'ase_mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename melsim_gene_ase_by_read.tsv \
	   --key gene \
	   --column ase_value \
       --out-basename ase_summary_by_read \
		{analysis_dir} \
		| tee {log}
        """

rule ase_summary_refalt:
    input:
        *samples_to_files('wasp_gene_ase_by_read.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'ase_summary_refalt.tsv')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename wasp_gene_ase_by_read.tsv \
	   --key gene \
       --refalt \
       --float-format %3.0f \
       --exclude-column chrom \
       --exclude-column-value dmel_X \
       --exclude-samples melXsim_cyc14C_rep3 simXmel_cyc14C_rep2 \
	   --column ref_counts \
       --out-basename ase_summary_refalt \
		{analysis_dir}
        """

rule ase_summary_wasp:
    input:
        *samples_to_files('wasp_gene_ase_by_read.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'wasp_summary_by_read.tsv')
    log:
        path.join(analysis_dir, 'ase_mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename wasp_gene_ase_by_read.tsv \
	   --key gene \
	   --column ase_value \
       --out-basename wasp_summary_by_read \
		{analysis_dir} \
		| tee {log}
        """

rule cds_ase_summary:
    input:
        *samples_to_files('melsim_cds_ase_by_read.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'cds_ase_summary_by_read.tsv')
    log:
        path.join(analysis_dir, 'cds_mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename melsim_cds_ase_by_read.tsv \
	   --key gene \
	   --column ase_value \
       --out-basename cds_ase_summary_by_read \
		{analysis_dir} \
		| tee {log}
        """

rule exons_ase_summary:
    input:
        *samples_to_files('melsim_exon_ase_by_read.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'exon_ase_summary_by_read.tsv')
    log:
        path.join(analysis_dir, 'exon_mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename melsim_exon_ase_by_read.tsv \
	   --key gene \
	   --column ase_value \
       --out-basename exon_ase_summary_by_read \
		{analysis_dir} \
		| tee {log}
        """

rule psi_summary:
    input:
        *samples_to_files('psi.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'psi_summary.tsv')
    log:
        path.join(analysis_dir, 'psi_mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename psi.tsv \
	   --key exon_id \
	   --column psi \
       --out-basename psi_summary \
		{analysis_dir} \
		| tee {log}
        """

rule velvetant_summary:
    input:
        *samples_to_files('velvetant.tsv')(),
        map_stats=path.join(analysis_dir, 'map_stats.tsv'),
        sentinel=path.join(analysis_dir, 'retabulate'),
    output:
        path.join(analysis_dir, 'velvetant_summary.tsv')
    log:
        path.join(analysis_dir, 'velvetant_mst.log')
    shell: """
    {module}; module load fraserconda;
    python MakeSummaryTable.py \
	   --strip-low-reads 1000000 \
	   --strip-on-unique \
	   --strip-as-nan \
	   --mapped-bamfile assigned_dmelR.bam \
	   --strip-low-map-rate 52 \
	   --map-stats {input.map_stats} \
	   --filename velvetant.tsv \
	   --key 0 \
	   --column pref_index \
       --out-basename velvetant_summary \
		{analysis_dir} \
		| tee {log}
        """


rule sort_bam:
    input: "{sample}.bam"
    output: "{sample}.sorted.bam"
    log: "{sample}.sorted.log"
    threads: 4
    shell: """{module}; module load samtools
    samtools sort -o {output} --threads {threads} {input}
    """

rule index_bam:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    log: "{sample}.bam.bai_log"
    shell: "{module}; module load samtools; samtools index {input}"

rule dedup:
    input: "{sample}.bam"
    output: ("{sample}_dedup.bam")
    log: "{sample}_dedup.log"
    shell: """{module}; module load picard
    picard MarkDuplicates \
        SORTING_COLLECTION_SIZE_RATIO=.05 \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        MAX_RECORDS_IN_RAM=2500000 \
		READ_NAME_REGEX=null \
		REMOVE_DUPLICATES=true \
		DUPLICATE_SCORING_STRATEGY=RANDOM \
		INPUT={input} OUTPUT={output} METRICS_FILE={log}
        """

rule get_sra:
    output:
        "sequence/{srr}_1.fastq.gz",
        "sequence/{srr}_2.fastq.gz"
    #log: "sequence/{srr}.log"
    resources: max_downloads=1
    shell: """{module}; module load sra-tools

    fastq-dump --gzip --split-3 --outdir sequence {wildcards.srr}
    """

rule star_map:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        genome=path.join(hybrid_dir, 'masked', 'Genome'),
        genomedir=path.join(hybrid_dir, 'masked')
    params:
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
    output: "{sample}/assigned_dmelR.bam"
    threads: 6
    log: "{sample}/assigned_dmelR.log"
    shell: """{module}; module load STAR
    STAR --parametersFiles Parameters/STAR_params.in \
    --genomeDir {input.genomedir} \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix {wildcards.sample}/ \
    --outSAMattributes MD NH \
    --outSAMtype BAM SortedByCoordinate \
    --clip5pNbases 6 \
    --readFilesIn {params.r1s} {params.r2s}
    if [ -s  {wildcards.sample}/Aligned.sortedByCoord.out.bam ]; then
        mv {wildcards.sample}/Aligned.sortedByCoord.out.bam {output};
    fi
    """

rule star_map_mel:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        genome='Reference/dmel_prepend/Genome',
        genomedir='Reference/dmel_prepend',
    params:
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
    output: "{sample}/orig.bam"
    threads: 6
    shell: """{module}; module load STAR
    STAR --parametersFiles Parameters/STAR_params.in \
    --genomeDir {input.genomedir} \
    --outFileNamePrefix {wildcards.sample}/ \
    --outSAMattributes MD NH \
    --outSAMtype BAM SortedByCoordinate \
    --clip5pNbases 6 \
    --readFilesIn {params.r1s} {params.r2s}
    if [ -s  {wildcards.sample}/Aligned.sortedByCoord.out.bam ]; then
        mv {wildcards.sample}/Aligned.sortedByCoord.out.bam {output};
    fi
    """

def interleave_reads(wildcards):
    retval =  " ".join(" ".join([a, b])
            for a, b
            in zip(
                getreads_nowc(1)(path.basename(wildcards.sample)),
                getreads_nowc(2)(path.basename(wildcards.sample))
                )
            )
    return retval

rule makedir:
    output: "{prefix}.log", "{prefix}/"
    shell: "touch {wildcards.prefix}.log; mkdir -p {wildcards.prefix}"

rule kallisto_quant:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        index='Reference/dmel_5.57_kallisto',
        dir=ancient('{sample}/')
    priority: 50
    params:
        reads=interleave_reads
    output: "{sample}/abundance.h5", "{sample}/abundance.tsv"
    shell:"""
    mkdir -p {wildcards.sample}
    ~/Downloads/kallisto/kallisto quant \
        -i {input.index} \
        -o {wildcards.sample} \
        {params.reads}
    """


rule mel_star_ref:
    input:
        fasta='Reference/dmel_prepend.fasta',
        gtf=mel_gtf,
    output:
        outdir='Reference/dmel_prepend',
        outfile="Reference/dmel_prepend/Genome"
    log: path.join(hybrid_dir, "masked.log")
    priority: 50
    shell:""" {module}; module load STAR
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
	STAR --runMode genomeGenerate --genomeDir {output.outdir} \
		--outTmpDir {output.outdir}/_tmp/ \
		--genomeFastaFiles {input.fasta} \
		--sjdbGTFfile {input.gtf}
    """

rule masked_star_ref:
    input:
        fasta=path.join(hybrid_dir, "melsim_masked.fasta"),
        gtf=mel_gtf,
    output:
        outdir=path.join(hybrid_dir, "masked"),
        outfile=path.join(hybrid_dir, "masked/Genome")
    log: path.join(hybrid_dir, "masked.log")
    priority: 50
    shell:""" {module}; module load STAR
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
	STAR --runMode genomeGenerate --genomeDir {output.outdir} \
		--outTmpDir {output.outdir}/_tmp/ \
		--genomeFastaFiles {input.fasta} \
		--sjdbGTFfile {input.gtf}
    """

rule ref_genome:
    output: "prereqs/d{species}.fasta"
    log: "prereqs/d{species}.log"
    params:
        date=lambda wildcards: dates[wildcards.species],
        version=lambda wildcards: versions[wildcards.species]

    shell: """{module}; #module load wget
    mkdir -p prereqs
	wget -O {output}.gz ftp://ftp.flybase.org/releases/{params.date}/d{wildcards.species}_{params.version}/fasta/d{wildcards.species}-all-chromosome-{params.version}.fasta.gz
	gunzip --force {output}.gz
    """

rule ref_gff:
    output: "prereqs/d{species}.gff"
    log: "prereqs/d{species}.log"
    params:
        date=lambda wildcards: dates[wildcards.species],
        version=lambda wildcards: versions[wildcards.species]

    shell: """{module}; #module load wget
    mkdir -p prereqs
	wget -O {output}.gz ftp://ftp.flybase.org/releases/{params.date}/d{wildcards.species}_{params.version}/gff/d{wildcards.species}-all-{params.version}.gff.gz
	gunzip --force {output}.gz
    """

rule prepend_fasta:
    output: "Reference/d{species}_prepend.fasta"
    log: "Reference/d{species}_prepend.log"
    input:
        "prereqs/d{species}.fasta"
    shell: "perl -pe 's/>/>d{wildcards.species}_/' {input} > {output}"

rule all_gtf:
    input: "prereqs/d{species}.gff"
    output: "Reference/{species}_all.gtf"
    log: "Reference/{species}_all.log"
    shell: """{module}; module load cufflinks
    gffread {input} -C -E -T -o- \
        | awk '{{print "d{wildcards.species}_"$0}}' \
        > {output}
    """

rule good_gtf:
    input: "Reference/{species}_all.gtf"
    output: "Reference/{species}_good.gtf"
    log: "Reference/{species}_good.log"
    shell: """
	cat {input} \
		| grep -vP '(snoRNA|CR[0-9]{{4}}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		| grep 'gene_id' \
		> {output}
        """

rule bad_gtf:
    input: "Reference/{species}_all.gtf"
    output: "Reference/{species}_bad.gtf"
    log: "Reference/{species}_bad.log"
    shell: """
    echo {threads}
	cat {input} \
		| grep -P '(snoRNA|CR[0-9]{{4}}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> {output}
        """

rule process_variants:
    input:
        var_tab=path.join(analysis_dir, "on_{parent}", "{parent}{other}_variants.tsv"),
        ref_fasta="Reference/d{parent}_prepend.fasta",
    output:
        fasta=path.join(analysis_dir, "on_{parent}", "{parent}{other}_masked.fasta"),
        bed=path.join(analysis_dir, "on_{parent}", "{parent}{other}_variant.bed"),
    log: path.join(analysis_dir, "on_{parent}", "{parent}{other}_masking.log")
    shell: """
	python MaskReferenceFromGATKTable.py \
		--target-species {wildcards.parent} \
		--emit-bed {output.bed} \
		--outfasta {output.fasta} \
		{input.ref_fasta} \
		{input.var_tab}
    """


rule get_combined_variants:
    input:
        ref_fasta="Reference/d{parent}_prepend.fasta",
        parent_gvcf=path.join(analysis_dir, "on_{parent}", "{parent}_gdna_raw_variants_uncalibrated.p.g.vcf"),
        other_gvcf=path.join(analysis_dir, "on_{parent}", "{other}_gdna_raw_variants_uncalibrated.p.g.vcf")
    output: path.join(analysis_dir, "on_{parent}", "{parent}{other}_variants.tsv")
    log:    path.join(analysis_dir, "on_{parent}", "{parent}{other}_variants.log")
    params:
        dir=path.join(analysis_dir, "on_{parent}")
    shell: """
    {module}; module load java
	gatk -T GenotypeGVCFs \
		-R {input.ref_fasta} \
		-V {input.parent_gvcf} \
		-V {input.other_gvcf} \
		-o {params.dir}/{wildcards.parent}{wildcards.other}_variants_on{wildcards.parent}.gvcf
	gatk -T VariantsToTable \
		-R {input.ref_fasta} \
		-V {params.dir}/{wildcards.parent}{wildcards.other}_variants_on{wildcards.parent}.gvcf \
		-F CHROM -F POS -F REF -F ALT -F QUAL \
		-F HET -F HOM-REF -F HOM-VAR -F NCALLED \
		-GF GT \
		-o {output}


    """

rule fadict:
    input: "{file}.fasta"
    output: "{file}.dict"
    shell: "{module}; module load picard; picard CreateSequenceDictionary R={input} O={output}"

rule index_fasta:
    input: "{file}.fasta"
    output: "{file}.fasta.fai"
    shell: "samtools faidx {input}"

rule call_variants:
    input:
        ref_fasta="Reference/d{reference}_prepend.fasta",
        ref_fai="Reference/d{reference}_prepend.fasta.fai",
        ref_dict="Reference/d{reference}_prepend.dict",
        bam=path.join(analysis_dir, "on_{reference}", "{species}_gdna_bowtie2_dedup.bam"),
        bai=path.join(analysis_dir, "on_{reference}", "{species}_gdna_bowtie2_dedup.bam.bai"),
    output:
        path.join(analysis_dir, "on_{reference}", "{species}_gdna_raw_variants_uncalibrated.p.g.vcf")
    log:
        path.join(analysis_dir, "on_{reference}", "{species}_gdna_raw_variants_uncalibrated.log")
    threads: 4
    shell: """
    {module}; module load java
	gatk -T HaplotypeCaller \
		-R {input.ref_fasta} \
		-I {input.bam} \
		-nct 16 \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o {output}
    """

rule split_genome_to_regions:
    input: "Reference/d{reference}.chr.sizes"
    output: expand("Reference/d{reference}.split/10m_{subset}.bed", reference="{reference}", subset=range(num_mel_windows))
    log: dynamic("Reference/d{reference}.split/split.log")
    shell: "python SplitGenome.py \
        --size 10e6 \
        {input} \
        Reference/d{wildcards.reference}.split/10m"

rule call_all_region_variants:
    input:
        ref_fasta="Reference/d{reference}_prepend.fasta",
        variants = expand(path.join(analysis_dir, "on_{{reference}}",
                    "{{species}}_gdna_variants", "{subset}_raw_variants_uncalibrated.p.g.vcf"),
                    subset=range(num_mel_windows)),
    output: path.join(analysis_dir, "on_{reference}", "{species}_gdna_raw_variants_combined.p.g.vcf")
    log: path.join(analysis_dir, "on_{reference}", "{species}_gdna_raw_variants_combined.log")
    run:
        variant_files = " -V ".join(input.variants)
        shell("""java -cp /usr/share/java/gatk/GenomeAnalysisTK.jar \
                org.broadinstitute.gatk.tools.CatVariants \
                -R {input.ref_fasta} --outputFile {output} -assumeSorted -V """ + variant_files)

rule call_variants_region:
    input:
        ref_fasta="Reference/d{reference}_prepend.fasta",
        ref_fai="Reference/d{reference}_prepend.fasta.fai",
        ref_dict="Reference/d{reference}_prepend.dict",
        region="Reference/d{reference}.split/10m_{subset}.bed",
        bam=path.join(analysis_dir, "on_{reference}", "{species}_gdna_bowtie2_dedup.bam"),
        bai=path.join(analysis_dir, "on_{reference}", "{species}_gdna_bowtie2_dedup.bam.bai"),
    output:
        path.join(analysis_dir, "on_{reference}", "{species}_gdna_variants/{subset}_raw_variants_uncalibrated.p.g.vcf")
    params:
        outdir=path.join(analysis_dir, "on_{reference}", "{species}_gdna/")
    log:
        path.join(analysis_dir, "on_{reference}", "{species}_gdna/{subset}_raw_variants_uncalibrated.log")
    threads: 6
    shell: """
    mkdir -p {params.outdir}
	gatk -T HaplotypeCaller \
		-R {input.ref_fasta} \
		-I {input.bam} \
		-nct 8 \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
        --intervals {input.region} \
		-o {output}
    """

rule map_gdna:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        ancient(path.join(analysis_dir, "on_{reference}")+'/'),
        bt2_index="Reference/d{reference}_prepend.1.bt2"
    output:
        path.join(analysis_dir, "on_{reference}", "{sample}_bowtie2.bam")
    log:
        path.join(analysis_dir, "on_{reference}", "{sample}_bowtie2.log")
    params:
        index="Reference/d{reference}_prepend",
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
        outdir= lambda wildcards, output: path.dirname(output[0])
    threads: 12
    shell: """{module}; module load samtools/1.3 bowtie2
    bowtie2 \
		--very-sensitive-local \
		-p 11 \
		--rg-id {wildcards.sample} \
		--rg "SM:{wildcards.sample}" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x {params.index} \
		-1 {params.r1s} \
		-2 {params.r2s} \
		| samtools view -b \
		| samtools sort -o {output} -T {params.outdir}/{wildcards.sample}_bowtie2_sorting
        """

rule bowtie2_build:
    input: "{base}.fasta"
    output: "{base}.1.bt2"
    log:    "{base}.bt2.log"
    shell: "{module}; module load bowtie2; bowtie2-build --offrate 3 {input} {wildcards.base}"

rule sentinel_cufflinks:
    output: touch(path.join(analysis_dir, "recufflinks"))

rule sentinel_retabulate:
    output: touch(path.join(analysis_dir, "retabulate"))

rule sentinel_recalc_ase:
    output: touch(path.join(analysis_dir, "recalc_ase"))

rule sentinel_recalc_psi:
    output: touch(path.join(analysis_dir, "recalc_psi"))

ruleorder: ref_genome > melfasta > getfasta > combined_fasta > nonmel_fasta_from_bed

rule get_ase_stats:
    input:
       "analysis_godot/wasp_summary_by_read.tsv",
       "analysis_godot/summary.tsv",
       "prereqs/journal.pbio.1000590.s002",
       "prereqs/gene_map_table_fb_2016_01.tsv",
       "analysis/mel_deseq.tsv",
       "analysis/sim_deseq.tsv",
       'analysis/results/asepeak_genes.txt',
       'analysis/results/aselogist_genes.txt',
       'analysis/results/fd_peak.numpy',
       'analysis/results/fd_logist.numpy',
       'analysis/results/all_peak_r2s.csv',
       'analysis/results/all_logist_r2s.csv',
    output:
        expand('analysis/results/{fname}',
                fname=['maternal.txt', 'paternal.txt', 'mel_dom.txt', 'sim_dom.txt',
                        'mel_dom.svg', 'sim_dom.svg',
                        'me_zyg_lott_mat.svg',
                        'me_mat_lott_zyg.svg',
                        'stats.tex',
                ])
    shell: """ {module}
    module load fraserconda
    python GetASEStats.py"""

rule ase_figure:
    input:
        ase="analysis_godot/ase_summary_by_read.tsv",
        expression="analysis_godot/summary.tsv",
    output:
        expand("analysis/results/ase{fit}_{data}",
                fit=['logist', 'peak'],
                data=["genes.txt", "ase.svg", "ase.png",
                "ase_r1.svg", "ase_r1.png", "expr_hyb.svg",
                "expr.svg", "fits.png"],
                )
    shell: """{module}
    module load fraserconda
    python FitASEFuncs.py \
        --expression-file {input.expression} \
        --ase-file {input.ase} \
        --prefix ase \
        --male-samples melXsim_cyc14C_rep3 simXmel_cyc14C_rep2 \
        --cutoff-r2 0.45
        """

