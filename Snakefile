from os import path

analysis_dir = 'analysis_godot'
hybrid_dir = 'analysis_godot/on_melr5'
mel_gtf = 'Reference/mel_good.gtf'
mel_bad_gtf = 'Reference/mel_bad.gtf'
mel_fasta = 'Reference/dmel_prepend.fasta'

variants = path.join(hybrid_dir, 'melsim_variant.bed')

configfile: "Parameters/config.json"

mxs_r1 = ['RA{:02d}'.format(i) for i in range(4, 30)]
mxs_r2 = ['RD{:02d}'.format(i) for i in range(7, 33)] 
mxs_r3 = ['RE{:02d}'.format(i) for i in range(5, 32)] 
sxm_r1 = ['MA{:02d}'.format(i) for i in range(37, 62)] 
sxm_r2 = ['MB{:02d}'.format(i) for i in range(3, 30)] 
mxm_r1 = ['MC{:02d}'.format(i) for i in range(7, 34)] 
sxs_r1 = ['SB{:02d}'.format(i) for i in range(7, 34)] 

samples = mxs_r1 + mxs_r2 + mxs_r3 + sxm_r1 + sxm_r2 + mxm_r1 + sxs_r1


module = '''module () {
        eval `$LMOD_CMD bash "$@"`
        }'''

rule all:
    input:
        path.join(analysis_dir, "summary.tsv"),
        path.join(analysis_dir, "ase_summary_by_read.tsv")


rule print_config:
    run:
        print(config)

tfs = "OTF0070.1 hb  kni gt D FBgn0001325_4 tll hkb  FBgn0000251_3 FBgn0003448_3 twi OTF0532.1 OTF0478.1 OTF0181.1 FBgn0013753 hkb_FFS"
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
    'hkbFFS': 'hkb_FFS',
    'cic': 'FBgn0028386',
    'retn': 'FBgn0004795',

}

rule alignment_figure:
    input:
            aln="analysis/targets/{gene}/{A}_{B}.needleall",
            fa="analysis/targets/{gene}/{A}_{B}.fasta",
            bed="analysis/targets/{gene}/{A}.bed",
            fimoA="analysis/targets/{gene}/{A}/fimo.txt",
            fimoB="analysis/targets/{gene}/{B}/fimo.txt",
            code="PatserAlignToSVG.py",
    output: "analysis/targets/{gene}/{A}_{B}.needleall.svg"
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
        ref="Reference/dmelr5_prepend.fasta"
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
ruleorder: melfasta > getfasta

def samples_to_files(fname):
    # Returns a function that is callable by an input rule
    def retfun():
        return [path.join(analysis_dir, sample, fname) for sample in samples]
    return retfun

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
        bam="{sample}/assigned_dmelR_dedup.sorted.bam",
        bai="{sample}/assigned_dmelR_dedup.sorted.bam.bai",
        variants=variants,
        gtf=mel_gtf,
        sentinel=path.join(analysis_dir, 'recalc_ase')
    threads: 1
    output:
        "{sample}/melsim_gene_ase_by_read.tsv"
    shell: """ python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --id-name gene_name \
        --ase-function pref_index \
        {input.variants} \
        {input.gtf} \
        {input.bam}
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
    python MakeSummaryTable.py \
       --params Parameters/RunConfig.cfg \
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
    python MakeSummaryTable.py \
       --params Parameters/RunConfig.cfg \
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
        

rule sort_bam:
    input: 
        "{sample}.bam"
    output:
        "{sample}.sorted.bam"
    threads: 4
    shell: """{module}; module load samtools
    samtools sort -o {output} --threads {threads} {input}
    """

rule index_bam:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    shell: "samtools index {input}"

rule dedup:
    input: "{sample}.bam"
    output: temp("{sample}_dedup.bam")
    log: "{sample}_dedup.log"
    shell: """{module}; module load picard/2.8.1
    picard MarkDuplicates
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		READ_NAME_REGEX=null \
		INPUT={input} OUTPUT={output} METRICS_FILE={log}
        """
rule get_sra:
    output:
        "sequence/{srr}_1.fastq.gz",
        "sequence/{srr}_2.fastq.gz"
    resources: max_downloads=1
    shell: """{module}; module load sra-tools

    fastq-dump --gzip --split-3 --outdir sequence {srr}
    """

def getreads(readnum):
    def retfun(wildcards):
        return config['samples'][path.basename(wildcards.sample)]
    return retfun

rule star_map:
    input: 
        r1s=getreads(1),
        r2s=getreads(2),
        genome=path.join(hybrid_dir, 'masked', 'Genome'),
        genomedir=path.join(hybrid_dir, 'masked')
    output: "{sample}/assigned_dmelR.bam"
    shell: """{module}; module load STAR
    STAR --parametersFiles Parameters/STAR_params.in \
    --genomeDir {input.genome} \
    --outFileNamePrefix {wildcards.sample}/ \
    --outSAMattributes MD NH \
    --outSAMtype BAM SortedByCoordinate \
    --clip5pNbases 6 \
    --readFilesIn {input.r1s} {input.r2s}
    ln -s {wildcards.sample}/Aligned.sortedByCoord.out.bam {output}
    """


rule masked_star_ref:
    input: 
        fasta=path.join(hybrid_dir, "melsim_masked.fasta"),
        gtf=mel_gtf,
    output:
        outdir=path.join(hybrid_dir, "masked"),
        outfile=path.join(hybrid_dir, "masked/Genome")
    shell:""" {module}; module load STAR
    rm -rf {output.outdir}
	STAR --runMode genomeGenerate --genomeDir {output.outdir} \
		--outTmpDir {output.outdir}/_tmp/ \
		--genomeFastaFiles {input.fasta} \
		--sjdbGTFfile {input.gtf}
    """
