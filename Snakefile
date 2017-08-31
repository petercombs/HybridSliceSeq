#samples = (['RA{:02}'.format(i) for i in range(04, 30)]   # mXs r1
#           + ['RD{:02}'.format(i) for i in range(07, 33)] # mXs r2
#           + ['RE{:02}'.format(i) for i in range(05, 32)] # mXs r3
#           + ['MA{:02}'.format(i) for i in range(37, 62)] # sXm r1
#           + ['MB{:02}'.format(i) for i in range(03, 30)] # sXm r2
#           + ['MC{:02}'.format(i) for i in range(07, 34)] # mxm
#           + ['SB{:02}'.format(i) for i in range(07, 34)] # sXs
#          )
module = '''module () {
        eval `$LMOD_CMD bash "$@"`
        }'''

rule rnaseq_map:
        input:
        output: "tmp"
        shell: "mkdir {output}"


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
