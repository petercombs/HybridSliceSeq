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


tfs = "OTF0070.1 hb  kni gt D FBgn0001325_4 tll hkb  FBgn0000251_3 FBgn0003448_3 twi OTF0532.1 OTF0478.1 OTF0181.1"
tf_names = "bcd  hb  kni gt D Kr            tll hkb cad            sna           twi zen       Doc2      rho/pnt   "


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
            -v -X 100 --show-alignment \
            -x .75 --y-sep 60 --y-scale 3.5 \
            --comp1 {wildcards.A} --comp2 {wildcards.B} \
            --match-dim 0.7 \
            --meme-suite \
            --tf {tfs} --tf-names {tf_names} \
            --sequence --fasta {input.fa} \
            --needleall {input.aln} \
            --coordinates-bed {input.bed} \
            analysis/targets/{wildcards.gene}/
            """

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
