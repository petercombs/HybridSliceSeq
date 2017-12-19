from Bio import SeqIO, Seq
from Bio.Seq import reverse_complement as rc
from progressbar import ProgressBar, UnknownLength
import pandas as pd


# IUPAC ambiguous base pairs
ambigs = {
    ('A', 'T'): 'W',
    ('C', 'G'): 'S',
    ('A', 'C'): 'M',
    ('G', 'T'): 'K',
    ('A', 'G'): 'R',
    ('C', 'T'): 'Y',
}


# Eventually maybe I want to not hard code hunchback and give the opportunity to
# work on generic genes, but I'm not there yet...
gene_coords = [('dmel_3R', 4516702, 4519894),
             ('dmel_3R', 4520178, 4520322)]
gene_strand = '-'

if __name__ == "__main__":
    chroms = {
        rec.id: rec.seq.tomutable()
        for rec in SeqIO.parse('Reference/dmelr5_prepend.fasta',
                               'fasta', Seq.Alphabet.generic_dna)
    }

    mel_copy = {
        rec: Seq.MutableSeq(str(chroms[rec]), Seq.Alphabet.generic_dna)
        for rec in chroms

    }

    sim_copy = {
        rec: Seq.MutableSeq(str(chroms[rec]), Seq.Alphabet.generic_dna)
        for rec in chroms

    }

    # Replace ambiguous bases
    variants = 'analysis_godot/on_melr5/melsim_variant_hb.bed'
    variants_in_gene = {}
    pb = ProgressBar(maxval=UnknownLength)
    for line in pb(open(variants)):
        chrom, pos, _, variant = line.split()
        pos = int(pos)
        for gene_chrom, low, high in gene_coords:
            if chrom == gene_chrom and low <= pos <= high:
                variant = variant.split('|')
                mel, sim = variant
                #assert chroms[chrom][pos] == 'N'
                chroms[chrom][pos] = ambigs[tuple(sorted(variant))]
                mel_copy[chrom][pos] = mel
                sim_copy[chrom][pos] = sim
                assert mel_copy[chrom][pos] != sim_copy[chrom][pos]
                if gene_strand == '-':
                    variant = (rc(variant[0]),
                                rc(variant[1]))
                variants_in_gene[pos] = variant
                break


    # Build the hunchback transcript
    gene_str = []
    mel_str = []
    sim_str = []
    gene_bps = []
    for chrom, start, stop in gene_coords:
        gene_str.append(str(chroms[chrom][start:stop]))
        mel_str.append(str(mel_copy[chrom][start:stop]))
        sim_str.append(str(sim_copy[chrom][start:stop]))
        gene_bps.extend(range(start, stop))
    gene_str = ''.join(gene_str)
    mel_str = ''.join(mel_str)
    sim_str = ''.join(sim_str)
    if gene_strand == '-':
        gene_str = rc(gene_str)
        mel_str = rc(mel_str)
        sim_str = rc(sim_str)
        gene_bps = list(reversed(gene_bps))
    print(''.join(gene_str))

    primers = pd.read_table('analysis/targets/hb/ASprimers.tsv')

    new_primers = pd.DataFrame(columns=['seq', 'orientation', '5p_end',
                                        'length', 'species', 'Tm', 'n_snps'],
                              index=primers.index,
                              )
    for i, primer in primers.iterrows():
        seq = primer['Primer Seq']
        length = primer.Length
        new_primer_seq = []
        n_snps = 1
        if primer.Orientation == 'FORWARD':
            pos = primer.Start - 1 # Fix 0/1 based indexing
            primer_bps = gene_bps[pos:pos+length]
            refalt = variants_in_gene[primer_bps[-1]].index(seq[-1])
            assert len(primer_bps) == len(seq)
            for bp, base in zip(primer_bps, seq):
                if base in 'ACGT':
                    new_primer_seq.append(base)
                else:
                    new_primer_seq.append(variants_in_gene[bp][refalt])
                    n_snps += 1
        elif primer.Orientation == 'REVERSE':
            pos = primer.Pos - 1
            primer_bps = list(reversed(gene_bps[pos: pos + length]))
            refalt = variants_in_gene[gene_bps[pos]].index(rc(seq[-1]))
            assert len(primer_bps) == len(seq)
            for bp, base in zip(primer_bps, seq):
                if base in 'ACGT':
                    new_primer_seq.append(base)
                else:
                    new_primer_seq.append(rc(
                        variants_in_gene[bp][refalt]
                    ))
                    n_snps += 1
            pos = pos + length


        else:
            raise ValueError("Primer direction: {}".format(primer.Orientation))
        new_primers.ix[i] = {
            'seq':''.join(new_primer_seq),
            'orientation':primer.Orientation,
            '5p_end':pos,
            'length':length,
            'species': ['mel', 'sim'][refalt],
            'Tm':primer.Tm,
            'n_snps': n_snps,
        }

    new_primers.to_csv('analysis/targets/hb/phased_primers.tsv', sep='\t')




