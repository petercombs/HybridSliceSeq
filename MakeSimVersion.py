import pandas as pd
from Bio import SeqIO, Seq
from os import path

def replace_variant(variant, gt, sequence):
    var = variant[gt].split('/')[0]
    pos = variant.POS - 1
    if var == '*':
        sequence.pop(pos)
        return -len(variant.REF)
    elif var == '.' or (sequence[pos:pos + len(variant.REF)]
                        == var):
        return 0
    elif len(var) == len(variant.REF):
        sequence[pos:pos+len(var)] = var
        return 0
    else:
        print(sequence[pos:pos+len(variant.REF)])
        for i, c in enumerate(variant.REF):
            sequence.pop(pos)
        for c in reversed(var):
            sequence.insert(pos, c)
        return len(variant.REF) - len(var)







if __name__ == "__main__":
    mel_genome = {
        rec.id: rec.seq
        for rec in SeqIO.parse('prereqs/dmel-all-chromosome-r6.06.fasta', 'fasta')
    }


    regions = pd.read_table('prereqs/dnase_acc_spatial.tsv', header=None)
    vars = pd.read_table('analysis/on_mel/melsim_variants.tsv')

    for region_id in regions.index:
        region = regions.ix[region_id, 1]
        name = regions.ix[region_id, 2]
        if name == pd.np.nan:
            name = region.replace(':', '_')
        print(name)
        chrom = region.split(':')[0]
        mel_chr = Seq.MutableSeq(str(mel_genome[chrom]))
        sim_chr = Seq.MutableSeq(str(mel_genome[chrom]))
        coords = region.split(':')[1].split('..')
        start = int(coords[0])
        stop = int(coords[1])
        mel_diff = 0
        sim_diff = 0
        region_vars = vars.query('CHROM == "dmel_{}" and {} < POS and POS <= {}'
                                 .format(chrom, start, stop))
        lowest = 1e10
        for var_ix in reversed(region_vars.index):
            var = region_vars.ix[var_ix]
            assert lowest > var.POS
            lowest = var.POS
            mel_diff += replace_variant(var, 'mel_gdna.GT', mel_chr)
            sim_diff += replace_variant(var, 'sim_gdna.GT', sim_chr)

        SeqIO.write([SeqIO.SeqRecord(mel_chr[start:stop+mel_diff], 'mel_'+region),
                     SeqIO.SeqRecord(sim_chr[start:stop+sim_diff], 'sim_'+region)],
                    path.join('analysis', 'results', name+'.fasta'),
                    'fasta')




