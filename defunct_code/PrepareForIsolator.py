import pandas as pd
from GetASEStats import slices_per_embryo
from Utils import startswith, fbgns, pd_kwargs


if __name__ == "__main__":
    expr = pd.read_table('analysis_godot/summary.tsv', **pd_kwargs)
    n_slices = slices_per_embryo(expr)
    n_slices.pop('---', None)
    n_chunks = 10

    dirs_by_chunk = {emb.split('_rep')[0]: [[] for i in range(n_chunks)]
                     for emb in n_slices}

    out = open('tmp/isolator.yml', 'w')
    for line in open('analysis_godot/mst.log'):
        if '=' not in line: continue
        line = line.split()
        dirname = line[0]
        slice_id = line[-1]
        emb, sl_num = slice_id.split('_sl')
        genotype = emb.split('_rep')[0]

        chunk_num = (int(sl_num)-1) * n_chunks // n_slices[emb]
        dirs_by_chunk[genotype][chunk_num].append((slice_id, dirname))

    for genotype, chunks in dirs_by_chunk.items():
        for i, dirs in enumerate(chunks):
            out.write('\n{}-chunk{:02}:\n'.format(genotype.replace('_', '-'), i+1))
            for id, dir in dirs:
                out.write('    {}: {}.bam\n'.format(id.replace('_', '-'), dir))

    out.close()











