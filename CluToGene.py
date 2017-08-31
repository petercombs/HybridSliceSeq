def get_exon_pairs(gtf='Reference/mel_r5_good.gtf', key='gene_name'):
    exon_pairs = {}
    last_type = ''
    last_stop = -1
    for line in open(gtf):
        data = line.split('\t')
        chrom = data[0]
        line_type = data[2]
        start = int(data[3])
        stop = int(data[4])
        if last_type == 'exon' and line_type=='exon':
            annot = {
                d.strip().split()[0]: d.strip().split(maxsplit=1)[1].strip('";')
                for d in data[-1].split(';')
                if d.strip()
            }
            exon_pairs[(chrom, last_stop, start)] = annot.get(key)
        last_stop = stop
        last_type = line_type
    return exon_pairs


