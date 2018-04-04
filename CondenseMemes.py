from collections import defaultdict

def parse_into_pwms(fname):
    in_motif = False
    in_file = open(fname)
    while not in_motif:
        line = next(in_file)
        in_motif = line.startswith('MOTIF')

    all_tfs = defaultdict(list)
    nsites = 0
    motif = []
    # It has been pointed out that the FlyFactorSurvey version of hkb does have
    # predicted binding in the interesting region of hunchback. So I'm very
    # hackily including that as well.
    is_special_hkb = 'hkb_NAR' in line
    tf = line.strip().split()[-1].split('_')[0].lower()
    for line in in_file:
        if line.startswith('MOTIF'):
            if is_special_hkb:
                all_tfs['hkb_FFS'].append((nsites, motif))
            if nsites:
                all_tfs[tf].append((nsites, motif))
                nsites = 0
            motif = []
            tf = line.strip().split()[-1].split('_')[0].lower()
            is_special_hkb = 'hkb_NAR' in line
        elif line.startswith('MEME'):
            if nsites:
                all_tfs[tf].append((nsites, motif))
            if is_special_hkb:
                all_tfs['hkb_FFS'].append((nsites, motif))
            motif = []
            nsites = 0
        elif 'nsites' in line:
            line_data = line.split()
            nsites = int(line_data[line_data.index('nsites=') + 1])
        motif.append(line)

    return all_tfs




if __name__ == "__main__":
    in_file = 'prereqs/all_meme.meme'
    all_tfs = parse_into_pwms(in_file)
    out_file = open('Reference/all_meme_filtered.meme', 'w')

    for line in open(in_file):
        if line.startswith('MOTIF'):
            break
        out_file.write(line)

    for tf in all_tfs:
        best_motif = max(all_tfs[tf])[1]
        out_file.writelines(best_motif)

    out_file.close()

