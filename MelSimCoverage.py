from progressbar import ProgressBar
try:
    import matplotlib.pyplot as mpl
    mpl.figure()
    mpl.plot([1,2])
    mpl.close()
    has_mpl = True
except:
    has_mpl = False

from pysam import Samfile

if __name__ == "__main__":
    window_size = int(10e3)
    mel = Samfile('analysis/on_mel/mel_gdna_bowtie2_dedup.bam')
    sim = Samfile('analysis/on_mel/sim_gdna_bowtie2_dedup.bam')

    mel_covs = []
    sim_covs = []
    prog = 0
    pbar = ProgressBar(max_value=sum(mel.lengths) + 1)
    for r, l in zip(mel.references, mel.lengths):
        for start in range(0, l, window_size):
            mel_covs.append(mel.count(r, start, start+window_size))
            sim_covs.append(sim.count(r, start, start+window_size))
        prog += l
        pbar.update(prog)
    pbar.finish()




