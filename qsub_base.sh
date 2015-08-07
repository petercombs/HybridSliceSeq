#!/bin/bash
#PBS -S /bin/bash
#PBS -N {job_name}
#PBS -l nodes=1:ppn=1 -l pmem=700mb -l mem=700mb
cd $PBS_O_WORKDIR

lastz Reference/dmel_masked/{id1}.fa Reference/dsim_masked/{id2}.fa O=400 E=30 K=2200 L=4000 H=2000 Y=3400 Q=prereqs/blastz.q > Reference/lav/{id1}-{id2}.lav
