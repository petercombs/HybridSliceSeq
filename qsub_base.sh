#!/bin/bash
#PBS -S /bin/bash
#PBS -N {job_name}
#PBS -o logs/{id1}-{id2}.o.log
#PBS -e logs/{id1}-{id2}.e.log
#PBS -l nodes=1:ppn=1 -l pmem=700mb -l mem=700mb
cd $PBS_O_WORKDIR

date
/home/pcombs/bin/lastz Reference/d{species1}_masked/{id1}.fa Reference/d{species2}_masked/{id2}.fa O=400 E=30 K=2200 L=4000 H=2000 Y=3400 Q=prereqs/blastz.q > Reference/lav/{id1}-{id2}.lav
date
