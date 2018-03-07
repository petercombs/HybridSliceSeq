#! /usr/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time=6-00:00:00
#SBATCH --mem=30G
#SBATCH -p owners
#SBATCH --workdir=/home/users/pcombs/HybridSliceSeq
#SBATCH -o /home/users/pcombs/HybridSliceSeq/logs/FDR_PSI.out
#SBATCH -e /home/users/pcombs/HybridSliceSeq/logs/FDR_PSI.err

echo "Starting Script"

module () {
                eval $($LMOD_CMD bash "$@")
                        [ $? = 0 ] && eval $(${LMOD_SETTARG_CMD:-:} -s sh)
                }


export MODULEPATH=/share/PI/hbfraser/modules/modules:${MODULEPATH}
module load fraserconda

echo "Modules loaded"

export PYTHONPATH=$HOME/fyrd:$PYTHONPATH:/home/users/pcombs/HybridSliceSeq

python FindAutocorrPSI.py
