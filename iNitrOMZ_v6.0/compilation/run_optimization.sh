#!/bin/csh -f
#$ -cwd
#  input           = /dev/null
#  output          = /u/scratch/d/danieleb
#$ -o /u/scratch/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/compilation/optimization_run01.joblog
#$ -j y
#$ -l h_data=4G,h_rt=30:00:00,highp
#$ -pe shared 12
#$ -M dbianchi@atmos.ucla.edu
#$ -m bea 
source /u/local/Modules/default/init/modules.csh
module load matlab
/u/scratch/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/compilation/FAOxAnox4_int_cf3nr_25k/optimize_cmaes/optimize_cmaes
