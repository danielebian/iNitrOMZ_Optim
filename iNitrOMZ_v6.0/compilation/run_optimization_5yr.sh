#!/bin/csh -f
#$ -cwd
#  input           = /dev/null
#  output          = /u/scratch/d/danieleb
#$ -o /u/home/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/compilation/optimization_run_5nr_10.joblog
#$ -j y
#$ -l h_data=4G,h_rt=30:00:00,highp
#$ -pe shared 12
#$ -M dbianchi@atmos.ucla.edu
#$ -m bea 
source /u/local/Modules/default/init/modules.csh
module load matlab
/u/home/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/compilation/FAOxAnox5_ni_cf3_yr_25k/optimize_cmaes/optimize_cmaes
