
 hoffman2

 qrsh -l highp,h_rt=48:00:00,h_data=4G -pe shared 12

 module load matlab
 cd $SCRATCH/NitrOMZ/iNitrOMZ_v6.0/
 matlab -nodesktop -nosplash

 addpath(genpath('/u/scratch/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/'))
 optimize_cmaes


