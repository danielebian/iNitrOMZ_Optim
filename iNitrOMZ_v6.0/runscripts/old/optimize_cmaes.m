% Here, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
 bgc1d_root='/data/project2/dbianchi/NitrOMZ/iNitrOMZ_v6.0/';
 addpath(genpath(bgc1d_root)); % adds root to MATLAB's search paths
% Adds CMAES subroutine:
 addpath('/data/project2/dbianchi/NitrOMZ/optimization/CMA_ES');

% Parameters to tune:
 remin = 0.08/86400;

% Matrix of parameters for optimization
% Format: 	name 			min_value		max_value
 AllParam = {
%		'wup_param',		0.5e-7,			15.0e-7; 	...
%		'poc_flux_top',		-30/86400,		-3/86400; 	...
		'KO2Rem',		0.01,			3.0;		...
		'KNO2No',		0.01,			0.5;		...
		'KO2Den1',		0.01,			8.0;		...
		'KO2Den2',		0.01,			8.0;		...
		'KO2Den3',		0.01,			8;		...
		'KNO3Den1',		0.01,			8;		...
		'KDen1',		remin/10,		remin;		...
		'KDen2',		remin/10,		remin;		...
		'KDen3',		remin/10,		remin;		...
		'KNO2Den2',		0.001,			0.5;		...
		'KN2ODen3',		0.001,			0.5;		...
		'KAx',			remin/10,		remin*2;	...
		'KNH4Ax',		0.01,			0.5;		...
		'KNO2Ax',		0.01,			0.5;		...
		'KO2Ax',		0.01,			8;		...
		};

 ParNames = AllParam(:,1);
 ParMin = [AllParam{:,2}]';
 ParMax = [AllParam{:,3}]';
 
 nPar = size(ParNames,2); % number of parameters
 
% Initialize final output structure
 Optim.ParNames = ParNames;
 Optim.ParMin = ParMin;
 Optim.ParMax = ParMax;

% NOTES: 
% (1) Parameters are normalized by subtracting the min and dividing by
%     a normalization factor, typically the range (so they are b/w 0-1)
%     This is done to allow the CMAES algorithm to work in the space [0 1]
% (2) If needed, remember to add the constraint:
%     Constraints: KDen1 + KDen2 + KDen3 = remin
%     This should be done in the cost function (bgc1d setup step)
%     as an ad-hoc constraint (removes one degree of freedom)
%     Remember to remove the corresponding K from the parameter pool!
%     (suggestion: remove KDen1, since first step drives everuthing)
% Calculates useful quantities for normalization, optimization, etc.
 ParMean = (ParMin + ParMax)/2';
 ParRange = ParMax - ParMin;
 ParNorm = ParRange;
 ParStart = (ParMean - ParMin) ./ ParNorm;
 ParSigma = ParRange./ParNorm/sqrt(12);

% Options
 optn.EvalParallel = '1';
 optn.LBounds = (ParMin - ParMin) ./ ParNorm;
 optn.UBounds = (ParMax - ParMin) ./ ParNorm;

% Enables parallelization
 if strcmp(optn.EvalParallel,'1')
    FunName = 'bgc1d_fc2minimize_cmaes_parallel';
    npar = 9;
    ThisPool = parpool(npar);
 else
   %FunName = 'bgc1d_fc2minimize_cmaes';
    FunName = 'bgc1d_fc2minimize_cmaes_parallel';
 end

% Runs the optimization
 FunArg.ParNames = ParNames;
 FunArg.ParMin = ParMin;
 FunArg.ParNorm = ParNorm;
 tic;
 [pvarout, pmin, counteval, stopflag, out, bestever] = cmaes(FunName,ParStart,ParSigma,optn,FunArg);

% Stops parallel pool
 if strcmp(optn.EvalParallel,'1')
    delete(ThisPool);
 end

% Fills in some output in final structure
 % Renormalized parameters
 Optim.ParOpt = ParMin + ParNorm .* pvarout;
 Optim.cmaes.options = optn; 
 Optim.cmaes.pvarout = pvarout; 
 Optim.cmaes.pmin = pmin; 
 Optim.cmaes.counteval = counteval; 
 Optim.cmaes.stopflag = stopflag; 
 Optim.cmaes.out = out; 
 Optim.cmaes.bestever = bestever; 
 Optim.RunTime = toc;

% Save ga output using today's date
 DateNow = bgc1d_getDate();
 save([bgc1d_root,'saveOut/ga_out_',DateNow,'.mat']);

