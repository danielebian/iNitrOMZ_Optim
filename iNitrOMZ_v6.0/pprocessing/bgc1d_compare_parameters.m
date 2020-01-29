function bgc1d_compare_parameters(bgc1,bgc2,imode)

 if nargin<3
    imode = 1;
 end

 switch imode
 case 1
    parnames = {'wup_param', ...
                'Kv_param', ...
                'b', ...
                'poc_flux_top', ...
                'Krem', ...
                'KAo', ...
                'KNo', ...
                'KDen1', ...
                'KDen2', ...
                'KDen3', ...
                'KAx', ...
                'KO2Rem', ...
                'KO2Ao', ...
                'KO2No', ...
                'KO2Den1', ...
                'KO2Den2', ...
                'KO2Den3', ...
                'KO2Ax', ...
                'KNH4Ao', ...
                'KNO2No', ...
                'KNO3Den1', ...
                'KNO2Den2', ...
                'KN2ODen3', ...
                'KNH4Ax', ...
                'KNO2Ax', ...
                'Ji_a', ...
                'Ji_b' ...
                };
 otherwise
    error(['Mode not supported']);
 end

 npar = length(parnames);
 allval = nan(npar,2);

 fprintf('%s %s %s',char(10),'----------------------------------------------------------',char(10));
 fprintf('%12s \t %8s \t %8s \t %8s','Variable','Run # 1 ','Run # 2 ','Ratio 2:1');
 fprintf('%s %s %s',char(10),'----------------------------------------------------------',char(10));
 for indp=1:npar
   thisp = parnames{indp};
   tratio = bgc2.(thisp)/bgc1.(thisp);
   sdisp = sprintf('%12s \t %0.3e \t %0.3e \t %0.3f', thisp, bgc1.(thisp), bgc2.(thisp),tratio); 
   disp(sdisp);
 end
 fprintf('%s %s','----------------------------------------------------------',char(10));
 


