 function Cost = Cost_quad_std_weighted(constraints_model,constraints_data,depth_weights,iplot)
% ========================================================================
%
% File     : Cost.m
% Date     : September 2017
% Author   : Simon Yang
% Function : compute individual cost function 
% ========================================================================
% Quadratic cost, normalized at each depth by the standard deviation or range of the data+model,
% overall cost normalized by number of points
% sum((z-z_data)^2/std^2)/n_data

 if nargin<4
    iplot = 0;
 end

% If no standard deviation, use 20% of value at that point

% Removes minimum and normalizes by range at each depth
% here uses range from model+data
%[norm_model, norm_data] = minmax(constraints_model, constraints_data);
% Removes minimum and normalizes by range at each depth
% here uses range from data
 [norm_model, norm_data] = minmax_data(constraints_model, constraints_data);
    
% Get a RMSE of normalized variables
 err_sq = (norm_model - norm_data.val).^2;
 data_mask = double(~isnan(norm_data.val));

% Calculates Cost, weighting square errors by depth-dependent weights
%Cost = nansum(err_sq.*depth_weights,2)./nansum(data_mask.*depth_weights,2);
 w_err_sq = bsxfun(@times,err_sq,depth_weights);
 w_data_mask = bsxfun(@times,data_mask,depth_weights);
 Cost = nansum(w_err_sq,2)./nansum(w_data_mask,2);

% Weights NaN values in model by an arbitrary large amounts (here 1000)
 inan = any(isnan(norm_model),2);
 Cost(inan) = 1000;

 % Some diagnostics
 if (iplot>=1)
    % 'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2' 'nstar'
    vplot = [1 2 4 5 7 9];
    zlev = [1:size(constraints_model,2)];
    if (iplot==2)
       ff = figure('visible','off');
    else
       ff = figure;
    end
    for indp=1:length(vplot)
       %--------------------
       subplot(2,3,indp)
       plot(norm_model(vplot(indp),:),-zlev,'.b-','linewidth',3,'markersize',3)
       hold on
       plot(norm_data.val(vplot(indp),:),-zlev,'.r-','linewidth',3,'markersize',3) 
       title([constraints_data.name{vplot(indp)} ' : ' num2str(Cost(vplot(indp)))]);
    end
    if iplot==2
       cost_pre = Cost;
       cost_pre(isnan(cost_pre)) = 0;
       %tcost = nansum((cost_pre.* constraints_data.weights').^2)/(nansum(Data.weights)+nansum(Data.rates.weights))^2;
       tcost = nansum(cost_pre.* constraints_data.weights')/(nansum(constraints_data.weights));
       DateNow = bgc1d_getDate();
       tmpname = [DateNow '_c_' num2str(tcost,3) '.jpg']; 
       mprint_fig('name',tmpname,'for','jpeg','sty','nor1','silent',1);
       close(ff); 
    end
 end

