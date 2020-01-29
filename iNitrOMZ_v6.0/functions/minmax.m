function [norm_model norm_data] = minmax(constraints_model,constraints_data)

        % Finds min and max of either data or model prediction, as vertical profiles
	mmax = repmat(nanmax(nanmax(constraints_model,[],2),nanmax(constraints_data.val,[],2)),1,size(constraints_model,2));
	mmin = repmat(nanmin(nanmin(constraints_model,[],2),nanmin(constraints_data.val,[],2)),1,size(constraints_model,2));
        % Note: the max range (data or model) is given by: (mmax-mmin)

        % Removes minimum and normalizes by range
	norm_model    = (constraints_model - mmin)./(mmax - mmin);
	norm_data.val = (constraints_data.val - mmin)./(mmax - mmin);
