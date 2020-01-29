 function y = profile_sigmoid_dz(z,zf,Hf)
% Gives the z-derivative of the sigmoid profile:
% Gives a continuous derivable profile formed by
% a sigmoid  that grows from 0 to 1 
% with flex at zf and scale Hf

 y = 8/Hf*(exp(-(z-zf)/(Hf/8)))./(1+exp(-(z-zf)/(Hf/8))).^2;


end



