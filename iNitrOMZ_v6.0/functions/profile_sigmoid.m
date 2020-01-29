
 function y = profile_sigmoid(z,zf,Hf)
% Gives a continuous derivable profile formed by
% a sigmoid  that grows from 0 to 1 
% with flex at zf and scale Hf

 y = 1./(1+exp(-(z-zf)/(Hf/8)));


end



