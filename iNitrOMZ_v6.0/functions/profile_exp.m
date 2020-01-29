 function y = profile_exponential(z,z0,zefold)
 % Produce an exponential profile passing through (z0,1) 
 % with e-folding scale zfold
 % Usage: y =  profile_exponential(z,z0,zefold)

 y = exp( (z-z0)/zefold );

 end


