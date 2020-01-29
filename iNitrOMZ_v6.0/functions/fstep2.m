  function step_function = fstep2(x,Xlim,Xwidth)

  % Approximate a step function with a hyperbolic tangent.
  % A true step function is not derivable and could cause some problems 
  % depending the equatuions
  
  if Xwidth~=0
     Px = max(0,x);
     step_function = (tanh((x-Xlim)/Xwidth)+1)/2; 
  else
     % True step fuunction
     step_function = stepfun(x,Xlim);
  end
