  function step_function = fstep(x,Xlim)

  % step function
  % this function is not differentiable at Xlim
  % and causes problems of convergence in the solver - should ot be used!

   Px = max(0,x);
   if Px>=Xlim
      step_function = 1; 
   else
      step_function = 0;
   end

  end
