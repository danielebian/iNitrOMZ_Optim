 function mcost = cost_function_problem3(pvar,Nd);
 
 % Defines the entire problem to be optimized, finding the set of pvar parameters
 % that minimize a cost function
 Nx = 50; 
%Nd = 4;
 tmp = load(['dout_data' num2str(Nd) 'd']);
 dout = tmp.dout;
 
 % Allows for parallelization of cost function, accepting pvar input of size NpxNm
 % Np: # parameters of the problem
 % Nm: # simultaneous calls

 [Np Nm] = size(pvar);

 disp(['Internal iteration #: ' num2str(Nm)]);

 mcost = nan(1,Nm);
 
 parfor indm=1:Nm
    tpvar = pvar(:,indm);
    % Solve model for current parameter vector
    out = nd_model(Nd,Nx,tpvar);
    % Calculates cost function
    tcost = cost_function(out,dout);
    mcost(1,indm) = tcost;
 end

 
