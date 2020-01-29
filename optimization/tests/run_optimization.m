 addpath('/data/project2/dbianchi/NitrOMZ/optimization/CMA_ES');
%opt_method = 'ga';
 opt_method = 'CMA';
%opt_method = 'CMA_parallel';

 Nd = 3;	% dimensionality
 Ny = 10;	% # points for data

 Np = (3*Nd + Nd*Nd)/2; 
 
 % Target solutions
 switch Nd
 case 2
    pvar = [0.5 -0.5 0.5 -0.5 0.5];
 case 3
    pvar = [0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5];
 case 4
    pvar = [0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5];
 otherwise
   error('Number of dimensions not supported yet');
 end

 dout = nd_model(Nd,Ny,pvar);
 save(['dout_data' num2str(Nd) 'd'],'dout');  

 tic;
 switch opt_method
 case 'ga'
    % Sets up the optimization : Genetic Algorithm
    % Sets up the cost function that is passed to optimization algorithm
    costfunc = @(x)cost_function_problem(x,Nd);
    passfunct = @(x)costfunc(x);
    options = optimoptions('ga','FunctionTolerance',1e-6,'MaxStallGenerations',50, 'MaxGenerations', 200, ...
                           'PlotFcn', @gaplotbestf,'UseParallel', false, 'UseVectorized', false,'Display','iter');
    pvarout = ga(passfunct,Np,[],[],[],[],[],[],[],options);
 case 'CMA'
    % Here uses the as cost function: cost_function_problem2, which has hard-wired the # dimensions
    x0 = zeros(1,Np);
    sigma = 0.7 + zeros(1,Np);
    optn.EvalParallel = '0';
    [pvarout, pmin, counteval, stopflag, out, bestever] = cmaes('cost_function_problem2', x0(:), sigma(:),optn,Nd)
 case 'CMA_parallel'
    % Here uses the as cost function: cost_function_problem3, which has hard-wired the # dimensions
    % and is designed to be parallelized
    ppool = parpool(11);
    x0 = zeros(Np,1);
    sigma = 0.7 + zeros(Np,1);
    optn.EvalParallel = '1';
    [pvarout, pmin, counteval, stopflag, out, bestever] = cmaes('cost_function_problem3', x0, sigma,optn,Nd)
 end

 runtime = toc;
 disp(pvarout);
 disp(['runtime : ' num2str(runtime)]);

