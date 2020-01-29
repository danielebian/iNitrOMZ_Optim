 function mcost = cost_function_problem(pvar,Nd);
 
 % Defines the entire problem to be optimized, finding the set of pvar parameters
 % that minimize a cost function
 Nx = 50; 
%Nd = 4;
 tmp = load(['dout_data' num2str(Nd) 'd']);
 dout = tmp.dout;
 
 % Solve model for current parameter vector
 out = nd_model(Nd,Nx,pvar);

 if (0)
    figure(1)
    clf;
    subplot(2,1,1)
    pcolor(out);shading flat
    subplot(2,1,2)
    pcolor(dout);shading flat
    drawnow
 end
 
 % Calculates cost function
 mcost = cost_function(out,dout);

 
