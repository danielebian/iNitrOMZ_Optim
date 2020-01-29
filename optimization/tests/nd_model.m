 function out = nd_model(Nd,Nx,pvar);

 % Evaluates a ND model given by a quadratic form:
 % f(xi,pj) = p1*x1 + p2*x2 + p3*x3 + ...
 %            p4*x1*x1 + 54*x2*x2 + p6*x3*x3 + ...
 %            p7*x1*x2 + p8*x2*x3 + p9*x2*x3            
 % For any given dimensionality Nd, the number of parameters Np=(3*Nd+Nd^2)/2
 % Nd -> linear terms
 % Nd -> quadratic terms
 % Nd*(Nd-1)/2 -> interaction terms
 % The variables xi are defined in the range [-1,1], discretized by Nx steps

 % Checks the size of Np is right for the problem
 Np = length(pvar);
 Np_expected = (3*Nd + Nd*Nd)/2;
 if Np~=Np_expected;
   error('Number of parameters does not match dimensionality of the problem');
 end

 xi = linspace(-1,1,Nx); 
 % Could be generalized but it's a pain...
 switch Nd
 case 2
   [x1 x2] = ndgrid(xi,xi); 
   quadform = pvar(1)*x1 + pvar(2)*x2 + ...
              pvar(3)*x1.^2 + pvar(4)*x2.^2 + ...
              pvar(5)*x1.*x2;
   taper = ((cos(x1*pi)+1)/2).^2 .* ((cos(x2*pi)+1)/2).^2;
   out = quadform .* taper;
 case 3
   [x1 x2 x3] = ndgrid(xi,xi,xi); 
   quadform = pvar(1)*x1 + pvar(2)*x2 + pvar(3)*x3 + ...
              pvar(4)*x1.^2 + pvar(5)*x2.^2 + pvar(6)*x3.^2 +...
              pvar(7)*x1.*x2 + pvar(8)*x1.*x3 + pvar(9)*x2.*x3;
   taper = ((cos(x1*pi)+1)/2).^2 .* ((cos(x2*pi)+1)/2).^2 .* ((cos(x3*pi)+1)/2).^2;
   out = quadform .* taper;
 case 4
   [x1 x2 x3 x4] = ndgrid(xi,xi,xi,xi); 
   quadform = pvar(1)*x1 + pvar(2)*x2 + pvar(3)*x3 + pvar(4)*x4 + ...
              pvar(5)*x1.^2 + pvar(6)*x2.^2 + pvar(7)*x3.^2 + pvar(8)*x4.^2 + ...
              pvar(9)*x1.*x2 + pvar(10)*x1.*x3 + pvar(11)*x1.*x4 + pvar(12)*x2.*x3 + pvar(13)*x2.*x4 + pvar(14)*x3.*x4;
   taper = ((cos(x1*pi)+1)/2).^2 .* ((cos(x2*pi)+1)/2).^2 .* ((cos(x3*pi)+1)/2).^2 .* ((cos(x4*pi)+1)/2).^2;
   out = quadform .* taper;
 otherwise
   error('Number of dimensions not supported yet');
 end
 

              

