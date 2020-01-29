 function  mcost = cost_function(out,dout);

 % Interpolates out on dout grid, and finds the RMSE
 Nd = ndims(out); 
 Nx = size(out,1); 
 Ny = size(dout,1); 

 xi = linspace(-1,1,Nx);
 yi = linspace(-1,1,Ny);
 
 switch Nd
 case 2
   [x1 x2] = ndgrid(xi,xi);
   [y1 y2] = ndgrid(yi,yi);
   iout = interpn(x1,x2,out,y1,y2);
   diff2 = (dout-iout).^2;
   mcost = mean(diff2(:)); 
 case 3
   [x1 x2 x3] = ndgrid(xi,xi,xi);
   [y1 y2 y3] = ndgrid(yi,yi,yi);
   iout = interpn(x1,x2,x3,out,y1,y2,y3);
   diff2 = (dout-iout).^2;
   mcost = mean(diff2(:)); 
 case 4
   [x1 x2 x3 x4] = ndgrid(xi,xi,xi,xi);
   [y1 y2 y3 y4] = ndgrid(yi,yi,yi,yi);
   iout = interpn(x1,x2,x3,x4,out,y1,y2,y3,y4);
   diff2 = (dout-iout).^2;
   mcost = mean(diff2(:)); 
 otherwise
   error('Number of dimensions not supported yet');
 end

