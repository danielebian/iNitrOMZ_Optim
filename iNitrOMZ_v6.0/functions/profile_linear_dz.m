 function y = profile_linear_dz(x,x0,x1,y0,y1)
 % Produce a linear interpolation between two points P1 and P2
 % Arguments:
 % P1:(x0,y0)
 % P2:(x1,y1)
 % P :(x,y) 
 % Usage: y = profile_linear(x,x0,x1,y0,y1)

 y =  ( (y1-y0)/(x1 - x0) )   * ones(size(x)) ;


