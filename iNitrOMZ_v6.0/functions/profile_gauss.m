 function c0 = set_initial_profile(z,max,loc,std)

 % usage
 % phyto = profile_gauss(z,100,0.06,0.04)

 c0 = max*exp( -(z-(z(1)+z(end)*loc)).^2/(((z(1)+z(end)*std)^2)));

 end


