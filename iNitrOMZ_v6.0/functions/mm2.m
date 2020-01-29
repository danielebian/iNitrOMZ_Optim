  function fmen2 = mmen2(var,Kvar);

   Pvar = max(0,var);
   fmen2 = Kvar./(Pvar+Kvar);

  end
