  function oxygen_function = f2o2(o2,ko2,expn);

   o2pos = max(0,o2);
   oxygen_function = ko2.^expn./(o2pos.^expn+ko2.^expn);

  end
