  function oxygen_function = fo2(o2o,ko2,expn);

   Po2o = max(0,o2o);
   oxygen_function = Po2o.^expn./(Po2o.^expn+ko2.^expn);

  end
