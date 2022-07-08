function pickmap(NTYPE)

if (NTYPE == 2)
  worldmap('north pole');

  setm(gca, 'maplatlimit', [0 90]);

  setm(gca, 'mlabelparallel', 10);
  setm(gca, 'mlabellocation', [-135,-90,-45,0,45,90,135,180]);
  setm(gca, 'plabelmeridian', 0);
  setm(gca, 'plabellocation', [30,60,90]);

  setm(gca, 'plinelocation',   15);

  tightmap;

elseif (NTYPE == 3)
  worldmap('south pole');

  setm(gca, 'maplatlimit', [-90 0]);

  setm(gca, 'mlabelparallel', -10);
  setm(gca, 'mlabellocation', [-135,-90,-45,0,45,90,135,180]);
  setm(gca, 'plabelmeridian', 180);
  setm(gca, 'plabellocation', [-90,-60,-30]);

  setm(gca, 'plinelocation',   15);

  tightmap;

elseif (NTYPE == 4)
  worldmap('north pole');

  setm(gca, 'maplatlimit', [50 90]);

% setm(gca, 'mlabelparallel', 10);
% setm(gca, 'mlabellocation', []);
  setm(gca, 'plabelmeridian', 45);
  setm(gca, 'plabellocation', [60,75,90]);

  setm(gca, 'plinelocation',   15);
  setm(gca, 'mlinelocation',   45);

  mlabel;
  setm(gca, 'origin', [90 -135 0]);
  tightmap;

elseif (NTYPE == 5)
  worldmap('north america');

  setm(gca, 'origin', [0 -150 0]);
  setm(gca, 'maplatlimit', [50 80]);
  setm(gca, 'maplonlimit', [-190 -110]);

  setm(gca, 'plabelmeridian', -120);
  setm(gca, 'plabellocation', [70,60]);
  setm(gca, 'plinelocation',  10);

  tightmap;

elseif (NTYPE == 6)
  usamap('conus');

elseif (NTYPE == 7)
  worldmap('north pole');

  setm(gca, 'maplatlimit', [30 90]);

  setm(gca, 'mlabelparallel', 10);
  setm(gca, 'mlabellocation', [-135,-90,-45,0,45,90,135,180]);
  setm(gca, 'plabelmeridian', 0);
  setm(gca, 'plabellocation', [45,60,75,90]);

  setm(gca, 'plinelocation',   15);

  tightmap;

else
  worldmap('world');

  if (NTYPE == 97)
    setm(gca, 'mapprojection', 'winkel');
    tightmap;
  elseif (NTYPE == 98)
    setm(gca, 'mapprojection', 'mollweid');
    tightmap;
  elseif (NTYPE == 99)
    setm(gca, 'mapprojection', 'mercator');
    tightmap;
  end
end
