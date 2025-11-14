function [ZZ,NN,SS] = sat2grd(grdlat, grdlon, satlat, satlon, zz)

[LA,LO] = meshgrid(grdlat, grdlon);
ZZ = zeros(size(LA));
NN = zeros(size(LA));
SS = zeros(size(LA));

for ii = 1:numel(zz)
  dLA = abs(grdlat - satlat(ii));
  dLO = abs(grdlon - satlon(ii));

  iLA = find(dLA == min(dLA),1);
  iLO = find(dLO == min(dLO),1);

  NN(iLO,iLA) = NN(iLO,iLA) + 1;

  del  = zz(ii) - ZZ(iLO,iLA);
  ZZ(iLO,iLA) = ZZ(iLO,iLA) + del/NN(iLO,iLA);

  del2 = zz(ii) - ZZ(iLO,iLA);
  SS(iLO,iLA) = SS(iLO,iLA) + del*del2;
end
SS = SS./(NN - 1);

ZZ(NN == 0) = NaN;
SS(NN == 0) = NaN;
