%OBSPACK_WRITE_DEI  Write ObsPack model comparisons into file compatible
%   with NOAA DEI utilities

% Author: Brad Weir
%
% Changelog:
% 2018/02/07	Adding ability to handle generic trace gas
%
% TODO:
%==============================================================================%

for ic = 1:NSITES
  fobs  = cell_fobs{ic};

  if (isempty(strfind(fobs, 'flask_1'))), continue; end

  gasmod = cell_gasmod{ic};
  iok = find(~isnan(gasmod));

  if (numel(iok) == 0), continue; end

  dnobs = cell_dnobs{ic};
  dvecs = datevec(dnobs);
  years = dvecs(:,1);
  decyr = years + (dnobs - datenum(years,01,01)) ./ ...
                  (datenum(years+1,01,01) - datenum(years,01,01));

  if (isempty(strfind(fobs, 'poc_shipboard')))
    fout = ['dei/', fobs(5:7), '_01D0_dat.co2'];
    fid  = fopen(fout, 'w');

    fprintf(fid, ' %7.4f   %7.4f\n', [decyr(iok), gasmod(iok)]');

    fclose(fid);
  else
    for latb = -45:5:45
      ibin = find(-2.5 <= cell_lat{ic} - latb & cell_lat{ic} - latb < 2.5);
      iuse = intersect(iok, ibin);

      if (numel(iuse) == 0), continue; end

      if (latb <  0), fhem = 's'; end
      if (latb == 0), fhem = '0'; end
      if (0 <  latb), fhem = 'n'; end

      fbin = [fhem, num2str(abs(latb),'%02u')];
      fout = ['dei/', fobs(5:7), fbin, '_01D0_dat.co2'];
      fid  = fopen(fout, 'w');

      fprintf(fid, ' %7.4f   %7.4f\n', [decyr(iuse), gasmod(iuse)]');

      fclose(fid);
    end
  end
end
