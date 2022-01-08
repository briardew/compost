%OBSPACK_WRITE_CSV  Write ObsPack model comparisons into CSV file for OCO-2
%   model intercomparison project

% Author: Brad Weir
%
% Changelog:
% 2021/02/18	Minor tweaks
%
% TODO:
%==============================================================================%

%fout = 'bweir_g5apr__co2_gvp_v21.csv';
%stag = 'obspack_co2_1_GLOBALVIEWplus_v2.1_2016-09-02';

fout = 'bweir_g5apr__co2_gvp_v31.csv';
stag = 'obspack_co2_1_GLOBALVIEWplus_v3.1_2017-10-18';

%fout = 'bweir_g5apr__co2_nrt_v32.csv';
%stag = 'obspack_co2_1_NRT_v3.2_2017-01-13';

%fout = 'bweir_g5apr__co2_nrt_v40.csv';
%stag = 'obspack_co2_1_NRT_v4.0_2017-09-08';

%fout = 'bweir_g5apr__co2_orcas.csv';
%stag = 'obspack_co2_1_ORCAS_v2.0_2017-04-05';

fid  = fopen(fout, 'w');
sout = [stag, '~%s~%d, %f\n'];

for ic = 1:NSITES
  fobs  = cell_fobs{ic};
  icut  = strfind(fobs, '.');

  dnobs = cell_dnobs{ic};
  nn0   = find(datenum(2014,01,01) <= dnobs, 1, 'first');

  gasmod = cell_gasmod{ic};

  for nn = nn0:numel(dnobs)
    fprintf(fid, sout, fobs(1:icut-1), cell_opnums{ic}(nn), gasmod(nn));
  end
end

fclose(fid);
