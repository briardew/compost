%OBSPACK_WRITE_CSV  Write ObsPack model comparisons into CSV file for OCO-2
%   model intercomparison project

% Author(s):	Brad Weir <brad.weir@nasa.gov>
%
% Changelog:
% 2021-02-18	Minor tweaks
%
% TODO:
%===============================================================================

%fin  = 'm2cc_ana.obspack_co2_nrt_v6.mat';
%fout = 'bweir_m2cc_co2_nrt_v611.csv';
%stag = 'obspack_co2_1_NRT_v6.1.1_2021-05-17';

fin  = 'm2cc_ana.obspack_co2_mip_v3.mat';
fout = 'bweir_m2cc_co2_oco2mip_v32.csv';
stag = 'obspack_co2_1_GLOBALVIEWplus_v6.1_2021-03-01';

load(fin);
fid  = fopen(fout, 'w');

for ic = 1:NSITES
  fobs  = cell_fobs{ic};
  icut  = strfind(fobs, '.');

  sout = stag;
  if strcmp(fobs, 'co2_aircorenoaa_aircore_1_allvalid.nc')
    sout = 'obspack_co2_1_AirCore_v4.0_2020-12-28';
  end
  if strcmp(fobs, 'co2_con_aircraft-flask_42_allvalid.nc') | ...
     strcmp(fobs, 'co2_con_aircraft-insitu_42_allvalid.nc')
    sout = 'obspack_co2_1_CONTRAIL_v1.0_2021-09-13';
  end
  if strcmp(fobs, 'co2_alf_aircraft-pfp_433_representative.nc')   | ...
     strcmp(fobs, 'co2_pan_aircraft-pfp_433_representative.nc')   | ...
     strcmp(fobs, 'co2_rba-b_aircraft-pfp_433_representative.nc') | ...
     strcmp(fobs, 'co2_san_aircraft-pfp_433_representative.nc')   | ...
     strcmp(fobs, 'co2_tef_aircraft-pfp_433_representative.nc')
    sout = 'obspack_co2_1_INPE_RESTRICTED_v2.0_2018-11-13';
  end
  if strcmp(fobs, 'co2_ah2_shipboard-insitu_20_allvalid.nc')  | ...
     strcmp(fobs, 'co2_ftw_shipboard-insitu_20_allvalid.nc')  | ...
     strcmp(fobs, 'co2_ftws_shipboard-insitu_20_allvalid.nc') | ...
     strcmp(fobs, 'co2_gw_shipboard-insitu_20_allvalid.nc')   | ...
     strcmp(fobs, 'co2_nc2_shipboard-insitu_20_allvalid.nc')  | ...
     strcmp(fobs, 'co2_px_shipboard-insitu_20_allvalid.nc')   | ...
     strcmp(fobs, 'co2_sk_shipboard-insitu_20_allvalid.nc')   | ...
     strcmp(fobs, 'co2_tf1_shipboard-insitu_20_allvalid.nc')  | ...
     strcmp(fobs, 'co2_tf5_shipboard-insitu_20_allvalid.nc')
    sout = 'obspack_co2_1_NIES_Shipboard_v3.0_2020-11-10';
  end
  if strcmp(fobs, 'co2_man_aircraft-insitu_1_allvalid.nc')
    sout = 'obspack_multi-species_1_manaus_profiles_v1.0_2021-05-20';
  end

  dnobs = cell_dnobs{ic};
  nn0   = find(datenum(2014,01,01) <= dnobs, 1, 'first');

  gasmod = cell_gasmod{ic};

  for nn = nn0:numel(dnobs)
    fprintf(fid, '%s~%s~%d, %f\n', sout, fobs(1:icut-1), ...
            cell_opnums{ic}(nn), gasmod(nn));
  end
end

fclose(fid);
