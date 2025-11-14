%MAKE_BIAS  Compare and plot model bias at NOAA ObsPack locations
%
%   I don't know

% Author: Brad Weir
%
% Changelog:
% 2018/04/06	Site is based on string now
%
% TODO:
%==============================================================================%

% Anomaly detection (El Nino/growth rate)
% ---------------------------------------
%sname = 'co2_mlo_surface-insitu';			% Mauna Loa		(MLO)
%sname = 'co2_smo_surface-insitu';			% American Samoa	(SMO)

% Anomaly detection (Droughts)
% ----------------------------
% CASA-GFED already does a pretty good job picking this up; apparent from
% exceptionally high trough in summer (inhibited drawdown)
%sname = 'co2_lef_tower-insitu_1_allvalid-396magl';	% Park Falls		(LEF)
%sname = 'co2_sgp_surface-insitu_64_allvalid-60magl';	% Lamont		(SGP)
%sname = 'co2_wbi_tower-insitu_1_allvalid-379magl';	% West Branch		(WBI)

% Anomaly detection (Fires)
% -------------------------
% Have to fly through the plume, probably need CO
%sname = 'co2_arctas_aircraft-insitu_428_allvalid-dc8';	% ARCTAS  DC8 (2008)
%sname = 'co2_seac4rs_aircraft-insitu_428_allvalid-dc8';% SEAC4RS DC8
%sname = 'co2_seac4rs_aircraft-insitu_428_allvalid-ER2';% SEAC4RS ER2

% TCCON bias analysis
% -------------------
%sname = 'co2_lef_tower-insitu_1_allvalid-396magl';	% Park Falls	(LEF)
%sname = 'co2_sgp_surface-insitu_64_allvalid-60magl';	% Lamont	(SGP)
%sname = 'co2_asc_surface-flask_1_representative';	% Ascension	(ASC)
%sname = 'co2_crz_surface-flask_1_representative';	% Reunion-ish	(CRZ)
%sname = 'co2_sey_surface-flask_1_representative';	% Reunion-ish	(SEY)

% General bias analysis
% ---------------------
%sname = 'co2_brw_surface-insitu';			% Barrow		(BRW)
%sname = 'co2_spo_surface-insitu';			% South Pole		(SPO)

% *** WARNING: THESE NUMBERS DEPEND ON THE OBSPACK VERSION ***
% *** YOU SHOULD USE THE ABOVE TO SEARCH FOR THE SITE NAME INSTEAD ***

% Northern Hemisphere
% -------------------
%ic =  19; % Alert		(ALT)
%ic =  58; % Barrow		(BRW)

%ic = 165; % Niwot Ridge	(NWR)
%ic = 123; % Izana		(IZO)
%ic = 154; % Minamitorishima	(MNM)

% Tropics(ish)
% ------------
%ic = 210; % Mauna Loa		(MLO)
%ic = 289; % American Samoa	(SMO)

% Southern Hemisphere
% -------------------
%ic =  91; % Cape Point		(CPT)
%ic =  22; % Amsterdamn Isl.	(AMS)

%ic = 233; % Syowa		(SYO)
%ic = 300; % South Pole		(SPO)

% Canada
% ------
%ic = 119; % Nortwest Terr.	(INU)
%ic =  88; % British Columbia	(ESP)
%ic =  89; % Alberta		(EST)
%ic = 139; % Alberta		(LLB)
%ic =  37; % Saskatchewan	(BRA) 
%ic =  91; % Saskatchewan	(ETL)
%ic =  83; % Ontario		(EGB)
%ic =  93; % Ontario		(FSD)
%ic = 243; % Ontario		(TPD)
%ic =  73; % Quebec		(CPS)
%ic = 260; % Nova Scotia	(WSA)

% Europe
% ------
%ic = 146; % Ireland		(MHD)
%ic = 250; % England		(WAO)
%ic =  96; % Spain		(GIC)
%ic = 189; % France		(PUY)
%ic = 141; % Netherlands	(LUT)
%ic = 185; % Italy		(PRS)
%ic = 124; % Switzerland	(JFJ)
%ic = 125; % Switzerland	(JFJ)
%ic = 227; % Germany		(SSL)
%ic = 103; % Germany		(HEI)
%ic = 126; % Poland		(KAS)

% Oregon
% ------
%ic = 168; % Fir		(OFR)
%ic = 170; % Mary's Peak	(OMP)
%ic = 171; % Metolius		(OMT)
%ic = 172; % Burns		(ONG)
%ic = 174; % Silverton Tower	(OSI)
%ic = 176; % Walton		(OWA)
%ic = 178; % Yaquina Head	(OYQ)

% Other
% -----
%ic = 102; % Utah		(HDP)
%ic = 191; % Arizona		(RBA)
%ic = 220; % Storm Peak CO	(SPL)
%ic =  25; % Boulder Tower	(BAO)
%ic = 208; % Oklahoma		(SGP)
%ic = 156; % Martha's Vineyard	(MVY)

%ic = 262; % Yonagunijima	(YON)
%ic = 197; % Ryori		(RYO)

%ic = 264; % Ny-Alesund		(ZEP)
%ic = 180; % Pallas Cont.	(PAL)
%ic = 181; % Pallas Marine	(PAL)
%ic = 182; % Pallas Non-loc.	(PAL)


%==============================================================================%
dnout = [datenum(2009,01,01):1/24:datenum(2018,01,01)]';

%arctic_names  = {'co2_alt_surface-insitu_6_allvalid',		...
%                 'co2_brw_surface-insitu_1_allvalid',		...
%                 'co2_cby_surface-insitu_6_allvalid',		...
%                 'co2_inu_surface-insitu_6_allvalid',		...
%                 'co2_oli_surface-insitu_64_allvalid-10magl',	...	% Super noisy
%                 'co2_pal_surface-insitu_30_nonlocal',		...
%                 'co2_zep_surface-insitu_56_allvalid'};		% Messed up, use flask
arctic_names  = {'co2_alt_surface-insitu_6_allvalid',		...
                 'co2_brw_surface-insitu_1_allvalid',		...
                 'co2_pal_surface-insitu_30_nonlocal',		...
                 'co2_zep_surface-flask_1_representative'};
%tccon_names   = {'co2_lef_tower-insitu_1_allvalid-396magl',	...
%                 'co2_sgp_surface-insitu_64_allvalid-60magl',	...
%                 'co2_asc_surface-flask_1_representative',	...
%                 'co2_crz_surface-flask_1_representative',	...
%                 'co2_sey_surface-flask_1_representative'};

% Summit Greenland??? (SUM)
arctic_codes  = cell(numel(arctic_names), 1);
arctic_sitens = zeros(numel(arctic_names), 1);

arctic_biases = zeros(numel(dnout), numel(arctic_names));

figure;
hold on;
for ii = 1:numel(arctic_names)
  name = arctic_names{ii};

% Find site number
% ----------------
  for ic = 1:NSITES
    kk = strfind(cell_fobs{ic}, name);
    if (~isempty(kk)), break; end
  end

  if (isempty(kk))
    error(['Could not find matching site name: ', name, ' ...']);
  end

% Determine quality flags
% -----------------------
  isok = find(cell_qcflags{ic}(1,:) == '.' & cell_qcflags{ic}(2,:) == '.');
  isok = intersect(find(~isnan(cell_gasmod{ic})), isok);

% Average repeat obs
% ------------------
  [dnfix,iin,iout] = unique(cell_dnobs{ic}(isok), 'stable');

  omf   = cell_gasobs{ic}(isok) - cell_gasmod{ic}(isok);
  omfix = accumarray(iout, omf, [], @mean); 

  bias = -fit.fourier(dnfix, omfix, [1:numel(dnfix)]');
  arctic_biases(:,ii) = interp1(dnfix, bias, dnout);
  arctic_codes{ii}    = upper(name(5:7));
  arctic_sitens(ii)   = ic;

  plot(dnfix, bias, '.');
end
hold off;
legend(arctic_codes);
grid on; datetick('x', 'keeplimits');
