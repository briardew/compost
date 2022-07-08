%MAKE_MAPLCONSTS  Save MAPL constants to file for easy reference

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/07/06	First version
%==============================================================================%

fname = ['maplconsts__', datestr(now,'yyyymmdd'), '.mat'];

% Constants defined in MAPL
% -------------------------
RADIUS  = 6371.0E3;              % m
GRAV    = 9.80665;               % m/s^2
MOLMDRY = 28.965;                % kg/Kmole
MOLMH2O = 18.015;                % kg/Kmole
MOLMO3  = 47.9982;               % kg/Kmole
RUNIV   = 8314.47;               % J/(Kmole K)
RDRY    = RUNIV/MOLMDRY;         % J/(kg K)
RVAP    = RUNIV/MOLMH2O;         % J/(kg K)
CPDRY   = 3.5*RDRY;              % J/(kg K)
KAPPA   = RDRY/CPDRY;            % (2.0/7.0)
EPS     = RVAP/RDRY-1.0;

MASSDRY = 5.1352;                % x 10^18 kg

% Molar masses/molecular weights
% ------------------------------
MOLMC   =  12.011;
MOLMO   =  15.9994;			% MOLMO3/3
MOLMH   =   1.0078;			% (MOLMH2O - MOLMO)/2, not so accurate (1.008)
MOLMCO  =  28.0104;			% MOLMC + 2*MOLMO
MOLMCO2 =  44.0098;			% MOLMC + 2*MOLMO
MOLMCH4 =  16.0422;			% MOLMC + 4*MOLMH
MOLMSF6 = 146.05;			% 32.05 + 6*19.00;

CTOCO2  = MOLMCO2/MOLMC;

disp(['Saving constants to ', fname, ' ...']);
save(fname);
