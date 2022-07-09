%CONSTANTS  Define atmospheric constants (following GEOS/MAPL values)

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/07/06	First version
%===============================================================================

RADIUS  = 6371.0E3;              % m
GRAV    = 9.80665;               % m/s^2
MASSDRY = 5.1352;                % x 10^18 kg

MOLMDRY = 28.965;                % kg/Kmole
MOLMH2O = 18.015;                % kg/Kmole
MOLMO3  = 47.9982;               % kg/Kmole

RUNIV   = 8314.47;               % J/(Kmole K)
RDRY    = RUNIV/MOLMDRY;         % J/(kg K)
RVAP    = RUNIV/MOLMH2O;         % J/(kg K)
CPDRY   = 3.5*RDRY;              % J/(kg K)
KAPPA   = RDRY/CPDRY;            % (2.0/7.0)
EPS     = MOLMH2O/MOLMDRY;

% for backup barometric equation
KAP1 = KAPPA + 1;
KAPR = 1/KAPPA;
BPOW = 8.31447*0.0065/(GRAV*0.0289644);
BCON = 288.15/0.0065;

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
