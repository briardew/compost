%FIT.FOURIER  Compute Fourier series fit to obs - model residuals
%
%   I don't know

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/04/04	New version based off of noaa make_bias script
%
% TODO:
%==============================================================================%
function bias = fourier(dnobs, omfin, isok)

% FFT needs evenly spaced data:
%   * For now, interpolate omfs to hourly with a 365.25 day year
%   * Probably should average the obs to 3hrs and interpolate to fill
npd = 24;							% Samples per day
dt  = 1/(npd*365.25);						% Timestep
dnfit = [dnobs(isok(1)):1/npd:dnobs(isok(end))];		% Time grid to do FFT on
NFIT  = numel(dnfit);						% Number of points to fit to

omf   = interp1(dnobs(isok), omfin(isok), dnfit);

% Fourier analysis
LL = 2*floor(numel(omf)/2);
k1 = (1/dt)*[0:(LL/2)]'/LL;					% 1-sided frequency

ff = fft(omf);
p2 = abs(ff/LL);						% 2-sided spectrum
p1 = p2(1:LL/2+1);						% 1-sided spectrum
p1(2:end-1) = 2*p1(2:end-1);

% Long-term secular bias
gg = ff;
%gg(7:end-6) = 0;
gg(22:end-21) = 0;
gg(1:6) = 0;
gg(end-5:end) = 0;
ltb = real(ifft(gg));

bias = interp1(dnfit, ltb, dnobs) + mean(omfin(isok));
