%DIAG_AVGAREA  Performs a series of tests to make sure avgarea is working
%   correctly
%
%   1. HI -> LO conserves area averages
%   2. LO -> HI has no ringing
%   3. HI -> LO has no ringing
%   4. A step function stays between 0 and 1
%      (Equivalent to positivity of smoothing matrices)
%
%   Note: There is no reason the expect non-nested grids to pass the
%   LO -> HI -> LO test since the area averaging has to decide how to
%   distribute values across cells.  This can be accomplished by making one
%   transformation the Moore-Penrose pseudoinverse of the other.  However,
%   this transformation has negative values, which will then not preserve
%   positivity.  AFAIK, there is no pseudoinverse with all positive entries,
%   but cannot prove it.

% Author: Brad Weir
%==============================================================================%

addpath('/discover/nobackup/bweir/matlab/nanstats');

areahi = globarea(lathi, lonhi);
arealo = globarea(latlo, lonlo);

ooshi = ones(numel(lonhi), numel(lathi));
ooslo = avgarea(lathi, lonhi, ooshi, latlo, lonlo);
oosdn = avgarea(latlo, lonlo, ooslo, lathi, lonhi);

[LA, LO] = meshgrid(latlo, lonlo);
peakslo  = peaks(LO/36, LA/18);
peakshi  = avgarea(latlo, lonlo, peakslo, lathi, lonhi);
[peaksup, smx, smy, areaqq] = ...
           avgarea(lathi, lonhi, peakshi, latlo, lonlo);

totlo = sum(sum( arealo.*peakslo ));
tothi = sum(sum( areahi.*peakshi ));
totup = sum(sum( arealo.*peaksup ));
totqq = sum(sum( areaqq.*peaksup ));

nanslo = peakslo;
inan   = find(LA + LO < 0);
nanslo(inan) = NaN;

[nanshi, smxhi, smyhi] ...
       = avgarea(latlo, lonlo, nanslo, lathi, lonhi);
nansup = avgarea(lathi, lonhi, nanshi, latlo, lonlo);

oosnan = ones(numel(lonlo), numel(latlo));
oosnan(inan) = NaN;
oosndn = avgarea(latlo, lonlo, oosnan, lathi, lonhi);

% A. avg(peakshi) = avg(peakslo)
disp('Percent difference between HI & LO averages ...');
disp((totup - tothi)/tothi);
disp((totqq - tothi)/tothi);

% B. ooshi = 1 and oosdn = 1
disp('Mean and standard deviation of resampled all 1s ...');
disp([mean(ooslo(:)), std(ooslo(:))]);
disp([mean(oosdn(:)), std(oosdn(:))]);
disp([nanmean(oosndn(:)), nanstd(oosndn(:))]);

% C. arealo = smy' * areahi * smx
disp('Area averaging test (should be close to zero) ...');
checklo = smy'   * areahi * smx;
checkhi = smyhi' * arealo * smxhi;
disp( max(abs(checklo(:) - arealo(:)))/mean(arealo(:)) );
disp( max(abs(checkhi(:) - areahi(:)))/mean(areahi(:)) );

% D. LO -> HI -> LO = I
disp('RMS difference between LO -> HI -> LO and identity ...');
disp( sqrt(mean((peakslo(:) - peaksup(:)).^2)) );

disp('RMS difference between LO -> HI -> LO and identity (with NaNs) ...');
disp( sqrt(nanmean((nanslo(:) - nansup(:)).^2)) );
