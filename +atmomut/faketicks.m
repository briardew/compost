%FAKETICKS  Hacky utility to throw on ticks for fake map plots

% Author: Brad Weir
%==============================================================================%

function faketicks(xtickloc, ytickloc)

lh = ishold;
hold on;

for nn = 1:numel(xtickloc)
  plot(xtickloc(nn)*[1 1], ylim, ':', 'color', [0.75 0.75 0.75]);
end

for nn = 1:numel(ytickloc)
  plot(xlim, ytickloc(nn)*[1 1], ':', 'color', [0.75 0.75 0.75]);
end

if (lh), hold on; end
