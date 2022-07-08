function arrayc = globcent(arrayn)

arrayc = arrayn;

% Average over lats
arrayc = 0.5*(arrayc(:,1:end-1) + arrayc(:,2:end));

% Average over lons
arrayc = [arrayc; arrayc(1,:)];
arrayc = 0.5*(arrayc(1:end-1,:) + arrayc(2:end,:));
