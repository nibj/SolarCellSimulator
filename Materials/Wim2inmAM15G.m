function val = Wim2inmAM15G(nmlambda,CSun)
% Input: wavelength (in nm)
%val = refractive index for wavelength lam0 (experimental data)
data = load('AM15G.mat');
inputnmlambda = data.data(:,1);
flux = CSun.*data.data(:,2);
val = interp1(inputnmlambda,flux,nmlambda);  % linear interpolation into the table
