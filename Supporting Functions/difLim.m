%Function calculates the diffraction limit

%Input: numerical aperture (NA) of the objective, and wavelength (wl) of exciting light (nm)

%Output: diffraction limit, Airy-disk radius, expected gaussian width,
%full width half maximum 

function [diffractionLimit, airyRadius,sigma,FWHM]= difLim(NA,wl)

wl = wl.*(10^-9); 
diffractionLimit = wl./(2.*NA).*(10^6) ; %um

sigma = 0.21.*(1/NA).*wl.*(10^6) ; %um

airyRadius = 1.22.*wl.*(1./(2.*NA)).*(10^6); %um

FWHM = 2.35.*sigma;
end 
