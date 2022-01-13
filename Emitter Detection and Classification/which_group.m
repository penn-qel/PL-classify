%Determines emitter classification using Rchi2, snr, width, allowed width
%range, and eccentricity (if applicable)

%Input: 
% Rchi2: goodness-of-fit parameter given by "BestFit" struct from either elliptical or symmetric Gaussian fit 
% snr: signal-to-noise ratio given by "BestFit" struct from either elliptical or symmetric Gaussian fit 
% width: given by "BestFit" struct from either elliptical or symmetric Gaussian fit 
% widthRange: microns. defined for given confocal setup and expected emission wavelength 
% eccentricity: given by "BestFit" struct from elliptical fit. Enables group D 

%Output: groupStr gives group classification such that: 
    %1 = group A
    %2 = group B
    %3 = group C
    %4 = group D

function groupStr = which_group(Rchi2,snr, width,widthRange,eccentricity)

groupStr = zeros(length(Rchi2),1);

%check if eccentricity is known according to elliptical fit. Otherwise
%assume symmetry and set eccentricity equal to zero. 
if isempty(eccentricity)
    eccentricity = zeros(length(Rchi2));
end 

%Check the width within a 1nm to width limits
tolerance = 0.001;
widthCheck = abs(width(width~=0)-widthRange(1))<tolerance | abs(width(width~=0)-widthRange(2))<tolerance  ;

%All groups are within the width limits 
for i = 1:length(eccentricity)

    %Group A: good Rchi2, bright 
    if widthCheck(i) == 0 && all(width(i,:)>widthRange(1)) && all(width(i,:)<widthRange(2)) && Rchi2(i)<1.5 && Rchi2(i)>0.8 && snr(i)>10   
        groupStr(i) = 1 ; 
    %Group B: good Rchi2, dim
    elseif widthCheck(i) == 0 && all(width(i,:)>widthRange(1)) && all(width(i,:)<widthRange(2)) && Rchi2(i)<1.5 && Rchi2(i)>0.8 && snr(i)<10 && snr(i)>2  
        groupStr(i) = 2;         
    %Group C: good Rchi2, dim, rail against width limits 
    elseif widthCheck(i) == 1 && snr(i)>2 && snr(i)<10 && Rchi2(i)<1.5 && Rchi2(i)>0.8   
        groupStr(i) = 3; 
   %Group D: Bright, good eccentricity, poor Rchi2
    elseif widthCheck(i) == 0 && all(width(i,:)>widthRange(1)) && all(width(i,:)<widthRange(2)) && snr(i)>10 && eccentricity(i)<=0.66 && eccentricity(i)>=0 && Rchi2(i)>1.5  
        groupStr(i) = 4; 
   end
end 
end 


