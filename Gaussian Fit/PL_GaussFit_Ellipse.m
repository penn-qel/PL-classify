%Elliptical Gaussian fit function that runs two fits. The first fit fixes the position of the emitter to
%estimate height, widths, and background. The second fit uses the results of the first
%fit as the input for height, widths, and background.

%Input: 
% emitters: struct containing CentroidXY and BoundingBoxXY for each detected object. 
% dataStruct: struct containing plScan, xCoords, yCoords.
% Vector iEvec: determines emitters to fit (ex: 1:length(emitters))
% widthLimit: 1x2 vector containing upper and lower width bounds
% SpotSize: scalar or 1x2 vector to estimate emitter width in microns. ex: 0.15
% bgval: initial background value for fit

%Output:
% BestFitNew: struct containing the results of the fit including position, widths,
% height, rotation in radians, signal-to-noise, reduced chi squared, background, eccentricity
% FitOptions: struct reflecting fit parameters
% Data: struct containing PL image and xCoords and yCoords

function [BestFitNew,FitOptions,Data] = PL_GaussFit_Ellipse(emitters,dataStruct,iEvec, widthLimit,BestFitSymmetric)

%Load PL image and binarize 
pl = dataStruct.data.plScan(:,:,1); 
res = mean(diff(dataStruct.data.xCoords)) ; 
thres = adaptthresh(mat2gray(pl), 0.3,'NeighborhoodSize',[51 51]);
B = imbinarize(mat2gray(pl),thres);
xCoords = dataStruct.data.xCoords;
yCoords = dataStruct.data.yCoords;

count=1;
%Loop for each region denoted by iEvec
for iE = iEvec

%Load region of interest
box = emitters(iE).BoundingBoxXY;
cropPl = imcrop(xCoords,yCoords,pl,box);
cropPlB = imcrop(xCoords,yCoords,B,box);

%define X and Y from box and set data
data.X = box(1):res:(box(1)+box(3)) + 4*res ; 
data.Y = box(2):res:(box(2)+box(4)) + 4*res ; 
data.Z = cropPl;

%Ensure bounding box does not exceed PL image bounds
[l, w] = size(cropPl) ; 
if length(data.X) ~= w 
    X = data.X;
    X = X(1:w);
    data.X = X;
end
if length(data.Y) ~= l
    Y = data.Y;
    Y = Y(1:l);
    data.Y = Y;
end 

%initialize emitter center, width, height from symmetric fit
clear startpoints startWidth bgval EstHts Posn Width Bkgnd Height 
startpoints = BestFitSymmetric(iE).Posn(:,1:2);
[r,c] = size(BestFitSymmetric(iE).Width(:)) ; 
StartWidth = 0.15.*ones(r,c); %[0.15, 0.15]; % Estimate for Gaussian st.dev (microns)
bgval = BestFitSymmetric(iE).Bkgnd;
EstHts = BestFitSymmetric(iE).Height(:) ;
FitOptions.StartHts = EstHts ; % starting amplitudes (scalar or vector)
np=BestFitSymmetric(iE).np; %emitters in the region

% Store starting parameters
FitOptions.npeaks = np ; % number of peaks to fit
FitOptions.StartPosns = startpoints; % starting positions (Nx2 array giving x,y)
FitOptions.StartWidths = StartWidth.*[1 1]; % starting width (scalar or vector)
FitOptions.StartBkgnd = bgval; % Starting background level (scalar)
FitOptions.StartTheta = 0 

% Parameter Limits:
FitOptions.PosnWindow = 0.3; % Peak positions are constrained to StartPosns±PosnWindow
FitOptions.LimHts = [0 , 5*max(EstHts,[],'all')]; % (min,max) limits on amplitude (sort in case any amplitudes are negative)
FitOptions.LimWidths = widthLimit; % (min,max) limits on width (all should be >0)
FitOptions.LimBkgnd = [1,3*bgval]
FitOptions.LimTheta = [0 pi]; 

% Fit options:
FitOptions.mode = 'FixedPosition'; % options are 'AllFree', 'FixedPosition', 'FixedWidth', 'FixedHeight','FixedTheta', or 'Arbitrary'
FitOptions.fittypes = true(np,4); % specifies freedom of (posn,width,ht) for each peak in 'Arbitrary' mode
FitOptions.CommonWidth = false; % use one width parameter common to all peaks (overrides fittypes above)
FitOptions.FixedBkgnd = false; % is background fixed or free?
FitOptions.CutoffDist = 3;  % To save computational time, gaussians are set to zero for distances greater than CutoffDist*std
FitOptions.PlotProgress = false; % Plots diagnostics of fmincon optimization
FitOptions.TrueHessian = false; % Whether to compute the true Hessian in order to calculate parameter confidence intervals
iE
[FitResult,BestFit,Errs] = Fit2dGaussiansEllipse(data,FitOptions);


%second fit
clearvars FitOptions FitResult Errs;

% Starting parameters (use best-fit results from previous constrained fit)
FitOptions.StartPosns = BestFit.Posn;
FitOptions.StartWidths = BestFit.Width;
FitOptions.StartHts = BestFit.Height;
FitOptions.StartBkgnd = BestFit.Bkgnd;

% Parameter Limits:
FitOptions.PosnWindow = 0.3; % Peak positions are constrained to StartPosns±PosnWindow
FitOptions.LimHts = [0,5*max(EstHts,[],'all')]; % (min,max) limits on amplitude (sort in case any amplitudes are negative)
FitOptions.LimWidths = widthLimit; % (min,max) limits on width (all should be >0)
FitOptions.LimBkgnd = [1,3*bgval];
FitOptions.LimTheta = [0 pi] ; %allows for rotation along axis

% Fit options:
FitOptions.mode = 'Arbitrary'; % options are 'AllFree', 'FixedPosition', 'FixedWidth', 'FixedHeight', 'FixedTheta', or 'Arbitrary'
FitOptions.fittypes = true(np,4); % specifies freedom of (posn,width,ht,theta) for each peak in 'Arbitrary' mode
% Fix a few parameters
%FixedPeakPosns = [11 19 20 24]; % constrain position for these peaks
%FitOptions.fittypes(FixedPeakPosns,1) = false;

FitOptions.CommonWidth = false;  % use one width parameter common to all peaks (overrides fittypes above)
FitOptions.FixedBkgnd = false; % is background fixed or free?
FitOptions.CutoffDist = 3; % To save computational time, gaussians are set to zero for distances greater than CutoffDist*std
FitOptions.PlotProgress = false; % Plots diagnostics of fmincon optimization
FitOptions.TrueHessian = false; % Whether to compute the true Hessian in order to calculate parameter confidence intervals

[FitResult,BestFit,Errs] = Fit2dGaussiansEllipse(data,FitOptions);

[chiSquareSum, ~,nPoints] = chi_square(data, BestFit);

%store results of fit --  including width, Rchi2, snr, eccentricity
BestFitNew(count).np = np;
BestFitNew(count).Posn = BestFit.Posn;
BestFitNew(count).TruePosn = BestFit.TruePosn;
BestFitNew(count).Width = BestFit.Width;
BestFitNew(count).Height = BestFit.Height;
BestFitNew(count).Theta = BestFit.theta
BestFitNew(count).Drift = BestFit.Drift;
BestFitNew(count).Bkgnd = BestFit.Bkgnd;
BestFitNew(count).Errs = Errs;
BestFitNew(count).ModelFunction = BestFit.ModelFunction;
BestFitNew(count).data = data; 
BestFitNew(count).Rchi2 = chiSquareSum./(nPoints-(6*np+1)); %Rchi2 for elliptical fit is not relevant to emitter classification. DOF: np*(height, Xwidth, Ywidth, xPos, yPos, theta) + background
BestFitNew(count).nPoints = l*w;
BestFitNew(count).snr = (BestFit.Height)./sqrt(BestFit.Bkgnd);

%Calculate eccentricity
Width = cat(1,BestFit.Width); 
eccentricity = zeros(1,length(Width(:,1)));
for i = 1:length(Width(:,1))
    eccentricity(i) = sqrt(1-((min(Width(i,:))).^2/(max(Width(i,:))).^2));
end 
BestFitNew(count).eccentricity = eccentricity; 


%Estimate background points according to maximum width from fit
maxWidth = max(BestFit.Width,[],'all'); 
anyCenter = BestFit.Posn(1,:); 
[X, Y] = meshgrid(data.X,data.Y) ;
radius = sqrt((X-anyCenter(1)).^2  + (Y-anyCenter(2)).^2); 
cropPlBinary  = radius<(2*maxWidth); 
data.cropPlB = cropPlBinary; 
Data(count) = data; %store data for output
count = count+1;

clear data    
end


end 

