function [FitResult,BestFit,Errs,bfpars] = Fit2dGaussiansEllipse(data,fitopts)



%Function to fit 2d data to a set of (not symmetric) Gaussian peaks
%
% The structure data must contain the field
% 'Z' (2d array): data to fit
% and the following fields are optional
% 'X', 'Y' (vectors or 2d arrays): 2d coordinates corresponding to 'Z' If
%   omitted, pixel indices are used as coordinates for 'Z'
% 'Weights' (2d array): gives weights to use for each point, e.g. inverse
%   of variance for statistical uncertainty
% 'Mask' (2d binary array): used to exclude points from fit (set weights to zero)

% The structure fitopts contains information to constrain the fit, and
% can contain the fields:
%
% 'npeaks' (scalar): number of Gaussian peaks to fit (note in principle these
%    could be dips with negative amplitude)
% 'StartPosns' (Nx2 array): initial positions (pixels or coordinates)
% 'StartHts' (scalar or Nx1 vector): initial amplitude for each Gaussian

% 'StartWidths' (scalar or Nx2 vector): initial x and y width for each Gaussian 

% 'StartBkgnd' (scalar): starting background value for fit
% 'PosnWindow' (scalar or Nx1 vector): size of window to constrain position
%   for each peak about 'startposns' in both x and y (square box)
% 'LimHts' (1x2 or Nx2 array): [min, max] limits for amplitude
% 'LimWidths' (1x2 or Nx2 array): [min, max] limits for widths
% 'LimBkgnd' (1x2 vector): [min, max] limits for background
% 'mode' (string): defines mode of fit following
%       'AllFree': All parameters free (default)
%       'FixedPosition': All positions fixed
%       'FixedWidth': All widths fixed
%       'Arbitrary': Freedom of parameters is determined by additional
%           field initstruc.fittypes 
% 'fittypes' (Nx4 logical array): Needed if 'mode'=='Arbitrary' to specify
%   freedom of each parameter, with rows giving [posn,width,ht,theta] and values
%   0==fixed, 1==free
% 'CommonWidth' (T/F): if true there is one common width parameter for
%   all peaks, otherwise all are independently defined.
% 'CommonDrift' (T/F): if true, the peak positions are allowed to shift
%   together by a common drift vector (delta_x,delta_y), even if the
%   relative peak positions are fixed.  If all positions are free to vary
%   and CommonDrift==T, the drift vector fully determines the shift in peak
%   #1 while others remain free to vary about the shifted relative
%   positions. Resulting position uncertainties are calculated using propagation that
%   assumes no correlation between drift and relative positions.
% 'StartDrift' (1x2 vector): Starting value for drift vector (default = [0,0])  
% 'DriftWindow' (scalar or 1x4 array): Determines bounds for drift vector
%   if "CommonDrift"==T.  A scalar value (a) bounds the vector to a box of
%   width [-a a] in both x and y about "StartDrift".  Alternatively,
%   arbitrary bounds can be specified by a 1x4 vector giving
%   [xmin,xmax,ymin,ymax].  
% 'FixedBkgnd' (T/F): Background level fixed if true (default ==F)
% 'CutoffDist' (scalar): In the model, Gaussians are set to zero for
%   distances greater than CutoffDist*std (default ==Inf)
% 'PlotProgress' (T/F): Whether to plot diagnostics of fitting progress
%   (default = F)
% 'TrueHessian' (T/F): Whether to compute the true Hessian in order to
%   estimate parameter uncertainties or use that from fmincon for speed
%   (default = F).  Requires hessian.m from DERIVEST suite.
% 'TrustWeights' (T/F): Whether to trust that the weights accurately
% 	reflect the variance in the data when calculating parameter confidence
% 	intervals.  Otherwise rescales the weights such that
% 	reduced-chi-squared equals one before calculating the covariance
% 	matrix.  Only used for weighted fits.  (Default = T).
%
% The optional bkgndopts structure may contain a suitable model to fit the background variation
% (e.g. a polynomial or sum of broad gaussians) with the fields
% model (str), startpars (vec), minpars (vec), maxpars (vec).  The bkgnd model
% should use parameters named bi (for integer i) to avoid confusion with pars
% in the peak-fitting function and the corresponding initialisation parameters should be
% listed in ascii dictionary order for these parameters (e.g.
% b1,b10,b11,...,b19,b2,b20,...).
%
% Outputs:
%
% FitResult: contains outputs of fmincon (fields are named after output variables)
%
% BestFit: Structure containing best-fit parameters
%   'np'
%   'Posn' (Results from fit not including drift vector, if any)
%   'TruePosn' (Including drift vector, equal to Posn if drift=0)
%   'Width'
%   'Height'
%   'Drift'
%   'Bkgnd'
%   'ModelFunction': Function handle that calculates best-fit model given (X,Y)  
    %LEAH Added theta for rotation 
%
% Errs: Structure containing parameter uncertainties at 68% confidence, with same fields as
% 'BestFit' (without Modelfunction).  These are calculated from the hessian
% of the objective function only if this output is requested.
%% Parse inputs

% Parse input data
Zdata = data.Z;
if isfield(data,'X')
    if all(size(data.X)==size(Zdata))
%         Xdata = data.X;
    elseif isvector(data.X) && length(data.X)==size(Zdata,2)
        data.X = repmat(data.X(:)',size(Zdata,1),1); % make into array of same size as Zdata
    else
        error('X coordinates must either be a vector or an array with the same dimensions as Z data');
    end
else
    data.X = repmat(1:size(Zdata,2),size(Zdata,1),1); % use indices if no coordinates are provided
end
if isfield(data,'Y')
    if all(size(data.Y)==size(Zdata))
%         Ydata = data.Y;
    elseif isvector(data.Y) && length(data.Y)==size(Zdata,1)
        data.Y = repmat(data.Y(:),1,size(Zdata,2)); % make into array of same size as Zdata
    else
        error('Y coordinates must either be a vector or an array with the same dimensions as Z data');
    end
else
    data.Y = repmat((1:size(Zdata,1))',1,size(Zdata,2));
end
if isfield(data,'Weights')
    if any(size(data.Weights)~=size(Zdata))
        error('Weights must be an array of the same size as Z');
    end
    Weighted = true;
else
    data.Weights = ones(size(Zdata));
    Weighted = false;
end
if isfield(data,'Mask')
    validateattributes(data.Mask,{'logical','numeric'},{'binary','size',size(Zdata)});
%     if any(size(data.Mask)~=size(Zdata))
%         error('Mask must be an array of the same size as Z');
%     end
    data.Weights = data.Weights.*data.Mask; % combine with Weights
end
if any(isnan(Zdata(find(data.Weights))))
    error('Data contains NaN values in unmasked region.');
else % any NaNs occur in masked region, so set to zero to avoid problems with computing weighted SSE
    data.Z(isnan(Zdata)) = 0; 
end


xstep = mean(diff(data.X(1,:)));
ystep = mean(diff(data.Y(:,1)));
meanstep = mean([xstep,ystep]);

% Parse fit structure to extract given parameters and apply defaults
if ~isfield(fitopts,'StartPosns') % This field is required so check if it is present
    error('Fit options structure must contain the field "startposns"');
else % record number of peaks
    np = size(fitopts.StartPosns,1);
end

% Note I need to use addParamValue() rather than addParameter() since I am
% still running R2012b on my laptop

%Leah: changed to addParameter
parse_fitopts = inputParser;
parse_fitopts.KeepUnmatched = true; % allows other parameters to be included that are not needed by this function
parse_fitopts.FunctionName = 'Fit2dGaussiansWidth';
addParameter(parse_fitopts,'npeaks',...
    np,...
    @(x)validateattributes(x,{'numeric'},{'scalar','integer'}));
if isfield(fitopts,'np') && fitopts.np ~=np % consistency check
    error('Number of rows in "StartPosns" does not match field "np"');
end

addParameter(parse_fitopts,'StartPosns',...
    NaN(np,2),... % default value, but should not happen since checked above to be sure was given
    @(x)validateattributes(x,{'numeric'},{'2d','ncols',2}));
addParameter(parse_fitopts,'StartHts',...
    mean(Zdata(:))*ones(np,1),... % default value = mean of response data
    @(x) assert(isnumeric(x) && (isscalar(x) || length(x)==np),...
    'Parameter must be a scalar or a vector of length "npeaks"'));
addParameter(parse_fitopts,'StartTheta',...
    0,... % default no angular shift
    @(x) assert(isnumeric(x) && (isscalar(x) || length(x)==np),...
    'Parameter must be a scalar or a vector of length "npeaks"'));
addParameter(parse_fitopts,'StartWidths',...
    10*meanstep,... % default width = 10 * mean step in X and y direction assuming reasonable (?) sampling of Gaussian peaks.  Obviously better to provide this value
    @(x) assert(isnumeric(x) && (isscalar(x) || length(x)==np || length(x) == 2*np),...
    'Parameter must be a scalar or a vector of length "npeaks" or 2*length "npeaks"'));
addParameter(parse_fitopts,'StartBkgnd',...
    0,...
    @(x) validateattributes(x,{'numeric'},{'nonempty','scalar'}));
addParameter(parse_fitopts,'PosnWindow',...
    Inf,... % default value
    @(x) assert(isnumeric(x) && (isscalar(x) || length(x)==np && all(x)>0),...
    'Parameter must be a scalar or a vector of length "npeaks"'));
addParameter(parse_fitopts,'LimHts',...
    [-Inf,Inf],... % default
     @(x) assert(isnumeric(x) && (all(size(x)==[1,2]) || all(size(x)==[np,2])),...
    'Parameter must be size 1x2 or Nx2 where N=number of peaks'));
addParameter(parse_fitopts,'LimWidths',...
    [0,Inf],... % default
     @(x) assert(isnumeric(x) && (all(size(x)==[1,2]) || all(size(x)==[np,2])),...
    'Parameter must be size 1x2 or Nx2 where N=number of peaks'));
addParameter(parse_fitopts,'LimBkgnd',...
    [-Inf,Inf],...
     @(x) validateattributes(x,{'numeric'},{'numel',2}));
addParameter(parse_fitopts,'mode',...
    'AllFree',... % default value
    @(x) any(strcmp(x,{'AllFree','FixedPosition','FixedWidth','FixedHeight','FixedTheta','Arbitrary'})));
addParameter(parse_fitopts,'fittypes',...
    true(np,4),... % default is all free (if not given)
    @(x) validateattributes(x,{'logical','numeric'},{'binary','size',[np,4]}));
addParameter(parse_fitopts,'CommonWidth',...
    false,... % default
    @(x) validateattributes(x,{'logical','numeric'},{'binary'}));
addParameter(parse_fitopts,'FixedBkgnd',...
    false,... % default
    @(x) validateattributes(x,{'logical','numeric'},{'binary'}));
addParameter(parse_fitopts,'CommonDrift',...
    false,... % default
    @(x) validateattributes(x,{'logical','numeric'},{'binary'}));
addParameter(parse_fitopts,'StartDrift',...
    [0,0],... % default value
    @(x) assert(isnumeric(x) && (isscalar(x) || all(size(x)==[1,2])),...
    'Parameter must be a scalar or a vector with two components'));
addParameter(parse_fitopts,'DriftWindow',...
    Inf,... % default value
    @(x) assert(isnumeric(x) && (isscalar(x) || length(x)==4),...
    'Parameter must be a scalar or a vector with four components'));
addParameter(parse_fitopts,'CutoffDist',...
    Inf,... % default
    @(x) validateattributes(x,{'numeric'},{'positive'}));
% addParamValue(parse_fitopts,'Weighted',...
%     false,... % default
%     @(x) validateattributes(x,{'logical','numeric'},{'binary'}));
addParameter(parse_fitopts,'PlotProgress',...
    false,... % default
    @(x) validateattributes(x,{'logical','numeric'},{'binary'}));
addParameter(parse_fitopts,'TrueHessian',...
    false,... % default
    @(x) validateattributes(x,{'logical','numeric'},{'binary'}));
addParameter(parse_fitopts,'TrustWeights',...
    true,... % default
    @(x) validateattributes(x,{'logical','numeric'},{'binary'}));

parse(parse_fitopts,fitopts);
fitpars = parse_fitopts.Results;

if nargin>2
    error('Arbitrary backgrounds are not yet supported by this function.');
end

%% Expand all arrays to have np rows (can be provided as a single value or pair to apply to all peaks)

if length(fitpars.StartWidths) == 1
    fitpars.StartWidths = repmat(fitpars.StartWidths,np,2);
end
if length(fitpars.StartHts) == 1
    fitpars.StartHts = repmat(fitpars.StartHts,np,1);
end
if size(fitpars.LimWidths,1) == 1
    fitpars.LimWidths = repmat(fitpars.LimWidths,np,1);
end
if size(fitpars.LimHts,1) == 1
    fitpars.LimHts = repmat(fitpars.LimHts,np,1);
end
if size(fitpars.PosnWindow,1) == 1
    fitpars.PosnWindow = repmat(fitpars.PosnWindow,np,1);
end
if length(fitpars.StartTheta) == 1
    fitpars.StartTheta = repmat(fitpars.StartTheta,np,1);
end

%% Parse fitting options to determine free and fixed parameters

switch fitpars.mode
    case 'AllFree'
        fitpars.fittypes = true(np,4); % override any input values and set all parameters to be free
    case 'FixedPosition'
        fitpars.fittypes(:,1) = false; % override any input values to fix positions
    case 'FixedWidth'
        fitpars.fittypes(:,2) = false; % override any input values to fix widths
    case 'FixedHeight'
        fitpars.fittypes(:,3) = false; % override any input values to fix heights
    case 'FixedTheta'
    fitpars.fittypes(:,4) = false; % override any input values to fix heights
    case 'Arbitrary'
        if ~isfield(fitopts,'fittypes') % check that fittypes was provided by user (otherwise was set to all free by default)
            error('To use "Arbitrary" fit mode, "fittypes" array must be provided in the fit options.');
        end
    otherwise
        error('Unrecognized fit mode parameter: "%s"',fitpars.mode);
end

if fitpars.CommonWidth % in this case the starting/limit values and fittype for peak 1 applies to all peaks
    fitpars.StartWidths = repmat(fitpars.StartWidths(1),np,1);
    fitpars.LimWidths = repmat(fitpars.LimWidths(1,:),np,1);
    fitpars.fittypes(2:end,2) = false; % only peak 1 is allowed to be free (value set by user)
end

if fitpars.CommonDrift % All peak positions are allowed to shift by a common drift vector
    if isscalar(fitpars.DriftWindow)
        fitpars.DriftWindow = fitpars.DriftWindow*[-1 1; -1 1]; % [xmin, xmax; ymin, ymax]
    else % provided as row vector [xmin,xmax,ymin,ymax]
        fitpars.DriftWindow = reshape(fitpars.DriftWindow,[2,2])';
    end
    
    if all(fitpars.fittypes(:,1)) % if all positions are free, need to constrain one to include drift vector
        fitpars.fittypes(1,1) = false;
        constrain_pk1_for_drift = true; 
        % Note that this will have an effect on the allowed bounds for
        % relative positions of other peaks.  Fixing this position
        % effectively reduces the allowed relative shifts between peak 1 and all
        % other peaks by the value of posn_windows(1).  But it would not be
        % correct to simply add that value to the posn_windows of other
        % peaks since that would allow greater relative shifts between
        % them.  In practice it is better for the user to determine which
        % peak to use as the reference and explictly fix it in the peak,
        % but this seems to be a reasonable default action.
    else
        constrain_pk1_for_drift = false;
%         drift_ref = find(~fitpars.fittypes(:,1),'first'); % drift referenced to first fixed peak
    end
end
    
%% Set up problem for fmincon

% Extract free parameters
freeposns = fitpars.StartPosns(fitpars.fittypes(:,1),:);
freewidths = fitpars.StartWidths(fitpars.fittypes(:,2),:);
freehts = fitpars.StartHts(fitpars.fittypes(:,3));
freetheta = fitpars.StartTheta(fitpars.fittypes(:,4));

if fitpars.CommonDrift
    freedrift = fitpars.StartDrift(:);
else
    freedrift = [];
end
freebg = fitpars.StartBkgnd(~fitpars.FixedBkgnd);

posn_windows = fitpars.PosnWindow(fitpars.fittypes(:,1));
if ~isempty(freeposns)
    lim_posns = repmat(freeposns(:),[1,2])+repmat(posn_windows*[-1 1],[2,1]) ;  % [lbx ubx; lby uby]
else
    lim_posns = double.empty(0,2); % 0x2 empty array
end

lim_windows = fitpars.LimWidths(fitpars.fittypes(:,2),:);

if ~isempty(freewidths)
    lim_widths = repmat(lim_windows,[2,1]) ; % [lbx ubx; lby uby]
else
    lim_widths= double.empty(0,2); % 0x2 empty array
end

if ~isempty(freetheta)
    lim_theta = ones(np,2).*[0 2*pi] ; 
end

%lim_widths


lim_hts = fitpars.LimHts(fitpars.fittypes(:,3),:);
if ~isempty(freedrift)
    lim_drift = repmat(fitpars.StartDrift(:),[1 2])+fitpars.DriftWindow;
else
    lim_drift = double.empty(0,2); % 0x2 empty array
end
lim_bg = fitpars.LimBkgnd(~fitpars.FixedBkgnd,:);

% Total number of free parameters
npars = length(freeposns)+length(freewidths)+length(freehts)+length(freedrift)+length(freebg) + length(freetheta);

x0 = [freeposns(:);... % [x;y] combined in column vector
    freewidths(:);... %[wX ; wY] combined in column vector
    freehts;... %h
    freetheta;...
    freedrift;... % drift vector [x;y]
    freebg]; % background parameters

% Create Nx2 vector giving [min,max] for each parameter in x0
all_bounds = [lim_posns;... % [x;y]
    lim_widths;... %[x;y]
    lim_hts;...
    lim_theta,...
    lim_drift;... %[x;y]
    lim_bg];

lb = all_bounds(:,1);
ub = all_bounds(:,2);

problem.solver = 'fmincon';
% problem.options = optimoptions('fmincon','Algorithm','interior-point');
problem.options = optimset('fmincon');
problem.options.Algorithm = 'interior-point';
problem.options.Display = 'final-detailed';
problem.options.MaxFunEvals = 1e4; 
problem.options.MaxIter = 1000;
if fitpars.PlotProgress
    problem.options.PlotFcns = {@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotstepsize};
end % else use default = []
problem.options.TolX = 1e-5; % might be a good idea to scale parameters such that the all have the same target precision
problem.options.TolFun = 1e-8; % want to hit parameter tolerance rather than function tolerance
problem.options.UseParallel = 'always';
TypicalX = x0;
% TypicalX(x0==0) = FitStruc.scales(problem.x0==0)*fitopts.targetprecision;
TypicalX(x0==0) = problem.options.TolX; % these cannot be zero, so just set to TolX
problem.options.TypicalX = TypicalX;

problem.x0 = x0;
problem.lb = lb;
problem.ub = ub;

problem.objective = @(x) ObjectiveFunction(x,data,fitpars);


%% Use fmincon to fit

% tic;
[FitResult.x,FitResult.SSE,FitResult.exitflag,FitResult.output,FitResult.lambda,FitResult.grad,FitResult.hessian] = fmincon(problem);
% toc;
% FitResult.problem = problem; % Store problem, e.g. to use ComputeConfidenceIntervals later on

%% Collate outputs


bfpars = Create2dGaussianModelWidth(FitResult.x,fitpars);

BestFit.np = np;
BestFit.Posn = [bfpars.x,bfpars.y]; % These positions DO NOT include the drift vector
BestFit.TruePosn = BestFit.Posn+repmat(bfpars.drift,[np,1]); % True positions including the drift vector
BestFit.Width = bfpars.w;
BestFit.Height = bfpars.a;
BestFit.Drift = bfpars.drift;
BestFit.Bkgnd = bfpars.bgpars;
BestFit.theta = bfpars.theta; 
BestFit.ModelFunction = @(X,Y) Sum2dGaussians(bfpars,X,Y); % function handle that will calculate model function using best-fit parameters

a = (cos(bfpars.theta).^2)./(2.*bfpars.w(1).^2) + (sin(bfpars.theta).^2)./(2.*bfpars.w(2).^2) ; 
c = (sin(bfpars.theta).^2)./(2.*bfpars.w(1).^2) + (cos(bfpars.theta).^2)./(2.*bfpars.w(2).^2) ; 
    
BestFit.eccentricity = c/a; 

%% Calculate uncertainties

if nargout>2 && FitResult.exitflag>0 % only calculate if requested and fit was successful

    if fitpars.TrueHessian % compute actual Hessian (O(6*n^2) operations)
        hess = hessian(problem.objective,FitResult.x);
        FitResult.TrueHessian = hess; % record for later use
    else
        hess = FitResult.hessian;
    end
    dof = sum(data.Weights(:)>0) - npars;
    
    if Weighted && fitpars.TrustWeights % if the weighted SSE is calculated (returning chi-squared), and we "trust" the weights, then this is probably more accurate
        ErrMat = inv(hess./2);
    else
        rescale = dof/FitResult.SSE; % rescaling factor such that chisquare = dof = Nvals-npars ... not ideal but should give the right order of magnitude
        ErrMat = inv(rescale*hess./2);
    end
    FitResult.sigma_x = sqrt(diag(ErrMat));
    
    err_struct = Create2dGaussianModelWidth(FitResult.sigma_x,fitpars); % will collate signa_x into structure (values of fixed parameters are starting positions, but get replaced by NaN below)
    
    Errs.Posn = [err_struct.x,err_struct.y];
    Errs.Posn(~fitpars.fittypes(:,1))= NaN; % set confidence intervals for fixed parameters to NaN
    if fitpars.CommonDrift % need to account for uncertainty in both positions and drift vector
        Errs.Posn = hypot(Errs.Posn,repmat(err_struct.drift,[np,1])); % sum uncertainties in quadrature, assuming no correlations
        if constrain_pk1_for_drift
            Errs.Posn(1,:) = err_struct.drift; % This peak position was intended to be free, but fixed to serve as a reference for the drift vector
        end
        Errs.Drift = err_struct.drift;
    else
        Errs.Drift = NaN(1,2);
    end
    
    Errs.Width = err_struct.w;
    Errs.Width(~fitpars.fittypes(:,2))= NaN;
    
    Errs.Height = err_struct.a;
    Errs.Height(~fitpars.fittypes(:,3))= NaN;
    
    if fitpars.FixedBkgnd
        Errs.Bkgnd = NaN;
    else
        Errs.Bkgnd = err_struct.bgpars;
    end
    
else
    Errs = struct([]);
end


%% Objective function

function SSE = ObjectiveFunction(x,data,fitoptions)
% Function to minimize -- returns weighted sum of squared error between
% model and fit

model = Create2dGaussianModelWidth(x,fitoptions);

SSE = WeightedSSE(Sum2dGaussians(model,data.X,data.Y),data.Z,data.Weights);


%% Function to create input for Sum2dGaussians from vector of free parameters and fit options

function model = Create2dGaussianModelWidth(x,fitoptions)

fittypes = fitoptions.fittypes;


npars = sum(fittypes(:)) + sum(fittypes(:,1)) + sum(fittypes(:,2)) + 2*fitoptions.CommonDrift + ~fitoptions.FixedBkgnd; % expected number of free parameters
if length(x)~=npars
    error('Length of input vector does not match expected number of free parameters');
end

lastix = 0; % last index in x vector

allpars = NaN(fitoptions.npeaks,6); % will hold all parameters in array giving [x,y,wX,wY,a, theta] as columns
% Parse positions
freexy = fittypes(:,1); % indices of peaks with free positions
nfree = sum(freexy);
if any(freexy)
    allpars(freexy,1) = x(lastix+(1:nfree));
    lastix = lastix+nfree;
    allpars(~freexy,1) = fitoptions.StartPosns(~freexy,1);
    allpars(freexy,2) = x(lastix+(1:nfree));
    lastix = lastix+nfree;
    allpars(~freexy,2) = fitoptions.StartPosns(~freexy,2);
else
    allpars(:,[1 2]) = fitoptions.StartPosns;
end
% Parse widths
freew = fittypes(:,2); % indices of free widths
nfree = sum(freew);
if any(freew)
    if fitoptions.CommonWidth
        if nfree>1
            error('Error parsing parameter vector. Expected only one free width parameter in "CommonWidth" mode.');
        end
        allpars(:,3:4) = x(lastix+1); % all widths are set to the same value
        lastix = lastix+nfree; % should be lastix+1
    else
        allpars(freew,3) = x(lastix+(1:nfree));
        lastix = lastix + nfree; 
        allpars(~freew,3) = fitoptions.StartWidths(~freew,1);
        
        %LEAH ADDED
         allpars(freew,4) = x(lastix+(1:nfree));
        lastix = lastix + nfree; 
        allpars(~freew,4) = fitoptions.StartWidths(~freew,1);
    end
else
    allpars(:,3:4) = fitoptions.StartWidths;
end

% Parse amplitudes  
freeht = fittypes(:,3);
nfree = sum(freeht);
if any(freeht)
    allpars(freeht,5) = x(lastix+(1:nfree));
    lastix = lastix+nfree;
    allpars(~freeht,5) = fitoptions.StartHts(~freeht);
else
    allpars(:,5) = fitoptions.StartHts;
end

%LEAH ADDED
freetheta = fittypes(:,4) ; 
nfree = sum(freetheta) ; 
if any(freetheta)
    allpars(freetheta,6) = x(lastix+(1:nfree));
    lastix = lastix+nfree;
    allpars(~freetheta,6) = fitoptions.StartTheta(~freetheta);
else
    allpars(:,6) = fitoptions.StartTheta;
end


% Drift vector
if fitoptions.CommonDrift
    drift_vec = x(lastix+(1:2))'; % [driftx,drifty]
    lastix = lastix+2;
else
    drift_vec = [0,0];
end

% Include background
if ~fitoptions.FixedBkgnd
    bgval = x(lastix+1); % should be background value
    lastix = lastix+1;
else
    bgval = fitoptions.StartBkgnd;
end

if lastix ~= npars
    error('Error parsing input vector of free parameters'); % should never happen but just in case
end

% Now construct model using these parameters
model.N = fitoptions.npeaks;
model.x = allpars(:,1);
model.y = allpars(:,2);
model.w = allpars(:,3:4);
model.a = allpars(:,5);
model.theta = allpars(:,6) ; 
model.drift = drift_vec;
model.bgpars = bgval;
model.bgfn = @(X,Y) bgval*ones(size(X)); % function that returns background value (planned future improvement)    
model.cutoff = fitoptions.CutoffDist;


%% Function that calculates sum of 2d Gaussians

function Z = Sum2dGaussians(model,X,Y)
% Calculates sum of 2d Gaussians with parameters defined in 'model' at
% points given in X, Y (scalars, vectors, or matrices)

% model needs to contain:
%   'N': number of gaussians
%   'x': Nx1 vector of x positions
%   'y': Nx1 vector of y positions
%   'w': Nx2 [Wx;Wy] vector of Gaussian standard deviation 
%   'a': Nx1 vector of Gaussian amplitudes
%   'drift': 2x1 vector giving common drift of all positions in [x,y]
%   'cutoff': Values a distance greater than cutoff*w from peak N at will be set to zero to save computational time 
%   'bgfn': Function handle that returns background given (X,Y)

np = model.N;
cutoff = model.cutoff; % used to define maximum radius (units of s.d.) beyond which value is set ==0 for speed
Z = model.bgfn(X,Y); % Initialize with background values at points X,Y
for pp=1:np
    xp = model.x(pp)+model.drift(1);
    yp = model.y(pp)+model.drift(2);
    w = model.w(pp,:);
    wX = w(1);
    wY = w(2);
    A = model.a(pp);
    theta = model.theta(pp);
    dist2 = (X-xp).^2 + (Y-yp).^2; % array of same size as (X,Y) giving square of distance from (xp,yp)
    if isfinite(cutoff)
        PtsInRange = (dist2 <= min((cutoff.*w).^2)); % points located within min of cutoff*w(x,y) of (xp,yp)
    else
        PtsInRange = true(size(X));
    end
    ThisPeak = zeros(size(X));
    
    xPoints = X(PtsInRange);
    yPoints = Y(PtsInRange);
    
    
    a = (cos(theta).^2)./(2.*wX.^2) + (sin(theta).^2)./(2.*wY.^2) ; 
    b = (-sin(2.*theta))./(4.*wX.^2) + (sin(2.*theta))./(4.*wY.^2) ;
    c = (sin(theta).^2)./(2.*wX.^2) + (cos(theta).^2)./(2.*wY.^2) ; 
   
    ThisPeak(PtsInRange) = A.*exp(-(a.*(xPoints-xp).^2 + 2*b.*(xPoints - xp).*(yPoints-yp) + c.*(yPoints-yp).^2)) ;
    

    Z = Z+ThisPeak ; 

end



