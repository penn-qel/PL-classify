%Function to compute chi squared for goodness-of-fit calculation 
%Input: data, BestFit structs from PlGaussFit or PlGaussFitEllipse
%Output: 
    %chiSquareSum: sum of chi^2 computed for each data point 
    %chiSquare: Array of chiSquare for each data point
    %nPoints: data points in the region of interest
    
function [chiSquareSum, chiSquare,nPoints] = chi_square(data, BestFit)

%Pull model function from the fit 
ModelFunction = BestFit.ModelFunction;

%Pull measured counts
X = data.X; %xCoords
Y = data.Y; %yCoords
Z = data.Z; %photoluminescence counts 
  

%Calculate fit values
[Xmat,Ymat] = meshgrid(X,Y);
functionValues = ModelFunction(Xmat,Ymat) ;
variance = Z; 

chiSquare = (Z-functionValues).^2./(variance);
chiSquareSum = sum(chiSquare(chiSquare~=Inf),'all');
nPoints = sum(chiSquare~=Inf,'all'); 

end 