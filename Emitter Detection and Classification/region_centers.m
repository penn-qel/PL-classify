%% Using built-in matlab regionprops, this function marks isolated emitters in a pl scan. 

%Input: struct "dataStruct" w/ "data" field containing plScan, and optionally
%x-coordinates (xCoords) and y-coordinates (yCoords) of the PL image
    %Example: dataStruct.data.yCoords = 1:1:10 

%Output: "emitters" struct  - each entry corresponds to an interconnected
%region detected by regionprops  

%Yields centroid position (in pixel and XY coords), area, bounding
%box (in pixel and XY coords), and other optional properties enabled by "regionprops" MATLAB function. 

 
function emitters = region_centers(dataStruct,binaryPl)
% define emitter struct and pull PL image from dataStruct
emitters = struct('Centroid',[],'Area',[],'CentroidXY',[],'BoundingBox',[],'regionalMax',[],'Orientation',[],'MeanIntensity',[],'Image',[],'Eccentricity',[],'BoundingBoxXY',[],'Extent',[],'MaxIntensity',[],'WeightedCentroidXY',[]);
pl = dataStruct.data.plScan ; 
pl = pl(:,:,1);

% Apply regionprops - built-in MATLAB function
r = regionprops(binaryPl,pl,'Area','Centroid','WeightedCentroid','BoundingBox','Eccentricity','Orientation','MajorAxisLength','MinorAxisLength','Perimeter','Extrema','Solidity','BoundingBox','Image','MaxIntensity','MeanIntensity');
data = dataStruct.data;


%% Emitter detection using pixel coordinates
if ~isfield(data, 'xCoords') || ~isfield(data, 'yCoords')

%Define rough lower/upper bounds for emitter size to filter small objects
[diffractionLim] = difLim(0.95,600) ;
    lowerLimit = (diffractionLim)^2 ;  
[diffractionLim] = difLim(0.95,850) ;
    upperLimit = (diffractionLim)^2 ;    


%initiate counter "k"     
k=1;

%Filter objects according to lower limit and expand bounding box if
%necessary
for i  = 1:length(r)
    if r(i).Area > lowerLimit && r(i).Area < upperLimit
       emitters(k).Centroid = r(i).Centroid;
       emitters(k).Area = r(i).Area;
       
   %If bounding box is too small, it will not fully encompass the emitter
   %and may skew the Gaussian fits. Expand the box if necessary
        box = r(i).BoundingBox;
        AiryLim  = 12;
        
        if box(3) < AiryLim || box(4) < AiryLim
         box =  r(i).BoundingBox;
         center = r(i).Centroid;
         
         if center(1)-5 > 0 
             box(1) = center(1) - 8;
         else
             box(1) = 0;
         end 
         
         if center(2)-5 > 0 
             box(2) = center(2) - 8;
         else
             box(2) = 0;
         end 
        
         box(3) = 20;
         box(4) = 20;
         r(i).BoundingBox = box;
        end 
        
       
       %Store properties
       emitters(k).MinorAxisLength = r(i).MinorAxisLength;
       emitters(k).MajorAxisLength = r(i).MajorAxisLength;
       emitters(k).Perimeter = r(i).Perimeter;
       emitters(k).Extrema = r(i).Extrema;
       emitters(k).BoundingBox = r(i).BoundingBox;
       
       emitters(k).Solidity = r(i).Solidity; 
       emitters(k).Image = r(i).Image; 
       emitters(k).MaxIntensity = r(i).MaxIntensity;
       emitters(k).WeightedCentroid = r(i).WeightedCentroid; 
       emitters(k).Eccentricity = r(i).Eccentricity;
       emitters(k).Orientation = r(i).Orientation;
       emitters(k).MeanIntensity = r(i).MeanIntensity;
       k = k+1;
    end 
end 
end 

%% Emitter detection using X and Y coordinates
if isfield(data,'xCoords') && isfield(data,'yCoords')
    
Xvec = dataStruct.data.xCoords;
Yvec = dataStruct.data.yCoords;
[nRows, nCols] = size(pl);


%Conversion from pixel coords to XY
graphMinX = Xvec(1);
graphMaxX = Xvec(end);
pixelWidthX = length(Xvec);

graphMinY = Yvec(1);
graphMaxY = Yvec(end);
pixelWidthY = length(Yvec);
res = mean(diff(Xvec)) ; 


%Define rough lower/upper bounds for emitter size to filter small objects
[diffractionLim,airyR] = difLim(0.95,600) ;
    lowerLimit = (diffractionLim/res)^2 ;
[diffractionLim,airyR] = difLim(0.95,850) ;
    upperLimit = (diffractionLim/res)^2;

%initiate counter "k"
k=1;

%Filter objects according to lower limit and expand bounding box if
%necessary
for i  = 1:length(r)
    if r(i).Area > lowerLimit  
       emitters(k).Centroid = r(i).Centroid;
       emitters(k).Area = r(i).Area;
       
   %If bounding box is too small, it will not fully encompass the emitter
   %and may skew the Gaussian fits. Expand the box if necessary
        box = r(i).BoundingBox;
        AiryLim  = 12;
        
         
        if box(3) < AiryLim || box(4) < AiryLim
         box =  r(i).BoundingBox;
         center = r(i).Centroid;
         
         if center(1)-5 > 0 
             box(1) = center(1) - 8;
         else
             box(1) = 0;
         end 
         
         if center(2)-5 > 0 
             box(2) = center(2) - 8;
         else
             box(2) = 0;
         end 
        
         box(3) = 14;
         box(4) = 14;
         r(i).BoundingBox = box;
        end 
       
       %Store properties in pixel coords
       emitters(k).MinorAxisLength = r(i).MinorAxisLength;
       emitters(k).MajorAxisLength = r(i).MajorAxisLength;
       emitters(k).Perimeter = r(i).Perimeter;
       emitters(k).Extrema = r(i).Extrema;
       emitters(k).BoundingBox = r(i).BoundingBox;
       
       emitters(k).Solidity = r(i).Solidity; 
       emitters(k).Image = r(i).Image; 
       emitters(k).MaxIntensity = r(i).MaxIntensity;
       emitters(k).WeightedCentroid = r(i).WeightedCentroid; 
       emitters(k).Eccentricity = r(i).Eccentricity;
       emitters(k).Orientation = r(i).Orientation;
       emitters(k).MeanIntensity = r(i).MeanIntensity;
       
       %Convert centroids, boundingboxes to XY coordinates
       centroidsPixels = r(i).Centroid;
       bBPixels = r(i).BoundingBox;
      
       centroidsXY(1) = centroidsPixels(1).*((graphMaxX-graphMinX)/pixelWidthX) + graphMinX ;
       centroidsXY(2) = centroidsPixels(2).*((graphMaxY-graphMinY)/pixelWidthY) + graphMinY ;
       
       centroidsPixels = r(i).WeightedCentroid;
       weightedcentroidsXY(1) = centroidsPixels(1).*((graphMaxX-graphMinX)/pixelWidthX) + graphMinX ;
       weightedcentroidsXY(2) = centroidsPixels(2).*((graphMaxY-graphMinY)/pixelWidthY) + graphMinY ;
       
       bBPixels(1) = bBPixels(1)*((graphMaxX-graphMinX)/pixelWidthX) + graphMinX ;
       bBPixels(2) = bBPixels(2)*((graphMaxY-graphMinY)/pixelWidthY) + graphMinY ;
       bBPixels(3:4) = bBPixels(3:4).*res;
       
       %Store XY coords 
       emitters(k).BoundingBoxXY = bBPixels;
       emitters(k).CentroidXY = centroidsXY; 
       emitters(k).WeightedCentroidXY = weightedcentroidsXY;
       k=k+1;
      
    end
end 
end 

%Keep PL image
emitters(1).pl = pl ; 
end





