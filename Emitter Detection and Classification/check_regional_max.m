
%Function uses MATLAB built-in imregionalmax to search for additional
%emitters missed with just adaptive thresholding
%Improves object detection for samples like HBN with uneven background and
%emitters of varying intensity. Cross-checks regional maximum with emitters
%detected by region_centers function
%Input: 
 %centersE: n*2 vector of initial emitter centers
 %EstHts: maximum intensity of each emitter
 %dataStruct: struct containing data.plScan field 
 %distance: threshold for distance of a new spot to an existing spot
 %intensity: minimum intensity for regional maximum considered valid spots
 
%Output: 
 %centersE: expanded vector of emitter  centers
 %EstHits: expanded vector of emitter intensity  

function [centersE,EstHts] = check_regional_max(centersE, EstHts, dataStruct, distance, intensity)

regionMax = imregionalmax(dataStruct.data.plScan);

pl = dataStruct.data.plScan;
rMiss = regionprops(regionMax,pl,'Centroid','MaxIntensity') ; 
centersR = cat(1,rMiss.Centroid) ; 

%convert to xy coords 
Xvec = dataStruct.data.xCoords;
Yvec = dataStruct.data.yCoords;

graphMinX = Xvec(1);
graphMaxX = Xvec(end);
pixelWidthX = length(Xvec);

graphMinY = Yvec(1);
graphMaxY = Yvec(end);
pixelWidthY = length(Yvec);


centersR(:,1) = centersR(:,1).*((graphMaxX-graphMinX)/pixelWidthX) + graphMinX ;
centersR(:,2) = centersR(:,2).*((graphMaxY-graphMinY)/pixelWidthY) + graphMinY ;

[B,I] = sort(cat(1,rMiss.MaxIntensity),'descend') ; 

centersR = centersR(I,:) ; 
intensities = B; 

%Check for duplicate centers
for iR = 1:length(centersR(:,1))
    %pull one center from rMiss
    checkR = centersR(iR,:);

        checkDistance = sqrt((checkR(1) - centersE(:,1)).^2 + (checkR(2) - centersE(:,2)).^2 )  ;
        
        %check against every center in emitters if within range of
        %distance limit to another center
        distanceLogic = sum(checkDistance>distance,'all') ; 
        
        %store intensity and position of new spot if it meets thresholds
        if distanceLogic == length(checkDistance) && intensities(iR) > intensity
           centersE = [centersE; checkR]  ;
           EstHts = [EstHts ; intensities(iR)] ;
        end    
end


end 
