%Sample script to detect and classify objects in PL images according to
%the groups defined in the paper "Efficient Analysis of Photoluminescence
%Images for the Identification of Single-Photon Emitters" 
%Total run time: 36 seconds

%% Begin by loading the PL image
dataStruct = load('1w_6um_circPol_highRes_072621.mat'); %load hBN PL image
pl = dataStruct.data.plScan(:,:,1);  
clockRate = dataStruct.data.clockRate;
pl = pl.*1000./clockRate;
pl = flipud(pl);
dataStruct.data.plScan = pl;
load('BlueSunset.mat') %load colormap 
%% Show the PL image
figure
imagesc(dataStruct.data.xCoords, dataStruct.data.yCoords,pl,[0 200])
xlabel('X Coordinates (\mum)')
ylabel('Y Coordinates (\mum)')
hold on 
c=colorbar;
c.Label.String = 'Photon Counts'
colormap(mymap)
res = mean(diff(dataStruct.data.xCoords));
title(['PL Image. HBN Sample. Resolution: ', num2str(100*res),'nm' ],'interpreter','none')
scalebar_v2(2.5,res,dataStruct.data.xCoords, dataStruct.data.yCoords,10,20,18,14) %scale bar


%% Adaptive object detection
thres = adaptthresh(mat2gray(pl), 0.3); %MATLAB built-in
binaryPl = imbinarize(mat2gray(pl),thres); %MATLAB built-in.

%%Alternatively, use universal threshold if adaptive threshold is unnecessarily slow/marks too many background points
%binaryPl = imbinarize(pl,20);

emitters = region_centers(dataStruct,binaryPl); %Employs MATLAB built-in function regionprops 

%% Apply Symmetric and Elliptical Gaussian Fits 
[BestFit,FitOptions,data] = PL_GaussFit(emitters,dataStruct,1:1:length(emitters),[0.1 0.35],0.15,5) ; 
[BestFitE,FitOptionsE,dataE] = PL_GaussFit_Ellipse(emitters,dataStruct,1:1:length(emitters),[0 1],BestFit) ;

%Comment lines 66-72 of PlGaussFit to disable regional maximum object
%detection. This step adds additional time and background sensitivity and
%is not necessary for well-isolated emitters (unlike hBN emitters). 

%% Plot object detection
centers = cat(1,BestFit.Posn)
p1 = plot(centers(:,1),centers(:,2),'*r')
title('Object Detection')

%% New figure for emitter marking 
figure
imagesc(dataStruct.data.xCoords, dataStruct.data.yCoords,pl,[0 200])
xlabel('X Coordinates (\mum)')
ylabel('Y Coordinates (\mum)')
hold on 
c=colorbar;
c.Label.String = 'Photon Counts'
colormap(mymap)
res = mean(diff(dataStruct.data.xCoords));
title(['PL Image. HBN Sample. Resolution: ', num2str(100*res),'nm' ],'interpreter','none')
scalebar_v2(2.5,res,dataStruct.data.xCoords, dataStruct.data.yCoords,10,20,18,14) %scale bar

%% Separate multiple emitters within one region of interest

%If multiple objects detected per cropped region of interest, it is necessary to unfurl eccentricity vector
%and align with Rchi2 from elliptical and symmetric fits
[l,w] = size(cat(1,BestFit.Width));
newChi = zeros(1,l);
newChiE = zeros(1,l);

if isfield(BestFitE,'eccentricity')
eccentricity = cat(2,BestFitE.eccentricity);
else
    eccentricity = []; 
end 

count = 1;
for i = 1:length(cat(1,BestFit.Rchi2))
    Rchi2 = BestFit(i).Rchi2; 
    Rchi2E = BestFitE(i).Rchi2; 
    e = BestFitE(i).eccentricity;
    np = length(e);
    for j = 1:np
    newChi(count) = Rchi2;
    count = count+1;
    end 
end

%% Emitter Classification
%set([p1,p1.Children],'Visible','off')

%Use results of fit to determine emitter groups
groups = which_group(newChi,cat(1,BestFit.snr),cat(1,BestFit.Width),[0.1 0.35], eccentricity);

%Label emitters with corresponding group 
centers = cat(1,BestFit.Posn);
centersE = cat(1,BestFitE.Posn); 

center1 =centers(groups==1,:);
for i  = 1:length(center1(:,1))
   text(center1(i,1)+0.2,center1(i,2),'A','Color',[1 0.64 0],'Fontsize',10,'FontWeight','bold');
end

center2 = centers(groups==2,:) ;
for i  = 1:length(center2(:,1))
  text(center2(i,1)+0.2,center2(i,2),'B','Color',[1 0.64 0],'FontSize',10,'FontWeight','bold');

end

center3 = centers(groups==3,:);
for i  = 1:length(center3(:,1))
   text(center3(i,1)+0.2,center3(i,2),'C','Color',[1 0.64 0],'FontSize',10,'FontWeight','bold');
end


center4 = centersE(groups==4,:);
for i  = 1:length(center4(:,1))
   text(center4(i,1)+0.2,center4(i,2),'D','Color',[1 0.64 0],'FontSize',10,'FontWeight','bold');
end

%% Mark only emitters in groups A-D

p1=plot(centers(groups~=0,1),centers(groups~=0,2),'x','MarkerSize',15,'Color',[1 0.64 0]);
hold on
p2=plot(centersE(groups==4,1),centersE(groups==4,2),'x','MarkerSize',15,'Color',[1 0.64 0]);

title('Emitter Classification - Annealed HBN Flake')



