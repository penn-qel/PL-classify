%Plots scale bar
%Input
 %l: length of the scale bar in microns
 %res: resolution of the image in microns
 %xCoords: vector of xCoordinates for the image
 %yCoords: vector of yCoordinates for the image
 
 %yBuffer: scalar buffer from the end of yCoords 
 %xBuffer: scalar buffer from the end of xCoords
 %textBuffer: scalar buffer from the end of xCoords for the text
 %separationBuffer: scalar buffer between the text and the scale bar

 
function scalebar_v2(l,res,xCoords,yCoords,yBuffer,xBuffer,textBuffer,separationBuffer)
scalebar = l ; 
x = floor(scalebar/res); 
plot([xCoords(end-xBuffer), xCoords(end-x-xBuffer)],[yCoords(end-yBuffer), yCoords(end-yBuffer)],'-w','LineWidth',3)
hold on 
text(xCoords(end-floor(x/2)-textBuffer-xBuffer), yCoords(end-separationBuffer-yBuffer),[num2str(scalebar),'\mum'],'Color','white','FontSize',14)
end 

 