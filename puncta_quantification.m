% Written By Sneha Shankar, edited by Dana Zemel

%% ==========================================
    % Change these to correct file name of your images
    greenFileName = 'FILENAME.tif';
    redFileName = 'FILENAME.tif';
   %==========================================
%% 
% read the files

%Load FingR labled or Channel 1 tif file: 
green = Tiff(greenFileName,'r+');  
greendata = read(green);
green = greendata;

%Load antibody stain labeled or Channel 2 tif file: 
red = Tiff(redFileName,'r+'); 
reddata = read(red);
red = reddata;

%% ==========================================
    %Set thershold intensity for each channel: - Change these based n your
    % images
    th_1 = 100 ; %Threshold for FingR (Ch1) image 
    th_2 =  100; %Threshold for Stain (Ch2) image 
 % ==========================================
%% locate and quantify puncta

% Blur images for cleaner results:
reddata =imgaussfilt(reddata,2);
greendata =imgaussfilt(greendata,2);

%Set values below puncta threshold to zero
greendata(greendata<=th_1)=0;
reddata(reddata<=th_2)=0;

% set above thershould pixels to 1 
greendata(greendata>th_1)=1;
reddata(reddata>th_2)=1;

SE = strel('disk',1);
greendata2= imopen(greendata,SE);
greendata2= imclose(greendata2,SE);
reddata2= imopen(reddata,SE);
reddata2= imclose(reddata2,SE);

%Find connected regions within threshold area to find individual puncta
green_cc=regionprops(logical(greendata2),'Area','PixelIdxList','Centroid');
green_cc = green_cc([green_cc(:).Area]<'SET PIXEL AREA LIMIT');

red_cc=regionprops(logical(reddata2),'Area','PixelIdxList','Centroid');
red_cc = red_cc([red_cc(:).Area]<'SET PIXEL AREA LIMIT');

%Find center coordinates of puncta for both channels
green_x=zeros(numel(green_cc),1);
green_y=zeros(numel(green_cc),1);
for i=1:numel(green_cc)
green_x(i,1)=green_cc(i).Centroid(1,1);
green_y(i,1)=green_cc(i).Centroid(1,2);
end


red_x=zeros(numel(red_cc),1);
red_y=zeros(numel(red_cc),1);
for i=1:numel(red_cc)
red_x(i,1)=red_cc(i).Centroid(1,1);
red_y(i,1)=red_cc(i).Centroid(1,2);
end

%Find number of colocalized puncta
colocalized=0; 
for i=1:length(green_x)
    for j=1:length(red_x)
     if sqrt((green_x(i,1)-red_x(j,1))^2+(green_y(i,1)-red_y(j,1))^2)<='SET PIXEL DISTANCE' 
        colocalized=colocalized+1;
     end 
    end 
end 

%Find distribution of puncta within each channel
green_distance=zeros(length(green_x),1);
for i=1:length(green_x)
green_distance(i,1) = sqrt((green_x(i,1)-green_x(1,1))^2+(green_y(i,1)-green_y(1,1))^2);
end

red_distance=zeros(length(red_x),1);
for i=1:length(red_x)
red_distance(i,1) = sqrt((red_x(i,1)-red_x(1,1))^2+(red_y(i,1)-red_y(1,1))^2);
end

fprintf('FingR (CH1) puncta: %d\n', numel(green_cc));
fprintf('Antibody stained (CH2) puncta: %d\n', numel(red_cc));
fprintf('overlap puncta: %d\n', colocalized);
fprintf('Precision/Selectivity: %f\n', 100*colocalized/(numel(green_cc)));
fprintf('Viral Efficiency: %f\n', 100*colocalized/(numel(red_cc)));
