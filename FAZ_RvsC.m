%% Circularity vs. Roundness validation for FAZ
% Created by: Jenna Grieshop
% Date: 6/17/2020 - 9/3/2020
% MCW AOIP Dr. Carroll
% This program creates simulated circles, ellipses, and polygons that can
% be measured the same way as an actual FAZ would be. It also reads in
% actual FAZ segmentations to measure them. Various metrics are produced
% including area, perimeter, roundness, circularity, etc.
%
% Notes:
%
% **Functions**
%
% Ellipse Family Function
%  - elliFam(maj, min, AR1, AR2, folder)
%  - maj: major axis
%  - min: minor axis
%  - AR1: first number of aspect ratio
%  - AR2: second number of aspect ratio
%  - folder: folder to store results in
% 
% Theoretical and Mathematical Calculations Function
%  - calc(sides,shape,radius, major, minor, coords)
%
%  *NOTE: if a parameter does not apply for shape, plug in 1)*
%
%  - sides: number of sides the shape has
%  - shape: number referencing the shape that is being analyzed
%  (1 - circle, 2 - ellipse, 3 - hexagon, 4 - octagon, 5 - actual FAZ)
%  - radius: radius of circle; use as sidelength for polygons
%  - major: major axis of ellipse
%  - minor: minor axis of ellipse
%  - coords: coordinates for the shape
%
% RegionProps and inscribed/circumscribed circles Function
%  - RP(coords, filename)
%  - coords: coordinates for the shape
%  - filename: filename for the mask to be saved as
%
% Connect to ImageJ Function
%  - imageJ(shape,rNum,major,minor,sL,name)
%
%  *NOTE: if a parameter does not apply for shape, plug in 1)*
%  
%  - shape: number referencing the shape that is being analyzed
%  (1 - circle, 2 - ellipse, 3 - hexagon, 4 - octagon, 5 - actual FAZ)
%  - rNum: radius
%  - major: major axis
%  - minor: minor axis
%  - sL: side length for polygons
%  - name: name of actual FAZ


clear all
clear global
clc
close all

ccPath = uigetdir('.','Select Circumscribed Circle directory');
addpath (ccPath);
folderPath = uigetdir('.','Select Results Location');
%need to addpath for MATLAB to find ImageJ in the scripts folder for the
%fiji app. If it doesn't have this it can't run ImageJ
fsPath = uigetdir('.','Select scripts folder in fiji-win64 -> Fiji.app');
addpath (fsPath);
folder = (folderPath);

%% Actual FAZs

% Batch Upload feature:
% User selects desired folder:
root_path = uigetdir('.','Select directory FAZ data');
root_dir = dir(root_path);
root_dir = struct2cell(root_dir)';
addpath(genpath(root_path));

% looks for the Batch Info spreadsheet that contains the ID initial, ID #, eye, axial length, and device #:
batch_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), 'Batch_Info')),:);
% Loads in the Batch Info file:
   batch_info = table2cell(readtable(fullfile(batch_dir{1,2},batch_dir{1,1})));
   
% looks for all the tifs
tif_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), '.tif')),:);
% looks for all the corresponding csvs
csv_dir = root_dir(...
    ~cellfun(@isempty, strfind(root_dir(:,1), '.csv')),:);

labelF = strings(length(csv_dir)-1,1); %Preallocate label array
labelF(1,1) = 'Image'; %Header for image labels for excel sheet
labelCount = 2; %set count to 2 so that Image does not get overwritten in loop
for i=1:size(batch_info, 1)
    
% Pulling information from batch info:
    Subject_ID = batch_info{i,1};
    eye = batch_info{i,2};
    axial = batch_info{i,3};
    device = batch_info{i,4};
       
% Pull corresponding scan and coord files based on subject ID and eye
    subj_tifs = tif_dir(...
    ~cellfun(@isempty, strfind(tif_dir(:,1), Subject_ID)),:);
    
    octa_img = char(subj_tifs(...
    ~cellfun(@isempty, strfind(subj_tifs(:,1), eye)), 1));
    
    subj_csv = csv_dir(...
    ~cellfun(@isempty, strfind(csv_dir(:,1), Subject_ID)),:);

    coordname = char(subj_csv(...
    ~cellfun(@isempty, strfind(subj_csv(:,1), eye)), 1));

%Input x and y dimensions of the image to provide a reference of the FAZ. Dimensions must be greater than the coordinates of the FAZ.
    if device == 1
        imageX = 304;
        imageY = 304;
        scanwidth = 3000; %scan width in microns
        assumedAL = 23.95; %assumed axial length of optovue
        scale = ((axial*scanwidth)/assumedAL)/imageX;
    elseif device == 2
        imageX = 1024;
        imageY = 1024;
        scanwidth = 3000; %scan width in microns
        assumedAL = 24.46; %assumed axial length of optovue
        scale = ((axial*scanwidth)/assumedAL)/imageX;
    else
        imageX = str2num(input('Please input image width in pixels: ', 's'));
        imageY = str2num(input('Please input image height in pixels: ','s'));
        scale = str2num(input('Please input image scale in um/pp: ','s'));
    end

for ii=1:size(octa_img,1)
        
    current_img = deblank(octa_img(ii,:));
    current_coords = deblank(coordname(ii,:));
    
%Read the coordinates from the .csv file.
coords = importdata(current_coords);
octa = imread(current_img);

fCalc = calc(scale,1,5,1,1,1,coords); %returned output from calculations function

filenameF = strcat('matlab_', current_img);
excelnameF = strcat('ActualFAZ_metrics.xls');

fRes = RP(scale,5, 100, coords, filenameF, folder); %returned output from RegionProps funtion

labelF(labelCount,1) = current_img; %saving image names for excel labels

IJFAZRes = imageJ(folder,5, 1, 1, 1, 1, filenameF); %returned output from imageJ function

labelCount=labelCount+1; 
end
ii = 1;
end

MIJ_R1 = {'Imj_Area', 'Imj_Perimeter', 'Imj_Major Axis', 'Imj_Minor Axis', 'Imj_Angle', 'Imj_Circularity', 'Imj_AR', 'Imj_Roundness', 'Imj_Solidity'};
ImageJResFAZ = [MIJ_R1;IJFAZRes]; %ccombine imageJ labels with data

fResNum = [fCalc, fRes]; %just numbers
header = {'Poly_Area','Coord_Perimeter','Calculated_Circularity','RP_Area','RP_Perimeter', 'RP_OldPerimeter', 'RP_Circularity', 'RP_Major Axis', 'RP_Minor Axis', 'RP_Roundness', 'CircleRatio_Roundness'};
FAZResults = [header; num2cell(fResNum)]; %combine labels with data
FAZResults = [FAZResults, ImageJResFAZ]; %combine both sets of data
FAZResults = [labelF, FAZResults]; %concatenate radius labels with data

xlswrite(fullfile(folder,excelnameF), FAZResults)

clear global;
close all;

%% Circle Family
%radius
scale = 1; %arbritrary
r = [2 3 4 5 6 7 8 9 10 11 14 17 20 23 26 29 32 35 38 41 44 47 50 75 100 125 150 175 200 300 400 500 600 700 800 900 1000]; %radius array to test circles
labelC = strings(length(r),1); %Preallocate label array
labelC(1,1) = 'Radius'; %header for radius labels for excel organization

for i=1:length(r) %for loop to go through each of the circles   
    
t = linspace(0, 2*pi, 300); % theta; vector of 300 evenly spaced points bw 0 and 2pi
R = string(r(i)); %converts r to a string for naming purposes
xy = [1650, 1650]; % center of mask
xx = r(i)*cos(t)+xy(1); 
yy = r(i)*sin(t)+xy(2); 

percirc = [xx;yy]; %2x300 vector containing the xx and yy values
coords = transpose(percirc); %rows->columns, columns->rows; changed to 300x2 vector 

cCalc = calc(scale,1,1,r(i),1,1,coords); %returned output from calculations function

filenameC = strcat('Circle_', R,'_FAZ_mask.tif');
excelnameC = strcat('CircleFamily_300ptsampling_FAZmetrics.xls');

cRes = RP(scale,1, r(i), coords, filenameC, folder); %returned output from RegionProps funtion

labelC(i+1,1) = R; %concatenate radii for excel label

IJCircRes = imageJ(folder,1, r(i), 1, 1, 1,1); % returned output from imageJ function
end

MIJ_R1 = {'Imj_Area', 'Imj_Perimeter', 'Imj_Major Axis', 'Imj_Minor Axis', 'Imj_Angle', 'Imj_Circularity', 'Imj_AR', 'Imj_Roundness', 'Imj_Solidity'};
ImageJResCirc = [MIJ_R1;IJCircRes]; %concatenate imageJ data with labels

cResNum = [cCalc, cRes]; %just numbers
header = {'Theo_Area','Theo_Perimeter','Theo_Circularity', 'Theo_Roundness', 'Poly_Area','Coord_Perimeter','RP_Area','RP_Perimeter', 'RP_OldPerimeter', 'RP_Circularity', 'RP_Major Axis', 'RP_Minor Axis', 'RP_Roundness', 'CircleRatio_Roundness'};
circResults = [header; num2cell(cResNum)]; %concatenate labels with data
circResults = [circResults, ImageJResCirc]; %concatenate both sets of data
circResults = [labelC, circResults]; %concatenate radius labels with data

xlswrite(fullfile(folder,excelnameC), circResults)

clear global;
close all;

%% Ellipse Families
AR10 = '10';
AR6 = '6';
AR4 = '4';
AR8 = '8';

%ellipse family for Aspect Ratio 10:6
maj10_6 = [2.581988897 3.872983346 5.163977795 6.454972244 7.745966692 9.036961141 10.32795559 11.61895004 12.90994449 14.20093894 18.07392228 21.94690563 25.81988897 29.69287232 33.56585567 37.43883901 41.31182236 45.18480571 49.05778905 52.9307724 56.80375574 60.67673909 64.54972244 96.82458366 129.0994449 161.3743061 193.6491673 225.9240285 258.1988897 387.2983346 516.3977795 645.4972244 774.5966692 903.6961141 1032.795559 1161.895004 1290.994449]; %major radius always larger than min
min10_6 = [1.549193338 2.323790008 3.098386677 3.872983346 4.647580015 5.422176685 6.196773354 6.971370023 7.745966692 8.520563362 10.84435337 13.16814338 15.49193338 17.81572339 20.1395134 22.46330341 24.78709342 27.11088342 29.43467343 31.75846344 34.08225345 36.40604345 38.72983346 58.09475019 77.45966692 96.82458366 116.1895004 135.5544171 154.9193338 232.3790008 309.8386677 387.2983346 464.7580015 542.2176685 619.6773354 697.1370023 774.5966692]; %minor axis radius

%ellipse family for Aspect Ratio 10:4
maj10_4 = [3.16227766 4.74341649 6.32455532 7.90569415 9.486832981 11.06797181 12.64911064 14.23024947 15.8113883 17.39252713 22.13594362 26.87936011 31.6227766 36.36619309 41.10960958 45.85302607 50.59644256 55.33985905 60.08327554 64.82669203 69.57010852 74.31352501 79.0569415 118.5854123 158.113883 197.6423538 237.1708245 276.6992953 316.227766 474.341649 632.455532 790.569415 948.6832981 1106.797181 1264.911064 1423.024947 1581.13883];
min10_4 = [1.264911064 1.897366596 2.529822128 3.16227766 3.794733192 4.427188724 5.059644256 5.692099788 6.32455532 6.957010852 8.854377448 10.75174404 12.64911064 14.54647724 16.44384383 18.34121043 20.23857703 22.13594362 24.03331022 25.93067681 27.82804341 29.72541001 31.6227766 47.4341649 63.2455532 79.0569415 94.86832981 110.6797181 126.4911064 189.7366596 252.9822128 316.227766 379.4733192 442.7188724 505.9644256 569.2099788 632.455532];

%ellipse family for Aspect Ratio 10:8
maj10_8 = [2.236067977 3.354101966 4.472135955 5.590169944 6.708203932 7.826237921 8.94427191 10.0623059 11.18033989 12.29837388 15.65247584 19.00657781 22.36067977 25.71478174 29.06888371 32.42298567 35.77708764 39.13118961 42.48529157 45.83939354 49.1934955 52.54759747 55.90169944 83.85254916 111.8033989 139.7542486 167.7050983 195.655948 223.6067977 335.4101966 447.2135955 559.0169944 670.8203932 782.6237921 894.427191 1006.23059 1118.033989];
min10_8 = [1.788854382 2.683281573 3.577708764 4.472135955 5.366563146 6.260990337 7.155417528 8.049844719 8.94427191 9.838699101 12.52198067 15.20526225 17.88854382 20.57182539 23.25510697 25.93838854 28.62167011 31.30495168 33.98823326 36.67151483 39.3547964 42.03807798 44.72135955 67.08203932 89.4427191 111.8033989 134.1640786 156.5247584 178.8854382 268.3281573 357.7708764 447.2135955 536.6563146 626.0990337 715.5417528 804.9844719 894.427191];

%call ellipse family function to create and get data for the ellipses
elliFam(maj10_6, min10_6, AR10, AR6, folder);
clear global;
close all;
elliFam(maj10_4, min10_4, AR10, AR4, folder);
clear global;
close all;
elliFam(maj10_8, min10_8, AR10, AR8, folder);
clear global;
close all;

%% Hexagon
scale = 1; %arbritrary
sideLength = 109.9636111; %side length of regular hexagon with area of the r=100 circle
sL = string(sideLength); %for naming purposes
labelP = {'Polygon'; strcat('Hex_',sL)};
hex = nsidedpoly(6,'SideLength', sideLength); %creates regular hexagon

[x,y] = boundary(hex); %get the coordinates from the nsidedpoly
xy = [1650, 1650]; %center of mask
x = x + xy(1);
y = y + xy(2);

coords = [x,y]; 

hCalc = calc(scale,6,3,sideLength,1,1,coords); %returned output from calculations function, pass sidelength as the radius

filenameH = strcat('Hexagon','_', sL ,'_FAZ_mask.tif');
excelnameP = strcat('Polygon_FAZmetrics.xls');

hRes = RP(scale,1,sideLength, coords, filenameH, folder); %returned output from RP funtion

IJHexRes = imageJ(folder,3,1,1,1,sL,1); %retuned output from ImageJ function

MIJ_R1 = {'Imj_Area', 'Imj_Perimeter', 'Imj_Major Axis', 'Imj_Minor Axis', 'Imj_Angle', 'Imj_Circularity', 'Imj_AR', 'Imj_Roundness', 'Imj_Solidity'};
ImageJResHex = [MIJ_R1;IJHexRes]; %concatenate imageJ labels with data

hResNum = [hCalc, hRes]; %just numbers
header = {'Theo_Area','Theo_Perimeter','Theo_Circularity', 'Poly_Area','Coord_Perimeter','RP_Area','RP_Perimeter', 'RP_OldPerimeter', 'RP_Circularity', 'RP_Major Axis', 'RP_Minor Axis', 'RP_Roundness', 'CircleRatio_Roundness'};
hexResults = [header; hResNum]; %concatenate other labels with other data
hexResults = [hexResults, ImageJResHex]; %concatenate both sets of data
hexResults = [labelP, hexResults]; %concatenate sidelength label with data

xlswrite(fullfile(folder,excelnameP), hexResults);

close all;

%% Octagon 
scale = 1; %arbritrary
sideLength = 80.66257759; %side length of regular octagon with area of the r=100 circle
sL = string(sideLength); %for naming purposes
labelP(3,1) = {strcat('Oct_',sL)};
oct = nsidedpoly(8,'SideLength', sideLength); %creates regular octagon

[x,y] = boundary(oct); %get coordinates from nsidedpoly
xy = [1650, 1650]; %center of mask
x = x + xy(1);
y = y + xy(2);

coords = [x,y]; 

oCalc = calc(scale,8,4,sideLength,1,1,coords); %returned output from calculations function, pass sideLength as the radius

filenameO = strcat('Octagon', '_', sL ,'_FAZ_mask.tif');
excelnameP = strcat('Polygon_FAZmetrics.xls');

oRes = RP(scale,1,sideLength,coords, filenameO, folder); %returned output from RP funtion

IJOctRes = imageJ(folder,4,1,1,1,sL,1); %returned output from the ImageJ function

MIJ_R1 = {'Imj_Area', 'Imj_Perimeter', 'Imj_Major Axis', 'Imj_Minor Axis', 'Imj_Angle', 'Imj_Circularity', 'Imj_AR', 'Imj_Roundness', 'Imj_Solidity'};
ImageJResOct = [MIJ_R1;IJOctRes]; %concatenate imageJ data with labels

oResNum = [oCalc, oRes]; %just numbers
header = {'Theo_Area','Theo_Perimeter','Theo_Circularity', 'Poly_Area','Coord_Perimeter','RP_Area','RP_Perimeter', 'RP_OldPerimeter', 'RP_Circularity', 'RP_Major Axis', 'RP_Minor Axis', 'RP_Roundness', 'CircleRatio_Roundness'};
octResults = [header; oResNum]; %concatenate other labels with other data
octResults = [octResults, ImageJResOct]; % concatenate both sets of data
octResults = [labelP, octResults]; %concatenate radius labels with data

xlswrite(fullfile(folder,excelnameP), octResults);

clear global;
close all;

%% Ellipse Family Function
function elliFam(maj, min, AR1, AR2, folder)
scale = 1; %arbritrary
labelE = strings(length(maj),1); %preallocate label array
labelE(1,1) = 'Radii'; %for excel organization

for i=1:length(maj)%for loop to go through each of the ellipses 

t = linspace(0, 2*pi, 300); % theta; vector of 300 evenly spaced points bw 0 and 2pi
major = string(round(maj(i),0)); %converts maj to string for naming purposes
minor = string(round(min(i),0)); 
majN = round(maj(i),0);
minN = round(min(i),0);
radii = strcat(major, '_', minor); %combining major and minor radius lengths for excel label
xy = [1650, 1650]; %center of mask
xx = maj(i)*cos(t)+xy(1); 
yy = min(i)*sin(t)+xy(2); 

perelli = [xx;yy]; %2x300 vector containing the xx and yy values
coords = transpose(perelli); %rows->columns, columns->rows; changed to 300x2 vector

eCalc = calc(scale,1,2,1,maj(i),min(i), coords); %returned output from calculations function   

filenameE = strcat('Ellipse_', 'Maj_', major, 'Min_', minor,'_FAZ_mask.tif');
excelnameE = strcat('EllipseFamily_', AR1,'_', AR2,'_300ptsampling_FAZmetrics.xls');

eRes = RP(scale,1,maj,coords, filenameE,folder); %returned output from RegionProps funtion

labelE(i+1,1) = radii; %concatenate radii for excel label

IJElliRes = imageJ(folder,2,1,majN, minN,1,1); %returned output from imageJ function
end

MIJ_R1 = {'Imj_Area', 'Imj_Perimeter', 'Imj_Major Axis', 'Imj_Minor Axis', 'Imj_Angle', 'Imj_Circularity', 'Imj_AR', 'Imj_Roundness', 'Imj_Solidity'};
ImageJResElli = [MIJ_R1;IJElliRes]; %concatenate imageJ data and labels

eResNum = [eCalc, eRes]; %just numbers
header = {'Theo_Area', 'Theo_Perimeter', 'Theo_Circularity', 'Theo_Roundness', 'Poly_Area','Coord_Perimeter','RP_Area','RP_Perimeter', 'RP_OldPerimeter', 'RP_Circularity', 'RP_Major Axis', 'RP_Minor Axis', 'RP_Roundness', 'CircleRatio_Roundness'};
elliResults = [header; num2cell(eResNum)]; %concatenate other data and labels
elliResults = [elliResults, ImageJResElli]; %concatenate both sets of data
elliResults = [labelE, elliResults]; %concatenate radius labels with data

xlswrite(fullfile(folder,excelnameE), elliResults)

close all;
end
%% Theoretical and Mathematical Calculations Function
function output3 = calc(scale,sides,shape,radius, major, minor, coords)%radius is used as sidelength for polygons
global output3;

xcoords = coords(:,2); 
ycoords = coords(:,1);

if(shape ==1) %circle
    %Calculating Theoretical Area
    area = pi*(radius^2);
    
    %Calculating Theoretical Perimeter
    perimeter = 2*pi*radius;
    
    %Calculating Theoretical Roundness
    round = radius/radius;
else 
    if(shape == 2) %ellipse
        %Calculating Theoretical area
        area = pi*major*minor;
        
        %Calculating Theoretical Perimeter
        top = ((major - minor)^2);
        bottom = ((major + minor)^2);
        h = top/bottom;
        perimeter = (pi *(major + minor)) * ((1 + ((3*h)/(10 + sqrt(4-(3*h))))));
        
        %Calculating Theoretical Roundness
        round = minor/major;
else
    if(shape == 3) %hexagon
        %Calculating Theoretical Area
        area = ((3*sqrt(3))/2)*(radius^2);
        
        %Calculating Theoretical Perimeter
        perimeter = sides*radius; 
        
        %Calculating Roundness **Cannot calculate: sidelength != radius
        %round = radius/radius;
else
    if(shape == 4) %octagon
        %Calculating Theoretical Area
        area = 2*(1+sqrt(2))*(radius^2);
        
        %Calculating Theoretical Perimeter
        perimeter = sides*radius; 
        
        %Calculating Roundness **Cannot calculate: sidelength != radius
        %round = radius/radius; 
        
    end
    end
    end
end
    
%Calculating Area via polyarea
if shape == 5
    polyArea = (polyarea(xcoords,ycoords))*(scale*scale/1000000);
else
    polyArea = polyarea(xcoords,ycoords);
end

%Calculating Perimeter via coordinates
xcoords_2 = [xcoords;xcoords(1)];
ycoords_2 = [ycoords;ycoords(1)];
x_diff = diff(xcoords_2); %calculates differences between adjacent elements
y_diff = diff(ycoords_2);
segment_lengths = sqrt(x_diff.^2 + y_diff.^2); %distance formula between points
if shape == 5
    coordPerimeter = (sum(segment_lengths))*(scale/1000);
else
coordPerimeter = sum(segment_lengths); %sum of distances
end


% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % Calculating Acircularity
% radius = sqrt((area/pi));
% PCPerimeter = (2*pi*radius);
% Acirc = perimeter/PCPerimeter;
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%Calculating Theoretical/Calculated Circularity
if(shape ==5)
    circ = ((4*pi*polyArea)/(coordPerimeter^2));
else
    circ = ((4*pi*area)/(perimeter^2));
end

%depending on the shape, different data has been calculated
if(shape == 1 || shape == 2)
    output(:,1) = area;
    output(:,2) = perimeter;
    output(:,3) = circ;
    output(:,4) = round;
    output(:,5) = polyArea;
    output(:,6) = coordPerimeter;
else
    if(shape == 3 || shape == 4)
    output(:,1) = area;
    output(:,2) = perimeter;
    output(:,3) = circ;
    output(:,4) = polyArea;
    output(:,5) = coordPerimeter;  
else
    output(:,1) = polyArea;
    output(:,2) = coordPerimeter;
    output(:,3) = circ;
    end
end
output3 = [output3;num2cell(output)];

end
%% RegionProps and inscribed/circumscribed circles Function
function output2 = RP(scale,shape, r, coords, filename, folder)
global output2;

%used in mask for size of the circle in rows and columns
if shape == 5
    imageX = 304; 
    imageY = 304; %originally here for FAZ compatibility - changed to 3300
    %because of 1000 as max radius
else
    imageX = 3300;
    imageY = 3300;
end

%mask set up
mask = poly2mask(coords(:,1), coords(:,2), imageY, imageX);
    
xcoords = coords(:,2);
ycoords = coords(:,1);
xycoord = [ycoords,xcoords];

%outer and inner circle
[rad,~,~] = ExactMinBoundCircle(xycoord);
rad = rad * (scale/1000);%don't need scale - it originally converted to
%mm
outerCirc = (pi*(rad^2));
mask_2 = ~mask;
ContourImage = (255*mask_2);
[rad2,~,~] = max_inscribed_circle(ContourImage);
rad2 = rad2 * (scale/1000); %don't need scale - it originally converted to
%mm
innerCirc = (pi*(rad2^2));

circRatio = innerCirc/outerCirc;

% RegionProps
%Creating a mask to use for regionprops
imwrite(mask*255, fullfile(folder,filename)); %mask*255 for Black and White %writes mask to tif image

%Regionprops Validation
metrics = regionprops(mask,'Area','Perimeter','PerimeterOld', 'Circularity','MajorAxisLength','MinorAxisLength');
if shape == 5
    RP_Area = ((getfield(metrics, 'Area')) * (scale*scale/1000000));
    RP_Perimeter = (getfield(metrics, 'Perimeter'))*(scale/1000);
    RP_OldPerimeter = (getfield(metrics, 'PerimeterOld'))*(scale/1000);
    Major_Axis = (getfield(metrics, 'MajorAxisLength'))*(scale/1000);
    Minor_Axis = (getfield(metrics, 'MinorAxisLength'))*(scale/1000);
else
    RP_Area = ((getfield(metrics, 'Area'))); 
    RP_Perimeter = ((getfield(metrics, 'Perimeter'))); 
    RP_OldPerimeter = (getfield(metrics, 'PerimeterOld'));
    Major_Axis = (getfield(metrics, 'MajorAxisLength'));
    Minor_Axis = (getfield(metrics, 'MinorAxisLength'));
end
RP_Circularity = ((getfield(metrics, 'Circularity')));
Roundness = Minor_Axis/Major_Axis;

% Save Data
%Create matrix with metrics
output(:,1) = RP_Area;
output(:,2) = RP_Perimeter;
output(:,3) = RP_OldPerimeter;
output(:,4) = RP_Circularity;
output(:,5) = Major_Axis;
output(:,6) = Minor_Axis;
output(:,7) = Roundness;
output(:,8) = circRatio;

output2 = [output2;num2cell(output)];

end

%% Connect to ImageJ Function
% For actual FAZs need to convert data from pixel to mm after in data
% analysis. Area * (scale*scale/1000000), else * (scale*scale/1000)
function imjOutput = imageJ(folder,shape,rNum,major,minor,sL,name)
global imjOutput;

%paths needed for MIJI
javaaddpath 'C:\Program Files\MATLAB\R2019b\java\mij.jar' 
javaaddpath 'C:\Program Files\MATLAB\R2019b\java\ij-1.53c.jar'

MIJ.start; %starts ImageJ

filename1 = strcat(folder, '\');
%depending on the shape, the file names to use are different
if(shape == 1)%circle
        filename2 = sprintf('Circle_%d_FAZ_mask.tif', rNum);

else
    if(shape == 2)%ellipse
        filename2 = sprintf('Ellipse_Maj_%dMin_%d_FAZ_mask.tif', major, minor); 
else
    if(shape == 3)%hexagon
        filename2 = sprintf('Hexagon_%s_FAZ_mask.tif' ,sL);
else
    if(shape == 4)%octagon
        filename2 = sprintf('Octagon_%s_FAZ_mask.tif' ,sL);  
else
    if(shape == 5) %actual FAZ
        filename2 = sprintf(name);
    end
    end
    end
    end
end
%opens the image
MIJ.run('Open...', strcat('path=[',filename1,filename2,']'));

MIJ.run('Threshold...');
MIJ.run('Convert to Mask'); %Apply button
MIJ.run('Convert to Mask'); %Apply button again to complete it
MIJ.run('Close'); %close Threshold

MIJ.run('Set Measurements...', 'area perimeter shape fit'); %set measurements to be recorded
% feret''s
MIJ.run('Analyze Particles...', 'show=Outlines display'); %run analyze particles


MIJ_R2 = MIJ.getResultsTable; %saves results
MIJ_R2 = num2cell(MIJ_R2); %need to convert to cell to combine them
imjOutput= [imjOutput; MIJ_R2]; %combines results with labels

MIJ.run('Clear Results');
MIJ.closeAllWindows(); %closes all windows
MIJ.exit; %closes imagej
end