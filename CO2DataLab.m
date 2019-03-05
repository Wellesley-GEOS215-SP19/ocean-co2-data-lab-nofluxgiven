%% THIS IS THE CORRECT ONE Grace Callahan & Shreya Parjan

% Instructions: Follow through this code step by step, while also referring
% to the overall instructions and questions from the lab assignment sheet.

%% 1. Read in the monthly gridded CO2 data from the .csv file
% The data file is included in your repository as �LDEO_GriddedCO2_month_flux_2006c.csv�
% Your task is to write code to read this in to MATLAB
% Hint: you can again use the function �readtable�, and use your first data lab code as an example.
filename = 'LDEO_GriddedCO2_month_flux_2006c.csv';
CO2data = readtable(filename);

%% 2a. Create new 3-dimensional arrays to hold reshaped data
%Find each unique longitude, latitude, and month value that will define
%your 3-dimensional grid
longrid = unique(CO2data.LON); %finds all unique longitude values
latgrid = unique(CO2data.LAT); %<-- following the same approach, find all unique latitude values
mongrid = unique(CO2data.MONTH); %<-- following the same approach, find all unique months

%Create empty 3-dimensional arrays of NaN values to hold your reshaped data
    %You can make these for any variables you want to extract - for this
    %lab you will need PCO2_SW (seawater pCO2) and SST (sea surface
    %temperature)
llmGrid = NaN*zeros(length(longrid),length(latgrid),length(mongrid));
pCO2Grid = NaN*zeros(length(longrid),length(latgrid),length(mongrid));
sstGrid = NaN*zeros(length(longrid),length(latgrid),length(mongrid));
%<--

%% 2b. Pull out the seawater pCO2 (PCO2_SW) and sea surface temperature (SST)
%data and reshape it into your new 3-dimensional arrays
for i = 1:length(CO2data.LON)
   m = CO2data.MONTH(i);
   j = find(longrid == CO2data.LON(i));
   k = find(latgrid == CO2data.LAT(i));
   l = find(mongrid == m);
   % Now we fill the grid
   pCO2Grid(j,k,l) =(CO2data.PCO2_SW(i));
   
   %SST 
   sstGrid(j,k,m) =(CO2data.SST(i));
   
end

%% 3a. Make a quick plot to check that your reshaped data looks reasonable
%Use the imagesc plotting function, which will show a different color for
%each grid cell in your map. Since you can't plot all months at once, you
%will have to pick one at a time to check - i.e. this example is just for
%January

%imagesc(sstGrid(:,:,1))

%imagesc(pCO2Grid(:,:,1))

%% 3b. Now pretty global maps of one month of each of SST and pCO2 data.
%I have provided example code for plotting January sea surface temperature
%(though you may need to make modifications based on differences in how you
%set up or named your variables above).

figure(1); clf
worldmap world
contourfm(latgrid, longrid, sstGrid(:,:,1)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January Sea Surface Temperature (^oC)')

%% Check that you can make a similar type of global map for another month
%and/or for pCO2 using this approach. Check the documentation and see
%whether you can modify features of this map such as the contouring
%interval, color of the contour lines, labels, etc.

figure(999); clf
worldmap world
contourfm(latgrid, longrid, sstGrid(:,:,7)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('July Sea Surface Temperature (^oC)')


%% 4. Calculate and plot a global map of annual mean pCO2
%<--
% Find the mean!
meanpCO2 = mean(pCO2Grid,3);
meanpsst = mean(sstGrid,3);
%disp(meanpCO2)

figure(2); clf
worldmap world
contourfm(latgrid, longrid, pCO2Grid(:,:,1)', 'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Annual Mean pCO2')

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2

% source = https://www.esrl.noaa.gov/gmd/ccgg/trends/gl_data.html
meanpCO22000 = 368.84;

% This takes the 2D matrix of annual pCO2 and subtracts the 2000 mean pCO2
mean_difference = meanpCO2' - meanpCO22000;

figure(3); clf
worldmap world
contourfm(latgrid, longrid, mean_difference, 'linecolor','none');
colorbar ;
geoshow('landareas.shp','FaceColor','black');
title('Difference between mean seawater pCO2 and atm pCO2')


%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle

% Normalize pCO2 around mean annual temp ---> pCO2 w/ biophysical effect
biopCO2 = pCO2Grid.*exp(0.0423*(repmat(meanpsst,[1 1 12])-sstGrid));

% Effect of Temperature on pCO2
    % calculating the amount of pCO2 change caused by temp change

temppCO2 = meanpCO2.*exp(0.0423*(sstGrid-repmat(meanpsst,[1 1 12])));

% T-B effects
% + means temp exceeds bio, - means bio exceeds temp

tb = mean(biopCO2-temppCO2);
tbAnnualMean = mean(biopCO2-temppCO2,3);

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

% BATS

% batsLoc = find(CO2data.LAT == 

% Station P


% Ross Sea


%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above




