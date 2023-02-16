%% topmain.m 
%THIS SCRIPT BUILDS A PLANAR OPTIMAL GEOMETRY AND MODIFIES IT TO BE
%MANUFACTURABLE WITH ADDITIVE METHODS
%IT AND IS DESCRIBED IN THE PUBLICATION:
%Optimal topology for additive manufacture: A method for enabling additive manufacture of support-free optimal structures
%Int J. Materials and Design 63 (2014) 678–690

%% PROGRAM STRUCTURE
% HOUSEKEEPING
% GLOBAL VARIABLES
% USER DEFINED VARIABLES
% CONVERSION OF OPTIMISED TOPOLOGY TO BINARY AND EXTRACTION OF BOUNDARY NODES
% SUPPORT MATERIAL OPTIMISATION
% INTERNAL SUPPORT
% OPTIMAL ORIENTATION FOR MINIMAL SUPPORT AREA
% EXPORT TXT FILE FROM OPTIMAL SHAPE
%% HOUSEKEEPING 
clear all 
close all
clc
%% GLOBAL VARIABLES
global nelx nely volfrac penal rmin binaryElem th ft x xPhys passive fixeddofs freedofs F ...
       grey red blue orange green index nodestruc ncel tmin m mdeg outflag OrigCompArea indexx alpha ...
       epsilon finalshape minArea contExt OptShape pArea Areas Height Base

epsilon = 1e-5;
grey    = [0.7 0.7 0.7];
red     = [230 0 10]./255;
blue    = [8 18 149]./255;
orange  = [255 130 0]./255;
green   = [51 204 0]./255;

%% USER DEFINED VARIABLES
prompt   = {'N° elem x','N° elem y','Volume fraction','Penalty coefficient','Min radius','Threshold','Filter type'}; %Dialog Prompt
def      = {'100','40','0.5','3','1.5','0.5','1'};   %Default Values in text box
dlgTitle = 'Insert optimization problem parameters'; %Prompt box title
lineNo   = 1; %Text Box Height
answer   = inputdlg(prompt,dlgTitle,lineNo,def);        ii=1;    %Displays Prompt Box
nelx     = sscanf(answer{ii},'%f');                     ii=ii+1; %Reads values for Number of Elements in Horizontal Direction inputted in Prompt Box
nely     = sscanf(answer{ii},'%f');                     ii=ii+1; %Reads values for Number of Elements in Horizontal Direction inputted in Prompt Box
volfrac  = sscanf(answer{ii},'%f');                     ii=ii+1; %Volume Fraction (target volume/initial volume)
penal    = sscanf(answer{ii},'%f');                     ii=ii+1; %Penalisation power
rmin     = sscanf(answer{ii},'%f');                     ii=ii+1; %Filter radius (Radius of filter sensitivity area)
th       = sscanf(answer{ii},'%f');                     ii=ii+1; %Density Threshold value used to set conditions for turning voxels density from either 0 or 1 
ft       = sscanf(answer{ii},'%f');                              %Filter type

quitflag  = 1; %Creates constant to disable dialog box
while quitflag
    choice = menu('Menu','99 line code','88 line code','Exit'); %Creates Menu for choice between 99 line code or 88 line code
    switch choice
        case 1   % 99 line code - SLOWER!
%           topall(nelx,nely,volfrac,penal,rmin)
            topCantiDist(nelx,nely,volfrac,penal,rmin) % Runs A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND function
            quitflag = 0;
        case 2   % 88 line code
            top88(nelx,nely,volfrac,penal,rmin,ft) % Runs Efficient topology optimization in MATLAB using 88 lines of code function
            quitflag = 0;
        case 3 % Exit
            quitflag = 0;
    end
end

%% CONVERSION OF OPTIMISED TOPOLOGY TO BINARY AND EXTRACTION OF BOUNDARY NODES
% Conversion to Binary
binaryElem = ConvertToBinary(x,th); % Converts voxel densitys (x in global variables) to either 1 or 0 (binary) based off userdefined threshold (th in global variables)
% Extraction of boundary nodes
bound = bwboundaries(binaryElem); % goes from bottom clockwise and gives x-y grid co-ordinates of boundaRY
ncel = size(bound,1); %Returns how many rows in bound which correspond to how many clockwise rotations (i.e 3x3 grid will have 2 rotations)
n = cell(ncel,1); %Creates cell named n with ncel rows and 1 column

for ii = 1:ncel
    n{ii,1} = nodes(bound{ii,1},ii); %Starting from bottom left moves around clockwise and gives each nodal co-ordinate on boundary referenced from grid
end
% Smoothing
[nodestruc,n] = smoothedge(n);
% Contour plot
figure
for jj = 1:ncel
    plot(n{jj,1}(:,2),-n{jj,1}(:,1),'.-','Color',grey)
    hold on
    plot(nodestruc{jj,1}(:,2),-nodestruc{jj,1}(:,1),'r','LineWidth',2)
    axis image
end

%% SUPPORT MATERIAL OPTIMISATION
prompt   = {'Feasible angle [deg]' 'tmin [dots]' 'Area penalty coefficient [-]'}; % Area penalty coefficient is related to the successive evaluation of
def      = {'45' '1' '0.3'};
dlgTitle = 'Insert support material optimisation parameters';
lineNo   = 1;
answer   = inputdlg(prompt,dlgTitle,lineNo,def);        ii=1;
mdeg     = sscanf(answer{ii},'%f');                     ii=ii+1;
tmin     = sscanf(answer{ii},'%f');                     ii=ii+1;
pArea    = sscanf(answer{ii},'%f'); 

m    = tan(deg2rad(mdeg));
ymax = max(nodestruc{1,1}(:,1));       % bottom

% Gradient assessment and plot
[ns]=gcolor(nodestruc);

for ii = 2:ncel
    todo{ii-1,1} = ns{ii,1}(:,:);
end

finalshape = cell(1,1);
OrigArea   = zeros(ncel,1);
OptArea    = zeros(ncel,1);
index      = 1;             % Internal loop index
minArea    = round(0.001*(nelx*nely));     % Area control parameter

% External support
contExt  = 0;
outflag  = 0;
newshape = ns{1,1};
while outflag == 0
    [D,newcontour] = buildextsupport(newshape);
    if isempty(D)==0 && polyarea(D(:,2),D(:,1)) >= minArea  % is NOT empty and Area is big enough
        todo = [todo; D];
    end
    newshape = newcontour;
    contExt = contExt+1;
end 
finalshape{1,1} = newshape;

%% INTERNAL SUPPORT
while sum(cellfun('isempty',todo))~=size(cellfun('isempty',todo),1)  % Internal support building (+ points/area control)
    
    [D1,D2,O,nonfeasdom] = buildsupport_EVO(todo{index,1});
    todo{index,1} = [];
    % if there are no unfeasible points --> O goes directly to finalshape
    if isempty(nonfeasdom)==1
        if polyarea(O(:,2),O(:,1)) >= minArea
            finalshape = [finalshape; O(:,1:2)];
        end
    else
        % O-LOOP
        if polyarea(O(:,2),O(:,1)) >= minArea && (sum(ismember(O(:,1:2),nonfeasdom(:,1:2),'rows'))==0)  % O big enough and without unfeasible points --> OK
            finalshape = [finalshape; O(:,1:2)];
        elseif polyarea(O(:,2),O(:,1)) >= minArea && (sum(ismember(O(:,1:2),nonfeasdom(:,1:2),'rows'))>0) % O big enough but contains unfeas points --> to do
            todo = [todo ; O];
        end
        % D1 LOOP
        if (isempty(D1) == 0) && (sum(ismember(D1(:,1:2),nonfeasdom(:,1:2),'rows'))>0) 
            todo = [todo; D1];
        elseif (isempty(D1) == 0) && (sum(ismember(D1(:,1:2),nonfeasdom(:,1:2),'rows'))==0) && polyarea(D1(:,2),D1(:,1)) >= minArea
            finalshape = [finalshape; D1(:,1:2)];
        end
        % D2 LOOP
        if (isempty(D2) == 0) && (sum(ismember(D2(:,1:2),nonfeasdom(:,1:2),'rows'))>0)
            todo = [todo; D2];
        elseif (isempty(D2) == 0) && (sum(ismember(D2(:,1:2),nonfeasdom(:,1:2),'rows'))==0) && polyarea(D2(:,2),D2(:,1)) >= minArea
            finalshape = [finalshape; D2(:,1:2)];
        end 
    end
    index = index+1;
end
clear todo
% Delete overlapped points --> intersection too close to the point (Useless??)
for ii = 1:size(finalshape,1)
    if abs(finalshape{ii,1}(end-1,2)-finalshape{ii,1}(end,2))<epsilon && abs(finalshape{ii,1}(end-1,1)-finalshape{ii,1}(end,1))<epsilon
        finalshape{ii,1}=[finalshape{ii,1}(1:end-2,:);finalshape{ii,1}(end,:)];
    end
end

% External contours and area assessment
figure
plot(finalshape{1,1}(:,2),-finalshape{1,1}(:,1),'k')
axis image
hold on
for jj = 2:size(finalshape,1)
    plot(finalshape{jj,1}(:,2),-finalshape{jj,1}(:,1),'k')
    OptVoidArea(jj-1) = polyarea(finalshape{jj,1}(:,2),finalshape{jj,1}(:,1));
end

for ii = 2:ncel 
    OrigVoidArea(ii) = polyarea(ns{ii,1}(:,2),ns{ii,1}(:,1)); 
end

SumOrigVoid  = sum(OrigVoidArea);
SumOptVoid   = sum(OptVoidArea);
OrigCompArea = polyarea(ns{1,1}(:,2),ns{1,1}(:,1))-SumOrigVoid;
OptCompArea  = polyarea(finalshape{1,1}(:,2),finalshape{1,1}(:,1))-SumOptVoid;
AreaIncrease = (OptCompArea-OrigCompArea)/OrigCompArea*100;

disp (['The increase in area compared to the original is: ', num2str(AreaIncrease), '%'])

% Buildable contours and unfeasible points FIGURE
figure
plot(finalshape{1,1}(:,2),-finalshape{1,1}(:,1),'Color',grey);
hold on
axis image
if sum(ns{1,1}(:,3)>=180-mdeg+1 & ns{1,1}(:,3)<=180+mdeg-1)>0  % +-1 just for edge points
    plot(ns{1,1}(ns{1,1}(:,3)>=180-mdeg+1 & ns{1,1}(:,3)<=180+mdeg-1 & ns{1,1}(:,1)<ymax,2),-ns{1,1}(ns{1,1}(:,3)>=180-mdeg+1 & ns{1,1}(:,3)<=180+mdeg-1 & ns{1,1}(:,1)<ymax,1),'.','MarkerSize',15,'Color',red);
end
if sum((ns{1,1}(:,3)>=180-mdeg-10 & ns{1,1}(:,3)<180-mdeg) | (ns{1,1}(:,3)>180+mdeg & ns{1,1}(:,3)<=180+mdeg+10))>0
    plot(ns{1,1}((ns{1,1}(:,3)>=180-mdeg-10 & ns{1,1}(:,3)<180-mdeg) | (ns{1,1}(:,3)>180+mdeg & ns{1,1}(:,3)<=180+mdeg+10),2),-ns{1,1}((ns{1,1}(:,3)>=180-mdeg-10 & ns{1,1}(:,3)<180-mdeg) | (ns{1,1}(:,3)>180+mdeg & ns{1,1}(:,3)<=180+mdeg+10),1),'.','MarkerSize',15,'Color',orange);
end
legend 'Buildable contours' 'Infeasible region' 'Almost infeasible region'
for jj = 2:size(finalshape,1)
    plot(finalshape{jj,1}(:,2),-finalshape{jj,1}(:,1),'Color',grey)
end
for ii = 2:size(ns,1)
    if sum((ns{ii,1}(:,3)>=0 & ns{ii,1}(:,3)<=mdeg) | (ns{ii,1}(:,3)>=360-mdeg & ns{ii,1}(:,3)<=360))>0
        plot(ns{ii,1}((ns{ii,1}(:,3)>=0 & ns{ii,1}(:,3)<=mdeg) | (ns{ii,1}(:,3)>=360-mdeg & ns{ii,1}(:,3)<=360),2),-ns{ii,1}((ns{ii,1}(:,3)>=0 & ns{ii,1}(:,3)<=mdeg) | (ns{ii,1}(:,3)>=360-mdeg & ns{ii,1}(:,3)<=360),1),'.','MarkerSize',15,'Color',red);
    end
    if sum((ns{ii,1}(:,3)>mdeg & ns{ii,1}(:,3)<=mdeg+10) | (ns{ii,1}(:,3)>=360-mdeg-10 & ns{ii,1}(:,3)<360-mdeg))>0
        plot(ns{ii,1}((ns{ii,1}(:,3)>mdeg & ns{ii,1}(:,3)<=mdeg+10) | (ns{ii,1}(:,3)>=360-mdeg-10 & ns{ii,1}(:,3)<360-mdeg),2),-ns{ii,1}((ns{ii,1}(:,3)>mdeg & ns{ii,1}(:,3)<=mdeg+10) | (ns{ii,1}(:,3)>=360-mdeg-10 & ns{ii,1}(:,3)<360-mdeg),1),'.','MarkerSize',15,'Color',orange);
    end   
end

% Filled final shape FIGURE
figure
plot(finalshape{1,1}(:,2),-finalshape{1,1}(:,1),'k');
hold on
axis image
fill(finalshape{1,1}(:,2),-finalshape{1,1}(:,1),grey)
for jj = 2:size(finalshape,1)
    fill(finalshape{jj,1}(:,2),-finalshape{jj,1}(:,1),[1 1 1])
    plot(finalshape{jj,1}(:,2),-finalshape{jj,1}(:,1),'k')
end

%% OPTIMAL ORIENTATION FOR MINIMAL SUPPORT AREA
alpha = [];
[alphaopt,Areasupp] = rotopt(nodestruc);

disp (['The increase in area compared to the original is: ', num2str(Areasupp), '%'])
disp (['The optimal inclination for manufacturing is: ', num2str(alphaopt), '°'])   

% Area-Height-Base plot
figure
plot([1:2:2*size(Areas,1)],Areas)
addaxis([1:2:2*size(Height,1)],Height,'.-')
addaxis([1:2:2*size(Base,1)],Base,'--')
addaxislabel(1,'- Increase in area [%]');
addaxislabel(2,'.- Max height [dots]');
addaxislabel(3,'-- Base width [dots]');

%% EXPORT TXT FILE FROM OPTIMAL SHAPE

xmin = min(OptShape{1,1}(:,2));
xmax = max(OptShape{1,1}(:,2));
ymin = min(-OptShape{1,1}(:,1));
ymax = max(-OptShape{1,1}(:,1));
xdim = abs(max(OptShape{1,1}(:,2))-min(OptShape{1,1}(:,2)));
ydim = abs(max(OptShape{1,1}(:,1))-min(OptShape{1,1}(:,1)));

figure
plot(finalshape{1,1}(:,2),-finalshape{1,1}(:,1),'k')
axis image
axis([0-0.01*nelx nelx+0.01*nelx -nely-0.01*nely 0+0.01*nely])
hold on
for jj = 2:size(finalshape,1)
    plot(finalshape{jj,1}(:,2),-finalshape{jj,1}(:,1),'k')
end
axis off

figure
plot(OptShape{1,1}(:,2),-OptShape{1,1}(:,1),'k')
axis image
axis([xmin-0.01*xdim xmax+0.01*xdim ymin-0.01*ydim ymax+0.01*ydim])
hold on
for jj = 2:size(OptShape,1)
    plot(OptShape{jj,1}(:,2),-OptShape{jj,1}(:,1),'k')
end
axis off

quitflagexp  = 1;

while quitflagexp
    expchoice = menu('Export file','IGES format','Text format','Exit');
    switch expchoice
        case 1 %%% Export finalshape and OptShape as IGES files
            finalshape{1,1}=finalshape{1,1}(:,1:2);
            for ii=1:size(finalshape,1)
                finalshape{ii,1}(:,3) = 0;
            end
            OptShape{1,1}=OptShape{1,1}(:,1:2);
            for ii=1:size(OptShape,1)
                OptShape{ii,1}(:,3) = 0;
            end
            igesout(finalshape,'Default45t1')
            igesout(OptShape,'Opt45t1')
        
        case 2   %%% Coordinates exported as .txt file
            for jj = 1:size(OptShape,1)
               OptShape{jj,1}(:,3)=0;
               filename = (['shape',num2str(jj),'.txt']);
               dlmwrite(filename,1);
               xlswrite(filename,OptShape{jj,1}(:,1:3));
               comma2point(filename)
            end
        case 3 % Exit
          quitflagexp = 0;
    end
end

%% Save images: 
% print '-r600' '-djpeg' 'nomefile'
% scrsz = get(0,'ScreenSize')       % screen size