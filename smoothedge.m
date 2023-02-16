%% smoothedge.m
% Defines constrained and loaded edges from a cell array of checkered
% boundaries and return a cell array with smoothed boundaries, except for
% constrained and loaded ones.
%% FUNCTION STRUCTURE
% USER DEFINED VARIABLES
% NODE DEFINITIONS
% SMOOTHING
function [ns,n2] = smoothedge(n)

global nely fixeddofs F

%% USER DEFINED VARIABLES
prompt   = {'Average span'};
def      = {'5'};
dlgTitle = 'Insert smoothing parameters';
lineNo   = 1;
answer   = inputdlg(prompt,dlgTitle,lineNo,def);        ii=1;
span     = sscanf(answer{ii},'%f');

ns = cell(size(n));
nodecon = [];
nodeload = [];


%% NODE DEFINITIONS
% Constraint node definition

for ii = 1:length(fixeddofs)
    
    nod  = round(fixeddofs(ii)/2);
    xnod = ceil(nod/(nely+1))-1;
    ynod = (nod-1)-xnod*(nely+1);
    nodecon = [nodecon; ynod xnod];  % [y_node x_node]
    
end

clear ii

% Load node defintion

nonzeros = find(F);

for ii = 1:size(nonzeros,1)
    
    nodf  = round(nonzeros(ii)/2);
    xnodf = ceil(nodf/(nely+1))-1;
    ynodf = (nodf-1)-xnodf*(nely+1);
    nodeload = [nodeload; ynodf xnodf];
    
end

% Passive areas (????????????????????????)

% Add 3rd column to n, 1=nosmooth 0=smooth, --> join smoothed vector with
% non smoothed

nosmooth = [nodecon; nodeload];

clear ii

for ii = 1:size(n,1)
    for jj = 1:size(nosmooth,1)
        n{ii,1}((n{ii,1}(:,1) == nosmooth(jj,1) & n{ii,1}(:,2) == nosmooth(jj,2)),3) = 1;
    end
end

n2 = n;

%% SMOOTHING

clear ii

for ii = 1:size(n,1)
    
    c1 = 2;  % first and last element index
    c2 = 3;
    
    ns2{ii,1}(1,1:2)   = n{ii,1}(1,1:2);
    %     ns{ii,1}(end,1:2) = n{ii,1}(end,1:2);
    
    while c2 <= size(n{ii,1},1)
        
        if n{ii,1}(c2,3) == 1 && n{ii,1}(c2-1,3) == 0
            
            smoy = smooth(n{ii,1}(c1-1:c2,1),span);
            smox = smooth(n{ii,1}(c1-1:c2,2),span);
            ns2{ii,1}(c1:c2-1,1) = smoy(2:end-1);    % Store smoothed
            ns2{ii,1}(c1:c2-1,2) = smox(2:end-1);
            c1 = c2;
            
        elseif n{ii,1}(c2,3) == 0 && n{ii,1}(c2-1,3) == 1
            
            ns2{ii,1}(c1:c2-1,1) = n{ii,1}(c1:c2-1,1);  % Store NON smoothed
            ns2{ii,1}(c1:c2-1,2) = n{ii,1}(c1:c2-1,2);
            c1 = c2;
            
        end
        
        c2 = c2+1;
        
    end
    
    if n{ii,1}(end,3) == 0
        
        ns2{ii,1}(c1:c2-1,1) = smooth(n{ii,1}(c1:c2-1,1),span);  % Store last elements smoothed
        ns2{ii,1}(c1:c2-1,2) = smooth(n{ii,1}(c1:c2-1,2),span);
        
    else
        
        ns2{ii,1}(c1:c2-1,1) = n{ii,1}(c1:c2-1,1);  % Store last elements NON smoothed
        ns2{ii,1}(c1:c2-1,2) = n{ii,1}(c1:c2-1,2);
        
    end
    
    ns2{ii,1}(:,3) = n{ii,1}(:,3);
    ns3{ii,1} = zeros(size(ns2{ii,1},1),3);
    
    if mod(size(ns2{ii,1},1),2) == 0
        ns3{ii,1}(1:floor(length(ns2{ii,1})/2),:) = ns2{ii,1}(floor(length(ns2{ii,1})/2)+1:end,:);
        ns3{ii,1}(floor(length(ns2{ii,1})/2)+1:end,:) = ns2{ii,1}(2:floor(length(ns2{ii,1})/2)+1,:);
    else
        ns3{ii,1}(1:floor(length(ns2{ii,1})/2),:) = ns2{ii,1}(floor(length(ns2{ii,1})/2)+2:end,:);
        ns3{ii,1}(floor(length(ns2{ii,1})/2)+1:end,:) = ns2{ii,1}(2:floor(length(ns2{ii,1})/2)+2,:);
    end
end


for ii = 1:size(n,1)
    
    c1 = 2;  % first and last element index
    c2 = 3;
    
    ns4{ii,1}(1,1:2) = ns3{ii,1}(1,1:2);
    %     ns{ii,1}(end,1:2) = n{ii,1}(end,1:2);
    
    while c2 <= size(ns3{ii,1},1)
        
        if ns3{ii,1}(c2,3) == 1 && ns3{ii,1}(c2-1,3) == 0
            
            smoy = smooth(ns3{ii,1}(c1-1:c2,1),span);
            smox = smooth(ns3{ii,1}(c1-1:c2,2),span);
            ns4{ii,1}(c1:c2-1,1) = smoy(2:end-1);    % Store smoothed
            ns4{ii,1}(c1:c2-1,2) = smox(2:end-1);
            c1 = c2;
            
        elseif ns3{ii,1}(c2,3) == 0 && ns3{ii,1}(c2-1,3) == 1
            
            ns4{ii,1}(c1:c2-1,1) = ns3{ii,1}(c1:c2-1,1);  % Store NON smoothed
            ns4{ii,1}(c1:c2-1,2) = ns3{ii,1}(c1:c2-1,2);
            c1 = c2;
            
        end
        
        c2 = c2+1;
        
    end
    
    if ns3{ii,1}(end,3) == 0
        
        ns4{ii,1}(c1:c2-1,1) = smooth(ns3{ii,1}(c1:c2-1,1),span);  % Store last elements smoothed
        ns4{ii,1}(c1:c2-1,2) = smooth(ns3{ii,1}(c1:c2-1,2),span);
        
    else
        
        ns4{ii,1}(c1:c2-1,1) = ns3{ii,1}(c1:c2-1,1);  % Store last elements NON smoothed
        ns4{ii,1}(c1:c2-1,2) = ns3{ii,1}(c1:c2-1,2);
        
    end
    ns4{ii,1}(:,3) = ns3{ii,1}(:,3);
    ns{ii,1} = zeros(size(ns4{ii,1},1),3);
    if mod(size(ns4{ii,1},1),2) == 0
        ns{ii,1}(1:floor(length(ns4{ii,1})/2),:) = ns4{ii,1}(floor(length(ns4{ii,1})/2):end-1,:);
        ns{ii,1}(floor(length(ns4{ii,1})/2)+1:end,:) = ns4{ii,1}(1:floor(length(ns4{ii,1})/2),:);
    else
        ns{ii,1}(1:floor(length(ns4{ii,1})/2)+1,:) = ns4{ii,1}(floor(length(ns4{ii,1})/2):end-1,:);
        ns{ii,1}(floor(length(ns4{ii,1})/2)+2:end,:) = ns4{ii,1}(1:floor(length(ns4{ii,1})/2),:);
    end
end


