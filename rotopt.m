%%rotopt.m
%THIS FUNCTION TAKES USER DEFINED INPUTS FOR INITAL,STEP & FINAL
%ORIENTATION AND ROTATES PART TO BE MANUFACTURED WHICH THEN ALLOWS USER TO
%SEE OPTIMAL ORIENTATION FOR MINIMAL SUPPORT STRUCTURE DURING ADDITIVE
%MANUFACTURE

function [alphaopt,minAreaInc] = rotopt(nodes)
global ncel m tmin OrigCompArea grey outflag indexx alpha minArea OptShape epsilon Areas Height Base

Areas    = [];
Height   = [];
Base     = []; 
OptShape = [];
alphaopt = [];

prompt     = {'Initial orientation [deg]' 'Final orientation [deg]' 'Angle increase step [deg]'};
def        = {'0' '180' '5'};
dlgTitle   = 'Building angle optimisation parameters';
lineNo     = 1;
answer     = inputdlg(prompt,dlgTitle,lineNo,def);        ii=1;
alphamin   = sscanf(answer{ii},'%f');                     ii=ii+1;
alphamax   = sscanf(answer{ii},'%f');                     ii=ii+1;
alphastep  = sscanf(answer{ii},'%f');

 for alpha = alphamin:alphastep:alphamax
    for ii = 1:ncel        
        for jj = 1:size(nodes{ii,1},1)
        nodesRot{ii,1}(jj,1) = nodes{ii,1}(jj,1)*cos(deg2rad(alpha))- nodes{ii,1}(jj,2)*sin(deg2rad(alpha));
        nodesRot{ii,1}(jj,2) = nodes{ii,1}(jj,2)*cos(deg2rad(alpha))+ nodes{ii,1}(jj,1)*sin(deg2rad(alpha));       
        end
    end
    
% keyboard
    
    %%% Rotated config plot
    [nodesRot]=gcolor(nodesRot);    
    
    %%% Support material optimization
    clear todo
    for ii = 2:ncel
        todo{ii-1,1} = nodesRot{ii,1};
    end
    finalshape   = cell(1,1);
    OrigArea     = zeros(ncel,1);
    OptArea      = zeros(ncel,1);
    OptVoidArea  = [];
    OrigVoidArea = [];
    indexx       = 1;       % Initialize index for the inner loop
    %minArea      = 0.3;     % Area control parameter
    
    %%% External support
    outflag  = 0;
    newshape = nodesRot{1,1};
    while outflag == 0
        [D,newcontour] = buildextsupport(newshape);
        if isempty(D)==0 && polyarea(D(:,2),D(:,1)) >= minArea  % is NOT empty and Area is big enough
            todo = [todo; D];
        end
        newshape = newcontour;
    end
    finalshape{1,1} = newshape;

    %%% Internal Support
    while sum(cellfun('isempty',todo))~=size(cellfun('isempty',todo),1)  % Internal support building (+ points/area control)
        [D1,D2,O,nonfeasdom] = buildsupport_EVO(todo{indexx,1});
        todo{indexx,1} = [];
        % if no unfeasible points
        if isempty(nonfeasdom)==1
            if polyarea(O(:,2),O(:,1)) >= minArea
                finalshape = [finalshape; O(:,1:2)];
            end
        else
            % O LOOP
            if polyarea(O(:,2),O(:,1)) >= minArea && (sum(ismember(O(:,1:2),nonfeasdom(:,1:2),'rows'))==0)    % O big enough and without unfeasible points --> OK
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
        indexx = indexx+1;
    end
    
    %%% Area assessment
    for jj = 2:size(finalshape,1)
        OptVoidArea(jj-1) = polyarea(finalshape{jj,1}(:,2),finalshape{jj,1}(:,1));
    end
    SumOptVoid   = sum(OptVoidArea);
    OptCompArea  = polyarea(finalshape{1,1}(:,2),finalshape{1,1}(:,1))-SumOptVoid;
    AreaIncrease = (OptCompArea-OrigCompArea)/OrigCompArea*100;
    Areas        = [Areas; AreaIncrease];
    if AreaIncrease<=min(Areas)
        OptShape   = finalshape;
        alphaopt   = alpha;
        minAreaInc = AreaIncrease;
    end
    %%% Delta y and figure "basis" for plot
    Height = [Height; abs(max(finalshape{1,1}(:,1))-min(finalshape{1,1}(:,1)))];
    basew  = finalshape{1,1}(finalshape{1,1}(:,1)>max(finalshape{1,1}(:,1))-epsilon,2);
    Base   = [Base; max(basew)-min(basew)];
 end

%%% Delete overlapped points --> intersection too close to the original
%%% point (useless?)
for ii = 1:size(OptShape,1)
    if abs(OptShape{ii,1}(end-1,2)-OptShape{ii,1}(end,2))<epsilon && abs(OptShape{ii,1}(end-1,1)-OptShape{ii,1}(end,1))<epsilon
        OptShape{ii,1}=[OptShape{ii,1}(1:end-2,:);OptShape{ii,1}(end,:)];
    end
end

%%% Final plot
figure
hold on
axis image
for ii = 1:size(OptShape,1)
    plot(OptShape{ii,1}(:,2),-OptShape{ii,1}(:,1),'k')
end

%%% Filled final shape
figure
plot(OptShape{1,1}(:,2),-OptShape{1,1}(:,1),'k');
hold on
axis image
fill(OptShape{1,1}(:,2),-OptShape{1,1}(:,1),grey)
for jj = 2:size(OptShape,1)
    fill(OptShape{jj,1}(:,2),-OptShape{jj,1}(:,1),[1 1 1])
    plot(OptShape{jj,1}(:,2),-OptShape{jj,1}(:,1),'k')
end