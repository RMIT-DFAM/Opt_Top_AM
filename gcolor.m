%% gcolor.m
% THIS FUNCTION CREATES GRADIENT PLOTS FOR STRUCTURE BASED ON USER DEFINED
% FEASABLE ANGLE OF MANUFACTURE

function [ns] = gcolor(cs)
% Color gradient of cs cell array
global red blue orange green epsilon mdeg

ymax = max(cs{1,1}(:,1));
for ii = 1:size(cs,1) 
    
%     cs{ii,1}(1,3) = atan(-(cs{ii,1}(2,1)-cs{ii,1}(end-1,1))/(cs{ii,1}(2,2)-cs{ii,1}(end-1,2)));
%         if cs{ii,1}(1,3) < 0
%             cs{ii,1}(1,3) = cs{ii,1}(1,3)+2*pi;
%         end
%     cs{ii,1}(end,3) = cs{ii,1}(1,3);

    cs{ii,1}=[cs{ii,1}(end-1,:);cs{ii,1};cs{ii,1}(2,:)];
    for jj = 2:size(cs{ii,1},1)-1
        if cs{ii,1}(jj+1,1) <= cs{ii,1}(jj-1,1) && cs{ii,1}(jj+1,2) >= cs{ii,1}(jj-1,2)                          % y2<y1 & x2>x1
            cs{ii,1}(jj,3) = atan(-(cs{ii,1}(jj+1,1)-cs{ii,1}(jj-1,1))/(cs{ii,1}(jj+1,2)-cs{ii,1}(jj-1,2)));
        elseif cs{ii,1}(jj+1,1) >= cs{ii,1}(jj-1,1) && cs{ii,1}(jj+1,2) >= cs{ii,1}(jj-1,2)
            cs{ii,1}(jj,3) = atan(-(cs{ii,1}(jj+1,1)-cs{ii,1}(jj-1,1))/(cs{ii,1}(jj+1,2)-cs{ii,1}(jj-1,2)))+2*pi;
        elseif cs{ii,1}(jj+1,1) >= cs{ii,1}(jj-1,1) && cs{ii,1}(jj+1,2) <= cs{ii,1}(jj-1,2)
            cs{ii,1}(jj,3) = atan(-(cs{ii,1}(jj+1,1)-cs{ii,1}(jj-1,1))/(cs{ii,1}(jj+1,2)-cs{ii,1}(jj-1,2)))+pi;
        else
            cs{ii,1}(jj,3) = atan(-(cs{ii,1}(jj+1,1)-cs{ii,1}(jj-1,1))/(cs{ii,1}(jj+1,2)-cs{ii,1}(jj-1,2)))+pi;
        end
    end
    cs{ii,1}=cs{ii,1}(2:end-1,:);
    cs{ii,1}(:,3) = round(rad2deg(cs{ii,1}(:,3))); % ROUNDED to avoid steps on gradient values
end

clear ii jj

% Fourth column binary for feasible or unfeasible
cs{1,1}(:,4) = zeros(size(cs{1,1},1),1);
for kk = 1:size(cs{1,1},1)
    if cs{1,1}(kk,3)>180-mdeg && cs{1,1}(kk,3)<180+mdeg && cs{1,1}(kk,1)<ymax-epsilon
        cs{1,1}(kk,4) = 1;
    end
end

for ii = 2:size(cs,1)  
    cs{ii,1}(:,4) = zeros(size(cs{ii,1},1),1);
    for jj = 1:size(cs{ii,1},1)
        if ((cs{ii,1}(jj,3)>=0 && cs{ii,1}(jj,3)<mdeg)) || ((cs{ii,1}(jj,3)>360-mdeg && cs{ii,1}(jj,3)<=360))
            cs{ii,1}(jj,4) = 1;
        end
    end
end        
        
clear ii jj

%%% GRADIENT PLOT
figure
for ii = 2:size(cs,1)
    % Internal Green points (feasible)
    if sum((cs{ii,1}(:,3)>mdeg+10 & cs{ii,1}(:,3)<360-mdeg-10))>0
        plot(cs{ii,1}((cs{ii,1}(:,3)>mdeg+10 & cs{ii,1}(:,3)<=360-mdeg-10),2),-cs{ii,1}((cs{ii,1}(:,3)>mdeg+10 & cs{ii,1}(:,3)<=360-mdeg-10),1),'.','MarkerSize',15,'Color',green);
        hold all
    end
    % Internal Red points (infeasible)
    if sum((cs{ii,1}(:,3)>=0 & cs{ii,1}(:,3)<=mdeg) | (cs{ii,1}(:,3)>=360-mdeg & cs{ii,1}(:,3)<=360))>0
        plot(cs{ii,1}((cs{ii,1}(:,3)>=0 & cs{ii,1}(:,3)<=mdeg) | (cs{ii,1}(:,3)>=360-mdeg & cs{ii,1}(:,3)<=360),2),-cs{ii,1}((cs{ii,1}(:,3)>=0 & cs{ii,1}(:,3)<=mdeg) | (cs{ii,1}(:,3)>=360-mdeg & cs{ii,1}(:,3)<=360),1),'.','MarkerSize',15,'Color',red);
    end
    % Internal Orange points (almost infeasible)
    if sum((cs{ii,1}(:,3)>mdeg & cs{ii,1}(:,3)<=mdeg+10) | (cs{ii,1}(:,3)>=360-mdeg-10 & cs{ii,1}(:,3)<360-mdeg))>0
        plot(cs{ii,1}((cs{ii,1}(:,3)>mdeg & cs{ii,1}(:,3)<=mdeg+10) | (cs{ii,1}(:,3)>=360-mdeg-10 & cs{ii,1}(:,3)<360-mdeg),2),-cs{ii,1}((cs{ii,1}(:,3)>mdeg & cs{ii,1}(:,3)<=mdeg+10) | (cs{ii,1}(:,3)>=360-mdeg-10 & cs{ii,1}(:,3)<360-mdeg),1),'.','MarkerSize',15,'Color',orange);
    end   
end
plot(cs{1,1}(:,2),-cs{1,1}(:,1),'LineWidth',2,'Color',blue);
% External Red points (infeasible)
% if sum(cs{1,1}(:,3)>=180-mdeg+1 & cs{1,1}(:,3)<=180+mdeg-1)>0
%     plot(cs{1,1}(cs{1,1}(:,3)>=180-mdeg+1 & cs{1,1}(:,3)<=180+mdeg-1 & cs{1,1}(:,1)<ymax-epsilon,2),-cs{1,1}(cs{1,1}(:,3)>=180-mdeg+1 & cs{1,1}(:,3)<=180+mdeg-1 & cs{1,1}(:,1)<ymax-epsilon,1),'.','MarkerSize',15,'Color',red);
% end

plot(cs{1,1}(cs{1,1}(:,4)==1,2),-cs{1,1}(cs{1,1}(:,4)==1,1),'.','MarkerSize',15,'Color',red)

% for ii = 1:size(cs{1,1},1) %%% Ext red dots plot
%     if cs{1,1}(ii,4)==1 && sum(cs{1,1}(ii-1,4) & cs{1,1}(ii+1,4))>0    
%         plot(cs{1,1}(ii,2),-cs{1,1}(ii,1),'.','MarkerSize',15,'Color',red);
%     end
% end

% External Orange points (almost infeasible)
if sum((cs{1,1}(:,3)>=180-mdeg-10 & cs{1,1}(:,3)<180-mdeg) | (cs{1,1}(:,3)>180+mdeg & cs{1,1}(:,3)<=180+mdeg+10))>0
    plot(cs{1,1}((cs{1,1}(:,3)>=180-mdeg-10 & cs{1,1}(:,3)<180-mdeg) | (cs{1,1}(:,3)>180+mdeg & cs{1,1}(:,3)<=180+mdeg+10),2),-cs{1,1}((cs{1,1}(:,3)>=180-mdeg-10 & cs{1,1}(:,3)<180-mdeg) | (cs{1,1}(:,3)>180+mdeg & cs{1,1}(:,3)<=180+mdeg+10),1),'.','MarkerSize',15,'Color',orange);
end
axis image
legend ('Buildable','Unbuildable','Almost buildable')

ns = cs;





