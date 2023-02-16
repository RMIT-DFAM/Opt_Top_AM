function [D3,D4,O,infeasibleVec]=buildsupport_EVO(dom)
global minA grey index ncel ns m tmin indexx alpha nodesRot finalshape pArea
% Inside-cell version
% dom = [y_smooth x_smooth grad feas/nonfeas]
% nonfeasdom = [y x q1 q2 y_intersection x_intersection]
%%% EVO VERSION %%%
% Problems to fix:

% if alpha==15 && indexx==129 
%     keyboard
% end

clear xmindom
xmindom = min(dom(:,2));
row  = size(dom,1);
c1   = 0;
c2   = 0;
jj   = 1;
nonfeasdom  = [];
unfeasibles = [];

%%% Remove overlapped points
delta    = 0.01;
dom      = [zeros(1,4);dom;zeros(1,4)];
dom(:,5) = 0;
for idx = 2:size(dom,1)-1
   if dom(idx,1)<dom(idx-1,1)+delta && dom(idx,1)>dom(idx-1,1)-delta && dom(idx,2)<dom(idx-1,2)+delta && dom(idx,2)>dom(idx-1,2)-delta && dom(idx,4)==1 && dom(idx-1,4)==0
       dom(idx,5)=1;
   elseif dom(idx,1)<dom(idx+1,1)+delta && dom(idx,1)>dom(idx+1,1)-delta && dom(idx,2)<dom(idx+1,2)+delta && dom(idx,2)>dom(idx+1,2)-delta && dom(idx,4)==1 && dom(idx+1,4)==0
       dom(idx,5)=1;
   end
end
if isempty(dom(:,5))==0
    dom(dom(:,5)==1,:)=[];
end

dom = dom(2:end-1,1:4);

%%% Reorder dom
pointer = find(dom(:,1)==min(dom(dom(:,2)==xmindom,1)),1,'first');
if pointer==1
    dom = dom;
else
    dom = [dom(pointer:end,:) ; dom(2:pointer,:)];
end
%%%

%%% Cleaning infeasible vertex
if sum(dom(:,4))==1 && sum(dom(:,1)==min(dom(:,1)))==1 && dom(dom(:,4)==1,1)==min(dom(:,1)) % there's 1 inf point and is a vertex
    dom(dom(:,4)==1,4)=0;
end

%%% Infeasible region identification
if sum(dom(:,4))>0 && sum(dom(:,4))<3 && row<=6   % If only 1-2 infeas points and domain dimension < 6
    unfeasibles{1,1} = dom(dom(:,4)==1,1:2);
else
    
    % NEW INFEASIBLE REGION IDENTIFICATION LOOP
    domzeros = [zeros(1,4);dom;zeros(1,4)];
    for ii = 2:size(domzeros,1)-2
        if domzeros(ii,4)==1 && domzeros(ii-1,4)==0 && domzeros(ii+1,4)==0 % delete single points
           domzeros(ii,4)=0; 
        end
        if domzeros(ii,4)==1 && domzeros(ii-1,4)==0 && domzeros(ii+1,4)==1 && domzeros(ii+2,4)==1
            c1 = ii;
        end
        if domzeros(ii,4)==1 && domzeros(ii-1,4)==1 && domzeros(ii-2,4)==1 && domzeros(ii+1,4)==0
            c2 = ii;
            unfeasibles{jj,1} = domzeros(c1:c2,1:2);
            jj = jj+1;
        end
    end
    dom = domzeros(2:end-1,:);

end

if isempty(unfeasibles) == 1 % modded: nonfeasdom --> unfeasibles
    D3 = [];
    D4 = [];
    O = dom;
    infeasibleVec = [];
else
    % New definition of nonfeasdom
    if sum(dom(:,4))<3 && row<=6   % Little domains with no more than 2 inf points
        nonfeasdom = unfeasibles{1,1};
    else                           % Normal domains: start from higher regions
%         yaverages = zeros(1,size(unfeasibles,1)); %
%         for ii = 1:size(unfeasibles,1)
%             yaverages(ii) = mean(unfeasibles{ii,1}(:,1));
%         end
        highest = zeros(1,size(unfeasibles,1));
        for ii = 1:size(unfeasibles,1)
            highest(ii) = min(unfeasibles{ii,1}(:,1));
        end
        nonfeasdom = unfeasibles{highest==min(highest),1}(:,1:2);
    end
    
    % Optimal Area of the feasible domains (V-shape)
    D1   = [];
    D2   = [];
    d1   = [];
    d2   = [];
    O    = [];
    Area = [];
    fitnessVec = [];
    doublens = [dom; dom(2:end,:)];

    for ii = 1:size(nonfeasdom,1)

        nonfeasdom(ii,3) = nonfeasdom(ii,1)+m*nonfeasdom(ii,2); % q1 for line 1 (left-side)
        nonfeasdom(ii,4) = nonfeasdom(ii,1)-m*nonfeasdom(ii,2); % q2 for line 2 (right-side)

        % Create left side domain (D1) CCW    
        ind1 = find(ismember(doublens(1:end-1,1:2),[nonfeasdom(ii,1) nonfeasdom(ii,2)],'rows'),1,'last')-1;
        ind11 = ind1;
        x1 = doublens(ind1,2);
        y1 = doublens(ind1,1);
        
        while x1<((-y1+nonfeasdom(ii,3))/m) && x1<=nonfeasdom(ii,2) % && y1>=nonfeasdom(ii,1)
            ind11 = ind11-1;
            x1    = doublens(ind11,2); % First point over the line
            y1    = doublens(ind11,1);
        end

        mp1 = (y1 - doublens(ind11+1,1))/(x1 - doublens(ind11+1,2)); % segment coefficient
        if abs(mp1)~=Inf 
            qp1 = y1-mp1*x1; 
            xy1 = [1 m;1 -mp1]\[nonfeasdom(ii,3); qp1]; % A\b
            yp1 = xy1(1); % domain intersection point coordinates
            xp1 = xy1(2);
        else
            xp1 = x1;
            yp1 = -m*xp1+nonfeasdom(ii,3);
        end

        d1  = [yp1 xp1 0 0; doublens(ind11+1:ind1+1,:); yp1 xp1 0 0];

        % Create right side domain (D2) CW
        ind2 = find(ismember(doublens(:,1:2),[nonfeasdom(ii,1) nonfeasdom(ii,2)],'rows'),1,'first')+1;
        ind22 = ind2;
        x2 = doublens(ind2,2);
        y2 = doublens(ind2,1);

        while x2>(y2-nonfeasdom(ii,4))/m && x2>=nonfeasdom(ii,2) && ind22<size(doublens,1) % !!!! ERRORE FINIRE DA QUI (ind22<size..)!!! && y2>=nonfeasdom(ii,1) 
            ind22 = ind22+1;
            x2    = doublens(ind22,2); % First point over the line
            y2    = doublens(ind22,1);
        end

        mp2  = (y2 - doublens(ind22-1,1))/(x2 - doublens(ind22-1,2)); % segment coefficient
        if abs(mp2) ~= Inf
            qp2 = y2-mp2*x2; 
            xy2 = [1 -m;1 -mp2]\[nonfeasdom(ii,4); qp2]; % A\b
            yp2 = xy2(1); % domain intersection point coordinates
            xp2 = xy2(2);
        else
            xp2 = x2;
            yp2 = m*xp2+nonfeasdom(ii,4);
        end

        d2  = [doublens(ind2-1:ind22-1,:); yp2 xp2 0 0; doublens(ind2-1,:)];   

        A1   = polyarea(d1(:,2),d1(:,1));
        A2   = polyarea(d2(:,2),d2(:,1));
        A    = A1+A2;
        Area = [Area; A];
        minA = min(Area);

        infeasibleVec = dom(dom(:,4)==1,1:2); % all infeasible points 
        Coverage      = (sum(ismember(d1(:,1:2),infeasibleVec,'rows'))+sum(ismember(d2(:,1:2),infeasibleVec,'rows'))); % Coverage of infeasible regions [0-1]
        maxCoverage   = size(infeasibleVec,1);
        
%         pArea      = 0.3;
        fitness    = pArea/A+(Coverage/maxCoverage)/A*(1-pArea); %% ?
        fitnessVec = [fitnessVec; fitness];
        fitnessMax = max(fitnessVec);
        
        if fitness>=fitnessMax % Optimal outputs !!!NEXT STEP!!! Implement penalty for non covering all infeas region
            D1 = d1;
            qd1 = nonfeasdom(ii,3);
            D2 = d2;
            qd2 = nonfeasdom(ii,4);
            O  = [yp1 xp1 0 0; doublens(ind2-1,:); yp2 xp2 0 0; doublens(ind22:ind11,:); yp1 xp1 0 0];
        end       
    end

    % Generete offset lines (and volume D3 or overwrite vol D1??)

    qoff1 = qd1 - tmin/cos(atan(-m)); % y+mx
    qoff2 = qd2 - tmin/cos(atan(m));

    yinf1 = max(D1(:,1))+1;         % y = mx+q   !!! Ho inserito +1 per prendere porzione retta offset da intersecare esterna al dominio
    xinf1 = (yinf1 - qoff1)/(-m); % x = (ymax-qoff)/m
    ysup1 = min(D1(:,1))-1;
    xsup1 = (ysup1 - qoff1)/(-m);
    yinf2 = max(D2(:,1))+1;
    xinf2 = (yinf2 - qoff2)/m;
    ysup2 = min(D2(:,1))-1;
    xsup2 = (ysup2 - qoff2)/m;

    [xoffset1,yoffset1] = polyxpoly(D1(:,2),D1(:,1),[xinf1 xsup1],[yinf1 ysup1],'unique'); 
    [xoffset2,yoffset2] = polyxpoly(D2(:,2),D2(:,1),[xsup2 xinf2],[ysup2 yinf2],'unique'); 

    % Insert here volume dimension control (dim < 2; dim == 2; dim > 2)
    % Plot integrated
    
    if isempty(xoffset1) || isempty(yoffset1)
        D3 = [];
    else
        D3 = [max(yoffset1) xoffset1(yoffset1==max(yoffset1)) 0 0; D1(D1(:,1)<(-m*(D1(:,2))+qoff1),:); min(yoffset1) xoffset1(yoffset1==min(yoffset1)) 0 0; max(yoffset1) xoffset1(yoffset1==max(yoffset1)) 0 0];
    end

    if isempty(xoffset2) || isempty(yoffset2)
        D4 = [];
    else
        D4 = [min(yoffset2) xoffset2(yoffset2==min(yoffset2)) 0 0; D2(D2(:,1)<(m*(D2(:,2))+qoff2),:); max(yoffset2) xoffset2(yoffset2==max(yoffset2)) 0 0; min(yoffset2) xoffset2(yoffset2==min(yoffset2)) 0 0];  
    end
    
end
% plot(O(:,2),-O(:,1),'Color',grey)
% plot(nonfeasdom([1,end],2),-nonfeasdom([1,end],1),'r.','MarkerSize',20)

end
