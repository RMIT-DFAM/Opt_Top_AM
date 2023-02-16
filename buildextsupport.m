%% buildextsupport.m
% This function builds support structures for external boundaries
% extbound = matrix [y x grad bin]
% unfeasibles = matrix [y x]


function [D,newshape] = buildextsupport(extbound)



global tmin m outflag alpha indexx epsilon

rows    = size(extbound,1); 
xmin    = min(extbound(:,2));
xmax    = max(extbound(:,2));
ymin    = min(extbound(:,1));
ymax    = max(extbound(:,1));
c1      = 0;
c2      = 0;
jj      = 1;

% if alpha==60  % DEBUG
%     keyboard
% end

%%% Reorder extbound
pointer = find(extbound(:,1)==min(extbound(extbound(:,2)==xmin,1)),1,'first');
if pointer==1
    extbound = extbound;
else
    extbound = [extbound(pointer:end,:) ; extbound(2:pointer,:)];
end

extbound = [zeros(1,4) ; extbound ; zeros(1,4)];
unfeasibles = [];
newshape = [];

%%% Identification and storage of unfeasible areas %%%

for ii = 2:size(extbound,1)-2
   if extbound(ii,4)==1 && extbound(ii-1,4)==0 && extbound(ii+1,4)==1 && extbound(ii+2,4)==1
        c1 = ii;
   end
   if extbound(ii,4)==1 && extbound(ii-1,4)==1 && extbound(ii-2,4)==1 && extbound(ii+1,4)==0
        c2 = ii;
        unfeasibles{jj,1} = extbound(c1:c2,1:2);
        jj = jj+1;
   end
end
extbound = extbound(2:end-1,:);
doubleextbound = [extbound;extbound(2:end,:)];

if isempty(unfeasibles) % if empty --> exit and quit loop!
   D = [];
   newshape = extbound;
   outflag = 1;
else
    highest = zeros(1,size(unfeasibles,1));  % Start logic with the highest infeas region
    for ii = 1:size(unfeasibles,1)
        highest(ii) = min(unfeasibles{ii,1}(:,1));
    end
    UNFvector = unfeasibles{highest==min(highest),1};    
%     UNFvector = unfeasibles{1,1};

    %%% Structure building
    
    % RETTA DA DX VERSO SX
    ydx = UNFvector(1,1);
    xdx = UNFvector(1,2);
    qdx = ydx+m*xdx;
    
    % trova intersezione a sx
    indL  = find(ismember(doubleextbound(:,1:2),[ydx xdx],'rows'),1,'first')+1;
%     indLL = indL-1;
%     exitflag = 0;
    foundL   = 0;
    intersectionsL = [];
    
    % new loop for intersection seeking on left dom
    
    for ii = indL:indL+rows-3
       xLpre  = doubleextbound(ii-1,2);
       yLpre  = doubleextbound(ii-1,1);
       xLcurr = doubleextbound(ii,2);
       yLcurr = doubleextbound(ii,1);
       if xLcurr>=(-yLcurr+qdx)/m && xLpre<(-yLpre+qdx)/m && xLcurr<UNFvector(1,2) && yLcurr>UNFvector(1,1)
          intersectionsL = [intersectionsL; yLcurr xLcurr sqrt((UNFvector(1,2)-xLcurr)^2+(UNFvector(1,1)-yLcurr)^2)]; 
       end
    end
    if isempty(intersectionsL)==0 % intersection found!
        foundL = 1;
        xL    = intersectionsL(intersectionsL(:,3)==min(intersectionsL(:,3)),2);
        yL    = intersectionsL(intersectionsL(:,3)==min(intersectionsL(:,3)),1);    
        indLL = find(ismember(doubleextbound(2:end-1,1:2),[yL xL],'rows'),1,'first')+1;
    end    
    
    %%% old LOOP to be substituted
%     while exitflag == 0 
%        indLL = indLL+1;
%        xL = doubleextbound(indLL,2);
%        yL = doubleextbound(indLL,1);
%        if indLL >= rows  % Cerco intersezione fino al punto iniziale (SHOULD BE MODIFIED: almost 1 round or check ccw if in left side )
%            exitflag = 1;
%        elseif xL>(-yL+qdx)/m && xL<UNFvector(1,2) && yL>UNFvector(1,1) % trovato il punto giusto
%            exitflag = 1;
%            foundL   = 1;
%        end
%     end
    %%%
    
    if foundL % trovato intersezione con componente
        if indLL<indL-1 % intersection is previous --> Result domain is open!
            AreaL = Inf;
        else            % closed domain --> compute Area
            mpL  = (yL-doubleextbound(indLL-1,1))/(xL-doubleextbound(indLL-1,2)); % segment coefficient
            if abs(mpL) ~= Inf
                qpL = yL-mpL*xL; 
                xyL = [1 m;1 -mpL]\[qdx; qpL]; % A\b
                ypL = xyL(1); % domain intersection point coordinates
                xpL = xyL(2);
            else
                xpL = xL;
                ypL = -m*xpL+qdx;
            end
            dL    = [ypL xpL 0 0; doubleextbound(indLL-1:-1:indL-1,1:4); ypL xpL 0 0];        
            AreaL = polyarea(dL(:,2),dL(:,1));
        end 
    
    else   % non trovato intersezione --> cerco intersezione con asse x=0
        ygroundL = ymax;
        xgroundL = (-ygroundL+qdx)/m;
        if xgroundL<xmin || xgroundL>xmax
            AreaL = Inf;
        elseif isempty(doubleextbound(find(doubleextbound(2:rows,1)<ymax+epsilon&doubleextbound(2:rows,1)>ymax-epsilon&doubleextbound(2:rows,2)<xgroundL)+1,:))==0  % has following points with y=ymax
            dL    = [ygroundL xgroundL 0 0; doubleextbound(find(doubleextbound(2:rows,1)<ymax+epsilon&doubleextbound(2:rows,1)>ymax-epsilon&doubleextbound(2:rows,2)<xgroundL,1,'first')+1:-1:indL-1,1:4) ;  ygroundL xgroundL 0 0];
            AreaL = polyarea(dL(:,2),dL(:,1));
%         elseif max(doubleextbound(doubleextbound(1:rows,1)==ymax,2))<xgroundL
%             dL    = [ygroundL xgroundL 0 0; doubleextbound(find(doubleextbound(2:end,1)==ymax,1,'first')+1:-1:indL-1,1:4) ;  ygroundL xgroundL 0 0];
%             AreaL = polyarea(dL(:,2),dL(:,1));
        else
            AreaL = Inf;
        end
    end
    
    % RETTA DA SX VERSO DX
    ysx = UNFvector(end,1);
    xsx = UNFvector(end,2);
    qsx = ysx-m*xsx;
    
    % trova intersezione a dx
    indR     = find(ismember(doubleextbound(:,1:2),[ysx xsx],'rows'),1,'last')-1; % MODIFICATO (1:end-1) con (:)  
%     indRR    = indR+1;
%     exitflag = 0;
    foundR   = 0;
    intersectionsR = [];
    
    % new loop for intersection seeking in the right dom    
    for ii = indR:-1:indR-rows+3
       xRpre  = doubleextbound(ii+1,2);
       yRpre  = doubleextbound(ii+1,1);
       xRcurr = doubleextbound(ii,2);
       yRcurr = doubleextbound(ii,1);
       if xRcurr<=(yRcurr-qsx)/m && xRpre>(yRpre-qsx)/m && xRcurr>UNFvector(end,2) && yRcurr>UNFvector(end,1)
          intersectionsR = [intersectionsR; yRcurr xRcurr sqrt((UNFvector(end,2)-xRcurr)^2+(UNFvector(end,1)-yRcurr)^2)]; 
       end
    end
    if isempty(intersectionsR)==0
        foundR = 1;
        xR    = intersectionsR(intersectionsR(:,3)==min(intersectionsR(:,3)),2);
        yR    = intersectionsR(intersectionsR(:,3)==min(intersectionsR(:,3)),1);    
        indRR = find(ismember(doubleextbound(2:end-1,1:2),[yR xR],'rows'),1,'last')+1;
    end   

%     %%% old loop for rigth dom
%     while exitflag == 0 
%        indRR = indRR-1;
%        xR = doubleextbound(indRR,2);
%        yR = doubleextbound(indRR,1);
%        if indRR <= rows  % Cerco intersezione fino al punto iniziale
%            exitflag = 1;
%        elseif xR<(yR-qsx)/m && xR>UNFvector(end,2) && yR>UNFvector(end,1) % trovato il punto giusto
%            exitflag = 1;
%            foundR   = 1;
%        end
%     end
%     %%%
    
    if foundR % trovato intersezione con componente
        if indRR>indR+1 % intersection is previous --> Result domain is open!
            AreaR = Inf;
        else  
            mpR  = (yR-doubleextbound(indRR+1,1))/(xR-doubleextbound(indRR+1,2)); % segment coefficient
            if abs(mpR) ~= Inf
                qpR = yR-mpR*xR; 
                xyR = [1 -m;1 -mpR]\[qsx; qpR]; % A\b
                ypR = xyR(1); % domain intersection point coordinates
                xpR = xyR(2);
            else
                xpR = xR;
                ypR = m*xpR+qsx;
            end
            dR    = [doubleextbound(indR+1:-1:indRR+1,1:4); ypR xpR 0 0; doubleextbound(indR+1,1:4)];        
            AreaR = polyarea(dR(:,2),dR(:,1));
        end
    else    % non trovato intersezione --> cerco intersezione con asse x=0
        ygroundR = ymax;
        xgroundR = (ygroundR-qsx)/m;
        if xgroundR<xmin || xgroundR>xmax
            AreaR = Inf;
        elseif isempty(doubleextbound(find(doubleextbound(2:rows,1)<ymax+epsilon&doubleextbound(2:rows,1)>ymax-epsilon&doubleextbound(2:rows,2)>xgroundR)+1,:))==0  % has previous points with y=ymax
            dR    = [ygroundR xgroundR 0 0; doubleextbound(indR+1-rows:-1:find(doubleextbound(2:rows,1)<ymax+epsilon&doubleextbound(2:rows,1)>ymax-epsilon&doubleextbound(2:rows,2)>xgroundR,1,'last')+1,:) ;  ygroundR xgroundR 0 0];
            AreaR = polyarea(dR(:,2),dR(:,1));
            
%         elseif min(extbound(extbound(:,1)==ymax,2))>xgroundR
%             dR    = [doubleextbound(indR+1:-1:find(doubleextbound(:,1)==ymax,1,'last'),1:4) ; ygroundR xgroundR 0 0 ; doubleextbound(indR+1,1:4)];
%             AreaR = polyarea(dR(:,2),dR(:,1));
        else
            AreaR = Inf;
        end
    end
    
    %%% Collect minimum area external domains
    %%% Salvo le aree aggiuntive del contorno più esterno; il materiale
    %%% aggiunto all'esterno del dominio originario sarà dato dalla
    %%% differenza tra queste aree e i domini vuoti calcolati in seguito
    
    if AreaL < AreaR         
        qoffL = qdx - tmin/cos(atan(-m)); % y+mx
        yinfL = max(dL(:,1))+1;         % y = mx+q
        xinfL = (yinfL - qoffL)/(-m); % x = (ymax-qoff)/m
        ysupL = min(dL(:,1))-1;
        xsupL = (ysupL - qoffL)/(-m);
        
        [xoffsetL,yoffsetL] = polyxpoly(dL(:,2),dL(:,1),[xinfL xsupL],[yinfL ysupL],'unique'); 

        if isempty(xoffsetL) || isempty(yoffsetL)
            D = [];
        else 
            D = [max(yoffsetL) xoffsetL(yoffsetL==max(yoffsetL)) 0 0; dL(dL(:,2)<(-dL(:,1)+qoffL)/m,:); min(yoffsetL) xoffsetL(yoffsetL==min(yoffsetL)) 0 0; max(yoffsetL) xoffsetL(yoffsetL==max(yoffsetL)) 0 0]; 
        end
        
        % Costruzione contorno esterno aggiornato Left
        if foundL
           newshape = [extbound(1:indL-1,:) ; ypL xpL 0 0 ; extbound(indLL:end,:)]; 
        else
           newshape = [extbound(1:indL-1,:) ; ygroundL xgroundL 0 0; extbound(find(extbound(2:end,1)<ymax+epsilon&extbound(2:end,1)>ymax-epsilon&extbound(2:end,2)<xgroundL,1,'first')+1:end,:)]; % Modded
        end
         
    else % Meglio areaR
        
%         if alpha==135
%             keyboard  %%% DEBUG FROM HERE UPON, PROBLEM WITH dR=[]
%         end
        
        qoffR = qsx - tmin/cos(atan(m));
        yinfR = max(dR(:,1))+1;
        xinfR = (yinfR - qoffR)/m;
        ysupR = min(dR(:,1))-1;
        xsupR = (ysupR - qoffR)/m; 

        [xoffsetR,yoffsetR] = polyxpoly(dR(:,2),dR(:,1),[xsupR xinfR],[ysupR yinfR],'unique'); 

        if isempty(xoffsetR) || isempty(yoffsetR)
            D = [];
        else 
            D = [min(yoffsetR) xoffsetR(yoffsetR==min(yoffsetR)) 0 0; dR(dR(:,2)>(dR(:,1)-qoffR)/m,:); max(yoffsetR) xoffsetR(yoffsetR==max(yoffsetR)) 0 0; min(yoffsetR) xoffsetR(yoffsetR==min(yoffsetR)) 0 0];  
        end
        
        % Costruzione contorno esterno aggiornato Right
        if foundR
           newshape = [doubleextbound(rows:indRR,:) ; ypR xpR 0 0 ; doubleextbound(indR+1:end,:)]; 
        else
           newshape = [doubleextbound(rows:find(doubleextbound(2:rows,1)<ymax+epsilon&doubleextbound(2:rows,1)>ymax-epsilon&doubleextbound(2:rows,2)>xgroundR,1,'last')+rows,:) ; ygroundR xgroundR 0 0; doubleextbound(indR+1:end,:)]; % Modded 
        end
    end
end