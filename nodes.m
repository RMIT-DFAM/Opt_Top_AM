%% nodes.m
% THIS FUNCTION IS USED TO CREATE BOUNDARY NODE COORDINATE FROM A BINARY MATRIX
% INPUT STARTING FROM THE SOUTH WEST CORNER OF 2D PLANAR GRID
function [n] = nodes(b,ind)

global nelx nely binaryElem
% Input = binary matrix binaryElem
% Output = n vector with boundary nodes coordinates

x = b(:,2);
y = b(:,1);
ref = 1;
n = [];
jj = 1;
cazzodiflag = [0 0];

if ind == 1  % PERIMETER LOGIC
    
   ref = 0;                    % reference elem value (0 for perimeter, 1 for internal voids)
   n = [y(1) x(1)-1; y(1)-1 x(1)-1]; % [y x-1] Two starting nodes (first SW point and NW point)
   
    for ii = 1:size(b,1)
       
        if ii<size(b,1) && x(ii+1)>x(ii) && y(ii+1)<y(ii) % MOVING ABOVE-RIGHT
            if ii>1 && x(ii+1)==x(ii-1) && y(ii+1)==y(ii-1) % Corner element
                n = [n; y(ii) x(ii); y(ii) x(ii)-1; y(ii)-1 x(ii)-1; y(ii)-1 x(ii)];    
            else                                    % Simply diagonal     
                if x(ii)-1<=0 || binaryElem(y(ii),x(ii)-1)==ref        % control elem on the left
                    if sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)>0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)==0 % "left-side" node exist, no right-side node
                        n = [n ; y(ii)-1 x(ii)-1];
                    elseif sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)==0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)==0
                       cazzodiflag(1) = 1;  
                    end
                end
                if y(ii)-1<=0 || binaryElem(y(ii)-1,x(ii))==ref        % control elem above
                    if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)>0  % "left-side" node exist
                        n = [n ; y(ii)-1 x(ii)];
                    else
                        cazzodiflag(2) = 1;
                    end
                end
                if sum(cazzodiflag) == 2
                     n = [n ; y(ii) x(ii)-1 ; y(ii)-1 x(ii)-1 ; y(ii)-1 x(ii)];     
                end
                cazzodiflag = [0 0];
            end        
            
        elseif ii<size(b,1) && x(ii+1)>x(ii) && y(ii+1)>y(ii)  % MOVING BELOW-RIGHT
            if ii>1 && x(ii+1)==x(ii-1) && y(ii+1)==y(ii-1) % Corner element            
                n = [n; y(ii) x(ii)-1; y(ii)-1 x(ii)-1; y(ii)-1 x(ii); y(ii) x(ii)];
            else                                    % Simply diagonal
                if y(ii)-1<=0 || binaryElem(y(ii)-1,x(ii))==ref       % control elem above
                    if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)>0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))==0 % "left-side" node exist, no right-side node
                        n = [n ; y(ii)-1 x(ii)];
                    elseif sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)==0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))==0
                        cazzodiflag(1) = 1;
                    end
                end
                if x(ii)+1>=nelx+1 || binaryElem(y(ii),x(ii)+1)==ref  % control elem on the right
                    if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))>0  % "left-side" node exist
                        n = [n ; y(ii) x(ii)];
                    else
                        cazzodiflag(2) = 1;
                    end
                end
                if sum(cazzodiflag) == 2
                     n = [n ; y(ii)-1 x(ii)-1 ; y(ii)-1 x(ii) ; y(ii) x(ii)];     
                end
                cazzodiflag = [0 0];
            end
            
            
        elseif ii<size(b,1) && x(ii+1)<x(ii) && y(ii+1)>y(ii)  % MOVING BELOW-LEFT
            if ii>1 && x(ii+1)==x(ii-1) && y(ii+1)==y(ii-1) % Corner element            
                n = [n; y(ii)-1 x(ii)-1; y(ii)-1 x(ii); y(ii) x(ii); y(ii) x(ii)-1];
            else                                    % Simply diagonal
                if x(ii)+1>=nelx+1 || binaryElem(y(ii),x(ii)+1)==ref   % control elem on the right
                    if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))>0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii))==0 % "left-side" node exist, no right-side node
                        n = [n ; y(ii) x(ii)];
                    elseif sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))==0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii))==0
                        cazzodiflag(1) = 1;
                    end
                end
                if y(ii)+1>=nely+1 || binaryElem(y(ii)+1,x(ii))==ref  % control elem below
                    if sum(n(:,1)==y(ii) & n(:,2)==x(ii))>0  % "left-side" node exist
                        n = [n ; y(ii) x(ii)-1];
                    else
                        cazzodiflag(2) = 1;
                    end
                end
                if sum(cazzodiflag) == 2
                     n = [n ; y(ii)-1 x(ii) ; y(ii) x(ii) ; y(ii) x(ii)-1];     
                end
                cazzodiflag = [0 0];
            end
         
            
        elseif ii<size(b,1) && x(ii+1)<x(ii) && y(ii+1)<y(ii)  % MOVING ABOVE-LEFT
            if ii>1 && x(ii+1)==x(ii-1) && y(ii+1)==y(ii-1)  % Corner element            
                n = [n; y(ii)-1 x(ii); y(ii) x(ii); y(ii) x(ii)-1; y(ii)-1 x(ii)-1];
            else                                     % Simple diagonal
                if y(ii)+1>=nely+1 || binaryElem(y(ii)+1,x(ii))==ref  % control elem below
                    if sum(n(:,1)==y(ii) & n(:,2)==x(ii))>0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)==0  % "left-side" node exist, no right-side node
                        n = [n ; y(ii) x(ii)-1];
                    elseif sum(n(:,1)==y(ii) & n(:,2)==x(ii))==0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)==0
                        cazzodiflag(1) = 1;
                    end
                end    
                if x(ii)-1<=0 || binaryElem(y(ii),x(ii)-1)==ref        % control elem on the left
                    if sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)>0 % "left-side" node exist, no right-side node
                        n = [n ; y(ii)-1 x(ii)-1];
                    else
                        cazzodiflag(2) = 1;
                    end
                end
                if sum(cazzodiflag) == 2
                     n = [n ; y(ii) x(ii) ; y(ii) x(ii)-1 ; y(ii)-1 x(ii)-1];     
                end
                cazzodiflag = [0 0];
            end
            
            
        else  % NORMAL ELEMENT
            
            % Double loop
            while jj<=2
                if x(ii)-1<=0 || binaryElem(y(ii),x(ii)-1)==ref        % control elem on the left
                    if sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)>0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)==0 % "left-side" node exist
                        n = [n ; y(ii)-1 x(ii)-1];
                    end
                end
                if y(ii)-1<=0 || binaryElem(y(ii)-1,x(ii))==ref        % control elem above
                    if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)>0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))==0 % "left-side" node exist
                        n = [n ; y(ii)-1 x(ii)];
                    end
                end
                if x(ii)+1>=nelx+1 || binaryElem(y(ii),x(ii)+1)==ref        % control elem on the right
                    if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))>0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii))==0 % "left-side" node exist
                        n = [n ; y(ii) x(ii)];
                    end
                end
                if y(ii)+1>=nely+1 || binaryElem(y(ii)+1,x(ii))==ref        % control elem below
                    if sum(n(:,1)==y(ii) & n(:,2)==x(ii))>0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)==0 % "left-side" node exist
                        n = [n ; y(ii) x(ii)-1];
                    end
                end
                jj=jj+1;
            end
            jj=1;
        end
    end
        
else  % VOIDS LOGIC
    
    n = [y(1) x(1)-1]; % [y x-1] Starting node
    
    for ii = 1:size(b,1)
        % Double loop
        while jj<=2
            if binaryElem(y(ii),x(ii)-1)==ref        % control elem on the left
                if sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)>0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)==0 % "left-side" node exist
                    n = [n ; y(ii)-1 x(ii)-1];
                end
            end
            if binaryElem(y(ii)-1,x(ii))==ref        % control elem above
                if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii)-1)>0 && sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))==0 % "left-side" node exist
                    n = [n ; y(ii)-1 x(ii)];
                end
            end
            if binaryElem(y(ii),x(ii)+1)==ref        % control elem on the right
                if sum(n(:,1)==y(ii)-1 & n(:,2)==x(ii))>0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii))==0 % "left-side" node exist
                    n = [n ; y(ii) x(ii)];
                end
            end
            if binaryElem(y(ii)+1,x(ii))==ref        % control elem below
                if sum(n(:,1)==y(ii) & n(:,2)==x(ii))>0 && sum(n(:,1)==y(ii) & n(:,2)==x(ii)-1)==0 % "left-side" node exist
                    n = [n ; y(ii) x(ii)-1];
                end
            end
            jj=jj+1;
        end
        jj=1;
    end
end

n = [n; n(1,1) n(1,2)];
