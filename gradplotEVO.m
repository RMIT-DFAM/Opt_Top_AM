function gradplotEVO(NodeStruc,choice)

% global m red green blue orange

% NodeStruc = [y x gradient feasible/not R G B];
red    = [185 5 0];
orange = [225 205 20];
blue   = [8 18 149];
% purple = [229 23 229];
% green  = [51 204 0];
% Main colors
% grey    = [0.7 0.7 0.7];
% red     = [230 0 10];
% blue    = [8 18 149];
% orange  = [255 230 0];
% green   = [51 204 0];

ymax   = max(NodeStruc{1,1}(:,1));
delta  = 0.01;
maxang = 45;  % should be global

switch choice
    
case 1 %%% Gradient variation
    for ii = 1:size(NodeStruc,1)
        for jj = 1:size(NodeStruc{ii,1},1)
            if NodeStruc{ii,1}(jj,3)>=0 && NodeStruc{ii,1}(jj,3)<=90
                % Red:   orange(1)-->red(1)
                % green: orange(2)-->red(2)
                % Blue:  orange(3)-->red(3)
                dr = abs(orange(1)-red(1))/90;
                dg = abs(orange(2)-red(2))/90;
                db = abs(orange(3)-red(3))/90;
                NodeStruc{ii,1}(jj,5) = orange(1)-NodeStruc{ii,1}(jj,3)*dr; % red
                NodeStruc{ii,1}(jj,6) = orange(2)-NodeStruc{ii,1}(jj,3)*dg; % Green
                NodeStruc{ii,1}(jj,7) = orange(3)-NodeStruc{ii,1}(jj,3)*db; % Blue
            elseif NodeStruc{ii,1}(jj,3)>90 && NodeStruc{ii,1}(jj,3)<=180
                % Red:   red(1)-->blue(1)
                % green: red(2)-->blue(2)
                % Blue:  red(3)-->blue(3)   
                dr = abs(red(1)-blue(1))/90;
                dg = abs(red(2)-blue(2))/90;
                db = abs(red(3)-blue(3))/90;
                NodeStruc{ii,1}(jj,5) = red(1)-(NodeStruc{ii,1}(jj,3)-90)*dr; % Red
                NodeStruc{ii,1}(jj,6) = red(2)+(NodeStruc{ii,1}(jj,3)-90)*dg; % green
                NodeStruc{ii,1}(jj,7) = red(3)+(NodeStruc{ii,1}(jj,3)-90)*db; % Blue
            elseif NodeStruc{ii,1}(jj,3)>180 && NodeStruc{ii,1}(jj,3)<=270
                % Red:   blue(1)-->red(1)
                % green: blue(2)-->red(2)
                % Blue:  blue(3)-->red(3)   
                dr = abs(red(1)-blue(1))/90;
                dg = abs(red(2)-blue(2))/90;
                db = abs(red(3)-blue(3))/90;
                NodeStruc{ii,1}(jj,5) = blue(1)+(NodeStruc{ii,1}(jj,3)-180)*dr; % Red
                NodeStruc{ii,1}(jj,6) = blue(2)-(NodeStruc{ii,1}(jj,3)-180)*dg; % Green
                NodeStruc{ii,1}(jj,7) = blue(3)-(NodeStruc{ii,1}(jj,3)-180)*db; % Blue            

            elseif NodeStruc{ii,1}(jj,3)>270 && NodeStruc{ii,1}(jj,3)<=360
                % Red:   red(1)-->orange(1)
                % Green: red(2)-->orange(2)
                % Blue:  red(3)-->orange(3)
                dr = abs(orange(1)-red(1))/90;
                dg = abs(orange(2)-red(2))/90;
                db = abs(orange(3)-red(3))/90;
                NodeStruc{ii,1}(jj,5) = red(1)+(NodeStruc{ii,1}(jj,3)-270)*dr; % Red
                NodeStruc{ii,1}(jj,6) = red(2)+(NodeStruc{ii,1}(jj,3)-270)*dg; % Green
                NodeStruc{ii,1}(jj,7) = red(3)+(NodeStruc{ii,1}(jj,3)-270)*db; % Blue
            end
        end
    end

    %%% Gradient PLOT
    clear ii jj
    figure
    hold on
    for ii = 1:size(NodeStruc,1)
        for jj = 1:size(NodeStruc{ii,1},1)-1
            plot(NodeStruc{ii,1}(jj:jj+1,2),-NodeStruc{ii,1}(jj:jj+1,1),'LineWidth',3,'Color',[NodeStruc{ii,1}(jj,5) NodeStruc{ii,1}(jj,6) NodeStruc{ii,1}(jj,7)]./255);
        end
    end
    axis image
    colorbar

case 2 %%% Feasible and infeasible areas

    for ii = 1:size(NodeStruc,1)
        for jj = 1:size(NodeStruc{ii,1})
            if ii==1 % EXTERNAL BOND
                if NodeStruc{ii,1}(jj,3)>180-maxang && NodeStruc{ii,1}(jj,3)<180+maxang && NodeStruc{ii,1}(jj,1)<ymax-delta % (135)-(225) 
                    % All RED
                    NodeStruc{ii,1}(jj,5) = 230;  % Red
                    NodeStruc{ii,1}(jj,6) = 0;  % Green
                    NodeStruc{ii,1}(jj,7) = 10;  % Blue
                elseif NodeStruc{ii,1}(jj,3)>180-maxang-5 && NodeStruc{ii,1}(jj,3)<=180-maxang % (130)-135
                    % Red:   255-->230
                    % Green: 230-->0
                    % Blue:  0-->10
                    dr = (255-230)/5;
                    dg = 230/5;
                    db = 10/5;
                    NodeStruc{ii,1}(jj,5) = 255-(NodeStruc{ii,1}(jj,3)-(180-maxang-5))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 230-(NodeStruc{ii,1}(jj,3)-(180-maxang-5))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 0+(NodeStruc{ii,1}(jj,3)-(180-maxang-5))*db;   % Blue
                elseif NodeStruc{ii,1}(jj,3)>=180+maxang && NodeStruc{ii,1}(jj,3)<180+maxang+5 % 225-(230)
                    % Red: 230-->255
                    % Green: 0-->230
                    % Blue: 10-->0
                    dr = (255-230)/5;
                    dg = 230/5;
                    db = 10/5;
                    NodeStruc{ii,1}(jj,5) = 230+(NodeStruc{ii,1}(jj,3)-(180+maxang))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 0+(NodeStruc{ii,1}(jj,3)-(180+maxang))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 10-(NodeStruc{ii,1}(jj,3)-(180+maxang))*db;   % Blue
                elseif NodeStruc{ii,1}(jj,3)>180-maxang-10 && NodeStruc{ii,1}(jj,3)<=180-maxang-5 % (125)-130
                    % Red:   51-->255
                    % Green: 204-->230
                    % Blue:    0-->0
                    dr = (255-51)/5;
                    dg = (230-204)/5;
                    NodeStruc{ii,1}(jj,5) = 51+(NodeStruc{ii,1}(jj,3)-(180-maxang-10))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 204+(NodeStruc{ii,1}(jj,3)-(180-maxang-10))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 0;   % Blue                    
                elseif NodeStruc{ii,1}(jj,3)>=180+maxang+5 && NodeStruc{ii,1}(jj,3)<180+maxang+10 % 230-(235)                    
                    % Red:   255-->51
                    % Green: 230-->204
                    % Blue:    0-->0
                    dr = (255-51)/5;
                    dg = (230-204)/5;
                    NodeStruc{ii,1}(jj,5) = 255-(NodeStruc{ii,1}(jj,3)-(180+maxang+5))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 230+(NodeStruc{ii,1}(jj,3)-(180+maxang+5))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 0;   % Blue
                else
                    % gradient ok --> Green
                    NodeStruc{ii,1}(jj,5) = 51; % Red
                    NodeStruc{ii,1}(jj,6) = 204;  % Green
                    NodeStruc{ii,1}(jj,7) = 0;   % Blue
                end
                  
            else % INTERNAL BONDS
                if NodeStruc{ii,1}(jj,3)>=0 && NodeStruc{ii,1}(jj,3)<maxang % 0-(45) 
                    % All RED
                    NodeStruc{ii,1}(jj,5) = 230;  % Red
                    NodeStruc{ii,1}(jj,6) = 0;  % Green
                    NodeStruc{ii,1}(jj,7) = 10;  % Blue
                elseif NodeStruc{ii,1}(jj,3)>=360-maxang && NodeStruc{ii,1}(jj,3)<=360 % (315)-360 
                    % All RED
                    NodeStruc{ii,1}(jj,5) = 230;  % Red
                    NodeStruc{ii,1}(jj,6) = 0;  % Green
                    NodeStruc{ii,1}(jj,7) = 10;  % Blue
                elseif NodeStruc{ii,1}(jj,3)>=maxang && NodeStruc{ii,1}(jj,3)<maxang+5 % 45-(50)  
                    % Red: 230-->255
                    % Green: 0-->230
                    % Blue: 10-->0
                    dr = (255-230)/5;
                    dg = 230/5;
                    db = 10/5;
                    NodeStruc{ii,1}(jj,5) = 230+(NodeStruc{ii,1}(jj,3)-maxang)*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 0+(NodeStruc{ii,1}(jj,3)-maxang)*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 10-(NodeStruc{ii,1}(jj,3)-maxang)*db;   % Blue
                elseif NodeStruc{ii,1}(jj,3)>=360-maxang-5 && NodeStruc{ii,1}(jj,3)<360-maxang % 310-(315) 
                    % Red:   255-->230
                    % Green: 230-->0
                    % Blue:  0-->10
                    dr = (255-230)/5;
                    dg = 230/5;
                    db = 10/5;
                    NodeStruc{ii,1}(jj,5) = 255-(NodeStruc{ii,1}(jj,3)-(360-maxang-5))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 230-(NodeStruc{ii,1}(jj,3)-(360-maxang-5))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 0+(NodeStruc{ii,1}(jj,3)-(360-maxang-5))*db;   % Blue
                elseif NodeStruc{ii,1}(jj,3)>=maxang+5 && NodeStruc{ii,1}(jj,3)<maxang+10 % 50-(55)  
                    % Red:   255-->51
                    % Green: 230-->204
                    % Blue:    0-->0
                    dr = (255-51)/5;
                    dg = (230-204)/5;
                    NodeStruc{ii,1}(jj,5) = 255-(NodeStruc{ii,1}(jj,3)-(maxang+5))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 230+(NodeStruc{ii,1}(jj,3)-(maxang+5))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 0;   % Blue
                elseif NodeStruc{ii,1}(jj,3)>=360-maxang-10 && NodeStruc{ii,1}(jj,3)<360-maxang-5 % 305-(310) 
                    % Red:   51-->255
                    % Green: 204-->230
                    % Blue:    0-->0
                    dr = (255-51)/5;
                    dg = (230-204)/5;
                    NodeStruc{ii,1}(jj,5) = 51+(NodeStruc{ii,1}(jj,3)-(360-maxang-10))*dr; % Red
                    NodeStruc{ii,1}(jj,6) = 204-(NodeStruc{ii,1}(jj,3)-(360-maxang-10))*dg;  % Green
                    NodeStruc{ii,1}(jj,7) = 0;   % Blue
                else
                    % gradient ok --> Green
                    NodeStruc{ii,1}(jj,5) = 51; % Red
                    NodeStruc{ii,1}(jj,6) = 204;  % Green
                    NodeStruc{ii,1}(jj,7) = 0;   % Blue        
                end
            end
        end
    end

    %%% Feasible/Infeasible PLOT

    clear ii jj
    load('MyColormaps','gradientmap')

    figure;
    set(gcf,'Colormap',gradientmap)
    hold on
    for ii = 1:size(NodeStruc,1)
        for jj = 1:size(NodeStruc{ii,1},1)-1
            plot(NodeStruc{ii,1}(jj:jj+1,2),-NodeStruc{ii,1}(jj:jj+1,1),'LineWidth',3,'Color',[NodeStruc{ii,1}(jj,5) NodeStruc{ii,1}(jj,6) NodeStruc{ii,1}(jj,7)]./255);
        end
    end
    axis image
    labels = {'Feasible','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','Almost infeasible','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','Infeasible'};
    lcolorbar(labels);
end
end