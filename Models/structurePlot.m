function structurePlot(X,omega,constraints,camera,clearFig,cylinders,cableLabel)
% This function will plot the tensegrity structure. The inputs for this
% function are the nodal positions, Cable Connectivity matrix, Rod
% Connectivity matrix, and the constraint matrix. The nodes that have
% constraints will appear as a star with the color that corresponds to the
% constraint applied to it. The nodes with no constraints will appear as
% solid blue circle. The rods will appear as a black line and the cable
% will appear as a solid orange cable. This code will also check the
% correctness of the cable connectivity, rod connectivity, and constraint
% matrices.

%%
%future updates:
%add ability to label nodes, cables, rods
%%

Verts = X.p;
C = omega.C;
R = omega.R;

%figure(1)
if(clearFig)
    cla
end
grid on; hold on; xlabel('X'); ylabel('Y'); zlabel('Z');

color1 = omega.cables.paired(1:3:end,:);
color2 = omega.cables.paired(2:3:end,:);
color3 = omega.cables.paired(3:3:end,:);
color1 = reshape(color1,1,numel(color1));
color2 = reshape(color2,1,numel(color2));
color3 = reshape(color3,1,numel(color3));
for n=1:size(C,1) %For each row in C find the connected vertices and plot a cable
    index = find(C(n,:)~=0); %index of non-zero values in the nth row of C (MUST HAVE 2 ELEMENTS)
    % index correlates to the row in Verts to make a connection
    x = [Verts(index(1)*3-2),Verts(index(2)*3-2)];
    y = [Verts(index(1)*3-1),Verts(index(2)*3-1)];
    z = [Verts(index(1)*3),Verts(index(2)*3)];
    if(~any(omega.cables.passive==n))
        rope = line(x,y,z); rope.LineWidth = 2; rope.Color = 'k';
    elseif(any(color1==n))
        rope = line(x,y,z); rope.LineWidth = 2; rope.Color = 'm';
    elseif(any(color2==n))
        rope = line(x,y,z); rope.LineWidth = 2; rope.Color = 'c';
    elseif(any(color3==n))
        rope = line(x,y,z); rope.LineWidth = 2; rope.Color = 'k';
    else
        rope = line(x,y,z); rope.LineWidth = 1; rope.Color = 'k';
    end
    if(cableLabel)
        mid = (Verts(index(1)*3-2:index(1)*3,:)+...
            Verts(index(2)*3-2:index(2)*3,:))/2;
        text(mid(1),mid(2),mid(3),num2str(n),'Color','blue','FontSize',14);
    end
end

N = 20;
for n=1:size(R,1) %For each row in R find the connected vertices and plot a rod
    index = find(R(n,:)~=0); %index of non-zero values in the nth row of R (MUST HAVE 2 ELEMENTS)
    % index correlates to the row in Verts to make a connection
    x = [Verts(index(1)*3-2),Verts(index(2)*3-2)];
    y = [Verts(index(1)*3-1),Verts(index(2)*3-1)];
    z = [Verts(index(1)*3),Verts(index(2)*3)];
    if(cylinders)
        xyz = [x;y;z];
        [x,y,z]=cylinder2P(0.0125,N,xyz(:,1)',xyz(:,2)');
        if(mod(n,3)==0)
            surf(x,y,z,'FaceColor',[1 0 0],'EdgeAlpha',0.5);
        elseif(mod(n,3)==1)
            surf(x,y,z,'FaceColor',[0 1 0],'EdgeAlpha',0.5);
        else
            surf(x,y,z,'FaceColor',[0 0 1],'EdgeAlpha',0.5);
        end
        
    else
        pipe = line(x,y,z); pipe.Color = 'm'; pipe.LineWidth = 3;
    end
end
view(camera);


xCOM = mean(X.p(1:3:end));
yCOM = mean(X.p(2:3:end));
Xrange = 0.5;%range(X.p(1:3:end));
Yrange = 0.5;%range(X.p(2:3:end));
zlim([-0.1 1])
xlim([xCOM-Xrange*1 xCOM+Xrange*1])
ylim([yCOM-Yrange*1 yCOM+Yrange*1])

axis square;


%% Check Validity
if size(C,2) ~= length(Verts)/3
    miss = abs(size(C,2) - length(Verts))
    warning(['Warning: Cable Connectivity Matrix has incorrect Dimensions: ',num2str(miss),' Nodes do not agree'])
    warning('size(C) == [# Cables, #Nodes]');
end

if size(R,2) ~= length(Verts)/3
    miss = abs(size(R,2) - length(Verts))
    warning(['Warning: Rod Connectivity Matrix has incorrect Dimensions: ',num2str(miss),' Nodes do not agree'])
    warning('size(R) == [# Rods, #Nodes]');
end

if all(sum(C,2)==0) ~= 1 
    warning('Warning: The row sum of C does not equal 0');
end

if all(sum(R,2)==0) ~= 1 
    warning('Warning: The row sum of R does not equal 0');
end

end