function [omega,X,U,constraints,generalForces]=tensegrity_manipulator_model()
%% Tensegrity Manipulator Modeling

% This code should take in no arguments and should output the
% characteristic properties for the rods and cables, the vertices of the
% structure, the constraints applied to the structure, and the connectivity
% matrix.
% The output variables are of a structure type. This tensegrity structure
% has 7 rods and this code only produces one structure.

NumberRods=4;
NumberCables = 15;
%% Define Physical Parameters

%This section will define all of the physical parameters of the lower
%tensegrity structure (6 Bar). All units are given in SI units. Define all
%vectors as column vectors.

% Rods
omega.rods.linear_velocity = ones(NumberRods,1);
omega.rods.maxLength   = ones(NumberRods,1);    % meters
omega.rods.minLength     = 0.4*ones(NumberRods,1);
omega.rods.constrained = [];
omega.rods.Radius     = 0.02*ones(NumberRods,1); % meters
omega.M = 0.047*ones(NumberRods*2,1); %mass for each node
omega.rods.Stiffness  = 1e5*ones(NumberRods,1);  % N/m
%spacing between planar rod-pairs (6-bar topology specific) (Stolen from Brian's Code)
node_Separation = (-omega.rods.maxLength(1)+sqrt(omega.rods.maxLength(1).^2+4*omega.rods.maxLength(1).^2))/2;
planeTranslate = @(l,x) sqrt(l.^2 - x.^2); % Amount to translate the plane... This should output 1x3 vector [dx,dy,dz]

% Cables [Change cable matrix dimensions]
omega.cables.minLength     = 0.2*ones(NumberCables,1);  % meters
omega.cables.maxLength     = 0.9*ones(NumberCables,1);  % meters
omega.cables.linear_velocity = 0.1*ones(NumberCables,1);  % meters/second
omega.cables.pretension    = 20*ones(NumberCables,1);   % Newtons
omega.cables.stiffness     = 1e3*ones(NumberCables,1);  %N/m
omega.cables.constrained = [];
omega.cables.paired = [];
%% Nodal Locations
% This section defines the location for each of the nodes of the structure

%groundPlane=[0,0,0;1,0,0;0,1,0]; % Establish a Ground Plane NOT USED

baseShape   = 0;

Verts=[];
if baseShape == 0
    theta=deg2rad(60);

    %First Layer
    Verts(1:3,:)=[0,0,0;node_Separation,0,0;
        node_Separation*cos(theta),node_Separation*sin(theta),0];
    layer1=planeTranslate(omega.rods.maxLength(1:3),node_Separation*ones(1,3)');
    
    %Second Layer
    Verts(4:6,:)=[0,0,layer1(1);node_Separation,0,layer1(2);
        node_Separation*cos(theta),node_Separation*sin(theta),layer1(3)];
    
    %     for n=4:6
    %         Verts(n+3,:)=(Verts(n,:)+Verts(n+1,:))/2;
    %         if n==6
    %             Verts(n+3,:)=(Verts(n-2,:)+Verts(n,:))/2;
    %         end
    %     end
    %     layer2=planeTranslate(omega.rods.maxLength(4:6),(Verts(8,:)-Verts(9,:)));
    
    %Manipulator
    center = [node_Separation/2,((sqrt(3)/2)*node_Separation)/3,0];
    Verts(7,:) = [center(1:2),0.5*mean(layer1)];
    Verts(8,:) = [center(1:2),0.5*mean(layer1)+omega.rods.maxLength(4)];
    
    % Cable Connections
    omega.cables.pairs = [1,2;1,3;3,2;1,4;3,6;2,5;6,8;8,5;4,6;4,5;4,8;4,7;5,6;5,7;6,7];
    omega.rods.pairs = [1,5;2,6;3,4;7,8];
    %% Add this section for the second tier
    
%     %Third Layer
%     for n=10:12
%         Verts(n,:)=Verts(n-3,:);
%     end
%     Verts(10:12,:)=Verts(10:12,:)+[0,0,layer2(1);0,0,layer2(1);0,0,layer2(1)];
%     
%     %Manipulator
%     center = [node_Separation/2,((sqrt(3)/2)*node_Separation)/3,0];
%     Verts(13,:) = [center(1:2),mean(layer1+layer2)];
%     Verts(14,:) = [center(1:2),mean(layer1+layer2)+omega.rods.maxLength(7)];
%     
%     % Steaks
%     Verts(15:17,:) = Verts(1:3,:);
%     Verts(15,1) = Verts(15,1)-1; Verts(16,1) = Verts(16,1)+1;
%     Verts(17,2) =Verts(17,2)+1;
%     Verts(18:20,:) = Verts(7:9,:); Verts(18:20,3)=Verts(18:20,3)-layer1(3);
%     Verts(18,1:2) = Verts(18,1:2)+[0,-1];
%     Verts(19,1:2) = Verts(19,1:2)+[1,1];
%     Verts(20,1:2) = Verts(20,1:2)+[-1,1];
%     Verts(21:23,:) = Verts(18:20,:);
%     Verts(21,1:2) = Verts(21,1:2)+[0.15,0.5];
%     Verts(22,1:2) = Verts(22,1:2)+[-0.5,-0.75];
%     Verts(23,1:2) = Verts(23,1:2)+[0.5,-0.75];
    
    % Cable Connections
%     omega.cables.pairs = [1,2;1,3;3,2;1,4;3,6;2,5;6,8;8,5;6,9;4,9;4,7;5,7;8,7;8,9;
%         7,9;10,11;10,12;11,12;7,10;9,12;8,11;10,13;11,13;12,13;15,4;17,6;
%         16,5;18,10;19,11;20,12;21,14;22,14;23,14];
%     omega.rods.pairs = [1,5;2,6;3,4;7,12;8,10;9,11;13,14];
end
omega.cables.pairs = sortrows(omega.cables.pairs);
omega.rods.pairs = sortrows(omega.rods.pairs);
for n=1:length(omega.cables.pairs)
    if omega.cables.pairs(n,1)>omega.cables.pairs(n,2)
        a=omega.cables.pairs(n,2);
        omega.cables.pairs(n,2)=omega.cables.pairs(n,1);
        omega.cables.pairs(n,1)=a;
    end
end
% structurePlot(cable,rod,Verts);
% set(gcf,'Position',[641.6667   41.6667  638.6667  599.3333]);
%% Connectivity Matrix
%//////////////////////////////////////////////////////////////////////////
%CABLE CONNECTIVITY MATRIX/////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
% C([cable number],[nodeA#])=1;C([cable number],[nodeB#])=-1;
% ex: C(1,1)=1;C(1,5)=-1;

omega.C = zeros(length(omega.cables.pairs),length(Verts)); %(number of cables)x(number of nodes)
for n=1:size(omega.C,1) %down the row of C == Cable number
    omega.C(n,omega.cables.pairs(n,:)) = [1,-1];
end

%//////////////////////////////////////////////////////////////////////////
% ROD CONNECTIVITY MATRIX (EACH ROW IS A ROD)//////////////////////////////
%//////////////////////////////////////////////////////////////////////////
% i.e. R([rod number],[nodeA#])=1; R([rod number],[nodeB#])=-1;

% R-Matrix dimension: (number of rods)x(number of nodes)

omega.R=zeros(length(omega.rods.pairs),length(Verts));
for n=1:size(omega.R,1) %down the row of R == Rod number
    omega.R(n,omega.rods.pairs(n,:)) = [1,-1];
end

%% Save Verts to X structure
X.p = NaN*ones(size(Verts,1)*size(Verts,2),1);
X.p(1:3:end) =  Verts(:,1);
X.p(2:3:end) =  Verts(:,2);
X.p(3:3:end) =  Verts(:,3);

X.pDOT = zeros(size(X.p));

%% Plot Check
% structurePlot(X.p,omega.C,omega.R,[]);
% set(gcf,'Position',[641.6667   41.6667  638.6667  599.3333]);

%%  DEFINE INITIAL CABLE/ROD LENGTH INPUTS
%helper function handles
z = @(p,V) -kron(V'*V,eye(3))*p; %tension direction in p-basis
% z = @(p,V) -kron(V'*V,eye(1))*p;
separationDist = @(p,V) sqrt(z(p,V)'*z(p,V)/2);

C = omega.C;
for cable = 1:size(C,1) %Need to Fix
    x = separationDist(X.p,C(cable,:));
    U.RL(cable) = norm(x);
%     U.RL(cable) = separationDist(X.p,C(cable,:));
end
U.RL = U.RL'; %column vector

%initialize rod lengths
R = omega.R;
for rod = 1:size(R,1) %Need to Fix
    x = separationDist(X.p,R(rod,:));
    U.L(rod) = norm(x);
%     U.L(rod) = separationDist(X.p,R(rod,:));
end
U.L = U.L'; %column vector
%% DEFINE MODEL CONSTRAINTS////////////////////////////////////////////////
%constraints is a cell array, where constraints{i} has two fields:
%filename: the filename of the corresponding m-file in the
%Models/Constraints folder.
%args: any file-specific input arguments necessary to pass to the m-file
%which defines the constraints

constraints{1}.filename = 'RodConstraints';
constraints{1}.args = [];

constraints{2}.filename = 'StationaryConstraints';
constraints{2}.args.p0 = X.p;
constraints{2}.args.constrain = zeros(size(X.p));
constraints{2}.args.constrain(1:9) = 1; %constrain x,y,z of first 3 nodes


%% DEFINE GENERAL EXTERNAL FORCES/////////////////////////////////////////
%generalForces is a cell array, where generalForces{i} has two fields:
%filename: the filename of the corresponding m-file in the
%Dynamics/General Forces folder.
%args: any file-specific input arguments necessary to pass to the m-file
%which defines the constraints

generalForces{1}.filename = 'Gravity';
generalForces{1}.args.gravity = 9.81;

generalForces{2}.filename = 'GeneralXYZDamping';
generalForces{2}.args.damping = 1e-4;

% generalForces{4}.filename = 'FloorForceHorizontal_ViscousFriction';
% generalForces{4}.args.baseFloor = baseFloor;
% generalForces{4}.args.damping = 5e1;
% generalForces{4}.args.Beta = 1e-6;
% 
% baseFloor = min(X.p(3:3:end))-0.1;
% generalForces{1}.filename = 'FloorForceVertical';
% generalForces{1}.args.baseFloor = baseFloor;
% generalForces{1}.args.stiffness = 5e2;
% generalForces{1}.args.damping = 3e1;
% generalForces{1}.args.Beta = 1e-6;





%% Constraints (NEED TO FIX)
% % Contrain all of the pins on the ground. Constraint matrix should be the
% % same size as the connectivity matrix C. 
% 
%         
% constraint = zeros(size(Verts,1),3);
% for n=1:length(constraint)
%     if n<=3
%         constraint(n,:) = ones(3,1);
%     elseif n>=15 && n<=17
%         constraint(n,:) = [1,0,1];
%     elseif n>=18 & n<=20
%         constraint(n,:) = [1,0,0];
%     elseif n>=21 && n<=23
%         constraint(n,:) = [0,0,1];
%     end
% end
