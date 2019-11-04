function [omega,X,constraints,generalForces] = I_six_bar_model_CenterCableConnection()

%% DEFINE ROD/CABLE PHYSICAL PARAMETERS
omega.rods.linear_velocity = 1*ones(7,1);
omega.rods.maxLength = [0.661*ones(6,1);0.1];
omega.rods.minLength = [0.4*ones(6,1);1e-4];
omega.rods.constrained = []; %unused for now

%mass matrix for node
omega.M = [.100/2*ones(12,1);.05*ones(24,1);1.5*ones(2,1)]
omega.cables.stiffness = [1.0*400*ones(48,1);1.0*1200*ones(24,1)];%375*ones(24,1);%
omega.cables.pretension = [1.0*20*ones(48,1);1.0*30*ones(24,1)];%25*ones(24,1);%
omega.cables.linear_velocity = [1.0*0.30*ones(48,1);1.0*0.30*ones(24,1)];
%omega.cables.unactuated = zeros(24,1);
omega.cables.maxLength = 0.6*ones(72,1);
omega.cables.minLength = 0.05*ones(72,1);
omega.cables.paired = [];
% omega.cables.paired = [1,13;2,14;3,4;5,9;6,10;7,8;11,12;15,16;17,23;...
%    18,24;19,21;20,22]; %this is a paired-cable approach similar to one proposed by Adrian at NASA
% omega.cables.paired = [4,1;7,5;16,14;11,10;2,20;13,17;6,24;9,21;18,12;23,3;22,15;19,8]; %Squishy paired
% omega.cables.passive = [1,2,3,5,7,8,10,11,12,13,14,15,17,20,21,...
%     22,23,24]; %actuated: 9,18,4,16,6,19
% omega.cables.passive = [3,5,6,7,11,13,14,15,19,20,23,24]; 
% %forward (-X direction):1,4,21,9,12,18; backward: 17,16,10,22,8,2
% omega.cables.passive = [2,3,5,6,7,8,10,11,13,14,15,16,17,19,20,22,23,24]; 
%forward (-X direction):1,4,21,9,12,18;
omega.cables.passive = [];
omega.cables.passive = 1:48;

%% DEFINE NOMINAL CABLE/ROD CONNECTIVITY MATRICES
Cbar = ...
   [0     1     0     0     0    -1     0     0     0     0     0     0;
    0     1     0     0    -1     0     0     0     0     0     0     0;
    0     1     0     0     0     0     0     0     0    -1     0     0;
    0     1     0     0     0     0     0     0     0     0     0    -1;
    1     0     0     0     0     0    -1     0     0     0     0     0;
    1     0     0     0     0     0     0    -1     0     0     0     0;
    1     0     0     0     0     0     0     0     0    -1     0     0;
    1     0     0     0     0     0     0     0     0     0     0    -1;
    0     0     1     0     0     0    -1     0     0     0     0     0;
    0     0     1     0     0     0     0    -1     0     0     0     0;
    0     0     1     0     0     0     0     0     0     0    -1     0;
    0     0     1     0     0     0     0     0    -1     0     0     0;
    0     0     0     1     0    -1     0     0     0     0     0     0;
    0     0     0     1    -1     0     0     0     0     0     0     0;
    0     0     0     1     0     0     0     0     0     0    -1     0;
    0     0     0     1     0     0     0     0    -1     0     0     0;
    0     0     0     0     0     1     0     0     0    -1     0     0;
    0     0     0     0     0     1     0     0    -1     0     0     0;
    0     0     0     0     1     0     0     0     0     0     0    -1;
    0     0     0     0     1     0     0     0     0     0    -1     0;
    0     0     0     0     0     0     1     0     0     0     0    -1;
    0     0     0     0     0     0     1     0     0     0    -1     0;
    0     0     0     0     0     0     0     1     0    -1     0     0;
    0     0     0     0     0     0     0     1    -1     0     0     0];

Rbar = ...;
  [-1     1     0     0     0     0     0     0     0     0     0     0;
    0     0     1    -1     0     0     0     0     0     0     0     0;
    0     0     0     0    -1     1     0     0     0     0     0     0;
    0     0     0     0     0     0     1    -1     0     0     0     0;
    0     0     0     0     0     0     0     0    -1     1     0     0;
    0     0     0     0     0     0     0     0     0     0     1     -1];

%% DEFINE INITIAL NODAL POSITIONS/VELOCITIES

%Node XYZ positions
%current Rod Length - 0.6604
Verts = [
    0.1415
   -0.1884
   -0.3085
   -0.3785
   -0.0459
    0.0728
    0.3785
    0.0459
   -0.0728
   -0.1415
    0.1884
    0.3085
   -0.2290
    0.3049
   -0.0728
   -0.0924
   -0.2168
    0.3085
    0.0924
    0.2168
   -0.3085
    0.2290
   -0.3049
    0.0728
    0.2339
    0.0283
    0.3085
   -0.1495
   -0.3507
   -0.0728
    0.1495
    0.3507
    0.0728
   -0.2339
   -0.0283
   -0.3085
   ];

%Node XYZ Velocities
X.pDOT = zeros(38*3,1);

%% FORM C,R MATRICES

%calculate midpoints
midpoints_XYZ = kron(abs(Cbar),eye(3))*Verts/2;
X.p = [Verts;midpoints_XYZ]; %append new midpoint nodes

%cable connection matrix
omega.C = [];
for i = 1:size(Cbar,1)
    k = find(Cbar(i,:));
    appendVec = zeros(1,size(Cbar,1)); %number of added midpoint nodes is equal to the number of original cables
    appendVec(i) = -1;
    for j = 1:numel(k)
        leftVec = zeros(1,size(Cbar,2));
        leftVec(k(j)) = 1;
        omega.C = [omega.C;
            leftVec,appendVec];
    end
end


%add payload nodes
payloadLength = 0.005;
COM = [mean(X.p(1:3:end));mean(X.p(2:3:end));mean(X.p(3:3:end))];
X.p = [X.p;COM+[0,0,payloadLength/2]';COM-[0,0,payloadLength/2]'];

X.p(3:3:end) = X.p(3:3:end)-min(X.p(3:3:end))+0.1; %zero ground
omega.X0 = X.p;

omega.C = [omega.C,zeros(size(omega.C,1),2)];

omega.R = [Rbar,zeros(size(Rbar,1),size(Cbar,1)),zeros(size(Rbar,1),2)];
omega.R = [omega.R;zeros(1,size(omega.R,2))];
omega.R(end,end-1:end) = [1 -1];

for i = 1:size(Cbar,1)
    %append new row
    omega.C = [omega.C;zeros(1,size(omega.C,2))];
    if(mod(i,2)==0)
        omega.C(end,end) = -1;
    else
        omega.C(end,end-1) = -1;
    end
    omega.C(end,size(Cbar,2)+i) = 1;
end
    
        


%% DEFINE INITIAL CABLE/ROD LENGTH INPUTS

%helper function handles
z = @(p,V) -kron(V'*V,eye(3))*p; %tension direction in p-basis
separationDist = @(p,V) sqrt(z(p,V)'*z(p,V)/2);

%initialize cable restlengths (untensioned)
C = omega.C;
for cable = 1:size(C,1)
    X.RL(cable) = separationDist(X.p,C(cable,:));
end
X.RL = X.RL'; %column vector

%initialize rod lengths
R = omega.R;
for rod = 1:size(R,1)
    X.L(rod) = separationDist(X.p,R(rod,:));
end
X.L = X.L'; %column vector

omega.cableConstraintMatrix = eye(size(omega.C,1));
omega.rodConstraintMatrix = eye(size(omega.R,1));
%% DEFINE MODEL CONSTRAINTS////////////////////////////////////////////////
%constraints is a cell array, where constraints{i} has two fields:
%filename: the filename of the corresponding m-file in the
%Models/Constraints folder.
%args: any file-specific input arguments necessary to pass to the m-file
%which defines the constraints

omega.constraints{1}.filename = 'RodConstraints';
omega.constraints{1}.args = [];


%% DEFINE GENERAL EXTERNAL FORCES/////////////////////////////////////////
%generalForces is a cell array, where generalForces{i} has two fields:
%filename: the filename of the corresponding m-file in the
%Dynamics/General Forces folder.
%args: any file-specific input arguments necessary to pass to the m-file
%which defines the constraints

baseFloor = min(X.p(3:3:end))-0.1;
omega.generalForces{1}.filename = 'FloorForceVertical_XIncline';
omega.generalForces{1}.args.baseFloor = baseFloor;
omega.generalForces{1}.args.stiffness = 5e3; %5e2
omega.generalForces{1}.args.damping = 3e2; %3e1
omega.generalForces{1}.args.Beta = 1e-2;
omega.generalForces{1}.args.incline = 0;

omega.generalForces{2}.filename = 'Gravity';
omega.generalForces{2}.args.gravity = 9.81;

omega.generalForces{3}.filename = 'GeneralXYZDamping';
omega.generalForces{3}.args.damping = 1e-4;

omega.generalForces{4}.filename = 'FloorForceHorizontal_ViscousFriction';
omega.generalForces{4}.args.baseFloor = baseFloor;
omega.generalForces{4}.args.damping = 5e1;
omega.generalForces{4}.args.Beta = 1e-2;



end

