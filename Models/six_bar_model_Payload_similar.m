function [omega,X,U,constraints,generalForces] = I_six_bar_model_Payload_similar()

%% DEFINE ROD/CABLE PHYSICAL PARAMETERS
omega.rods.linear_velocity = 1*ones(7,1);
omega.rods.maxLength = 1.0*ones(7,1);
omega.rods.minLength = 0*ones(7,1);
omega.rods.constrained = [];

%mass vector for node
omega.M = [.125*ones(12,1);
    1.5;1.5];

omega.cables.stiffness = [200*ones(24,1);800*ones(12,1)];%375*ones(24,1);%
omega.cables.pretension = 30*ones(36,1);%25*ones(24,1);%
omega.cables.linear_velocity = 0.10*ones(36,1);
omega.cables.unactuated = zeros(36,1);
omega.cables.maxLength = 1.00*ones(36,1);
omega.cables.minLength = 0.05*ones(36,1);
omega.cables.paired = [,];
% omega.cables.paired = [1,13;2,14;3,4;5,9;6,10;7,8;11,12;15,16;17,23;...
%    18,24;19,21;20,22]; %this is a paired-cable approach similar to one
% %    proposed by Adrian at NASA
omega.cables.passive = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,...
    22,23,24];
% omega.cables.passive = [];
% omega.cables.similar = [,];
% omega.cables.paired = [25,28;26,27;29,32;30,31;33,36;34,35];
omega.cables.similar = [25,28;26,27;29,32;30,31;33,36;34,35];

%% DEFINE CABLE/ROD CONNECTIVITY MATRICES
omega.C = ...
   [0     1     0     0     0    -1     0     0     0     0     0     0     0     0;
    0     1     0     0    -1     0     0     0     0     0     0     0     0     0;
    0     1     0     0     0     0     0     0     0    -1     0     0     0     0;
    0     1     0     0     0     0     0     0     0     0     0    -1     0     0;
    1     0     0     0     0     0    -1     0     0     0     0     0     0     0;
    1     0     0     0     0     0     0    -1     0     0     0     0     0     0;
    1     0     0     0     0     0     0     0     0    -1     0     0     0     0;
    1     0     0     0     0     0     0     0     0     0     0    -1     0     0;
    0     0     1     0     0     0    -1     0     0     0     0     0     0     0;
    0     0     1     0     0     0     0    -1     0     0     0     0     0     0;
    0     0     1     0     0     0     0     0     0     0    -1     0     0     0;
    0     0     1     0     0     0     0     0    -1     0     0     0     0     0;
    0     0     0     1     0    -1     0     0     0     0     0     0     0     0;
    0     0     0     1    -1     0     0     0     0     0     0     0     0     0;
    0     0     0     1     0     0     0     0     0     0    -1     0     0     0;
    0     0     0     1     0     0     0     0    -1     0     0     0     0     0;
    0     0     0     0     0     1     0     0     0    -1     0     0     0     0;
    0     0     0     0     0     1     0     0    -1     0     0     0     0     0;
    0     0     0     0     1     0     0     0     0     0     0    -1     0     0;
    0     0     0     0     1     0     0     0     0     0    -1     0     0     0;
    0     0     0     0     0     0     1     0     0     0     0    -1     0     0;
    0     0     0     0     0     0     1     0     0     0    -1     0     0     0;
    0     0     0     0     0     0     0     1     0    -1     0     0     0     0;
    0     0     0     0     0     0     0     1    -1     0     0     0     0     0;
    1     0     0     0     0     0     0     0     0     0     0     0    -1     0;
    0     1     0     0     0     0     0     0     0     0     0     0    -1     0;
    0     0     1     0     0     0     0     0     0     0     0     0    -1     0;
    0     0     0     1     0     0     0     0     0     0     0     0    -1     0;
    0     0     0     0     1     0     0     0     0     0     0     0    -1     0;
    0     0     0     0     0     1     0     0     0     0     0     0    -1     0;
    0     0     0     0     0     0     1     0     0     0     0     0     0    -1;
    0     0     0     0     0     0     0     1     0     0     0     0     0    -1;
    0     0     0     0     0     0     0     0     1     0     0     0     0    -1;
    0     0     0     0     0     0     0     0     0     1     0     0     0    -1;
    0     0     0     0     0     0     0     0     0     0     1     0     0    -1;
    0     0     0     0     0     0     0     0     0     0     0     1     0    -1];
omega.C = sparse(omega.C);

omega.R = ...;
  [-1     1     0     0     0     0     0     0     0     0     0     0     0     0;
    0     0     1    -1     0     0     0     0     0     0     0     0     0     0;
    0     0     0     0    -1     1     0     0     0     0     0     0     0     0;
    0     0     0     0     0     0     1    -1     0     0     0     0     0     0;
    0     0     0     0     0     0     0     0    -1     1     0     0     0     0;
    0     0     0     0     0     0     0     0     0     0     1    -1     0     0;
    0     0     0     0     0     0     0     0     0     0     0     0     1    -1];
omega.R = sparse(omega.R);

%% DEFINE INITIAL NODAL POSITIONS/VELOCITIES

%Node XYZ positions
X.p =[
   -0.0280
   -0.3453
   -0.0662
   -0.1761
    0.1217
    0.2803
    0.1761
   -0.1217
   -0.2803
    0.0280
    0.3453
    0.0662
   -0.2850
    0.1969
   -0.0662
    0.1935
    0.0917
    0.2803
   -0.1935
   -0.0917
   -0.2803
    0.2850
   -0.1969
    0.0662
    0.3130
    0.1484
   -0.0662
   -0.0173
   -0.2134
    0.2803
    0.0173
    0.2134
   -0.2803
   -0.3130
   -0.1484
    0.0662
    ];

payloadLength = 0.005;
COM = [mean(X.p(1:3:end));mean(X.p(2:3:end));mean(X.p(3:3:end))];

%add Payload nodes
X.p = [X.p;COM+[0,0,payloadLength/2]';COM-[0,0,payloadLength/2]'];

%zero floor
X.p(3:3:end) = X.p(3:3:end) - min(X.p(3:3:end));

omega.X0 = X.p;
    
%Node XYZ Velocities
X.pDOT = zeros(42,1);

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


%define cable constraint matrices for {paired,similar,passive} constraints
%matrices are used to constrain LQR controller
constrainedCables = [];
if(numel(omega.cables.paired)>0)
    pairedTranspose = omega.cables.paired';
    constrainedCables = [constrainedCables;pairedTranspose(:)];
end
if(numel(omega.cables.similar)>0)
    similarTranspose = omega.cables.similar';
    constrainedCables = [constrainedCables;similarTranspose(:)];
end

if(numel(unique([constrainedCables;omega.cables.passive']))<...
        numel([constrainedCables;omega.cables.passive']))
    %for now, we don't allow cables to have two simultaneous constraints
    error('Some cables are simultaneously paired/similar/passive')
end
cableConstraintMatrix = eye(size(omega.C,1));
for i = 1:size(omega.cables.paired,1)
    cableConstraintMatrix(omega.cables.paired(i,2),...
        omega.cables.paired(i,1)) = -1;
end
for i = 1:size(omega.cables.similar,1)
    cableConstraintMatrix(omega.cables.similar(i,2),...
        omega.cables.similar(i,1)) = 1;
end
%delete constraint cable columns for 2nd cable in each pair
cableConstraintMatrix(:,[constrainedCables(2:2:end);omega.cables.passive']) = [];
omega.cableConstraintMatrix = cableConstraintMatrix; %record in struct

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

baseFloor = min(X.p(3:3:end))-1e-3;
omega.generalForces{1}.filename = 'FloorForceVertical_XIncline';
omega.generalForces{1}.args.baseFloor = baseFloor;
omega.generalForces{1}.args.stiffness = 5e4; %5e2
omega.generalForces{1}.args.damping = 3e2; %3e1
omega.generalForces{1}.args.Beta = 1e-2;
omega.generalForces{1}.args.incline = 0;

omega.generalForces{2}.filename = 'Gravity';
omega.generalForces{2}.args.gravity = 9.81;

omega.generalForces{3}.filename = 'GeneralXYZDamping';
omega.generalForces{3}.args.damping = 1e-0;

omega.generalForces{4}.filename = 'FloorForceHorizontal_ViscousFriction';
omega.generalForces{4}.args.baseFloor = baseFloor;
omega.generalForces{4}.args.damping = 5e1;
omega.generalForces{4}.args.Beta = 1e-2;



end

