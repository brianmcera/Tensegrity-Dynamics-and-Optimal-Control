function [omega,X,U,constraints,generalForces] = I_six_bar_model()

%% DEFINE ROD/CABLE PHYSICAL PARAMETERS
omega.rods.linear_velocity = 1*ones(6,1);
omega.rods.maxLength = 1.0*ones(6,1);
omega.rods.minLength = 0.4*ones(6,1);
omega.rods.constrained = [];

%mass matrix for node
omega.M = .125*ones(12,1);

omega.cables.stiffness = 1.0*400*ones(24,1);%375*ones(24,1);%
omega.cables.pretension = 1.0*50*ones(24,1);%25*ones(24,1);%
omega.cables.linear_velocity = 2.0*0.20*ones(24,1);
omega.cables.unactuated = zeros(24,1);
omega.cables.maxLength = 0.6*ones(24,1);
omega.cables.minLength = 0.1*ones(24,1);
omega.cables.paired = [];
% omega.cables.paired = [1,13;2,14;3,4;5,9;6,10;7,8;11,12;15,16;17,23;...
%    18,24;19,21;20,22]; %this is a paired-cable approach similar to one
%    proposed by Adrian at NASA
% omega.cables.passive = [1,2,3,5,7,8,10,11,12,13,14,15,17,20,21,...
%     22,23,24]; %actuated: 9,18,4,16,6,19
omega.cables.passive = [];

%% DEFINE CABLE/ROD CONNECTIVITY MATRICES
omega.C = ...
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

omega.R = ...;
  [-1     1     0     0     0     0     0     0     0     0     0     0;
    0     0     1    -1     0     0     0     0     0     0     0     0;
    0     0     0     0    -1     1     0     0     0     0     0     0;
    0     0     0     0     0     0     1    -1     0     0     0     0;
    0     0     0     0     0     0     0     0    -1     1     0     0;
    0     0     0     0     0     0     0     0     0     0     1     -1];

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

%Node XYZ Velocities
X.pDOT = zeros(36,1);

%% DEFINE INITIAL CABLE/ROD LENGTH INPUTS

%helper function handles
z = @(p,V) -kron(V'*V,eye(3))*p; %tension direction in p-basis
separationDist = @(p,V) sqrt(z(p,V)'*z(p,V)/2);

%initialize cable restlengths (untensioned)
C = omega.C;
for cable = 1:size(C,1)
    U.RL(cable) = separationDist(X.p,C(cable,:));
end
U.RL = U.RL'; %column vector

%initialize rod lengths
R = omega.R;
for rod = 1:size(R,1)
    U.L(rod) = separationDist(X.p,R(rod,:));
end
U.L = U.L'; %column vector

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

