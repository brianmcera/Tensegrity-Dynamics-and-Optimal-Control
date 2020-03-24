function [omega,X,constraints,generalForces] = I_six_bar_model()

%% DEFINE ROD/CABLE PHYSICAL PARAMETERS
omega.rods.linear_velocity = 1*ones(6,1);
omega.rods.maxLength = 0.761*ones(6,1);
omega.rods.minLength = 0.4*ones(6,1);
omega.rods.constrained = [];

%mass matrix for node
omega.M = 0.3/2*ones(12,1);%.600/2*ones(12,1);

omega.cables.stiffness = 1.0*800*ones(24,1);
omega.cables.pretension = 1.0*30*ones(24,1);
omega.cables.linear_velocity = 0.10*ones(24,1);
omega.cables.unactuated = zeros(24,1);
omega.cables.maxLength = 0.6*ones(24,1);
omega.cables.minLength = 0.1*ones(24,1);

omega.cables.paired = [1,13;2,14;3,4;5,9;6,10;7,8;11,12;15,16;17,23;...
    18,24;19,21;20,22]; 
omega.cables.passive = [];
omega.cables.similar = [,];

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
omega.C = sparse(omega.C);

omega.R = ...;
  [-1     1     0     0     0     0     0     0     0     0     0     0;
    0     0     1    -1     0     0     0     0     0     0     0     0;
    0     0     0     0    -1     1     0     0     0     0     0     0;
    0     0     0     0     0     0     1    -1     0     0     0     0;
    0     0     0     0     0     0     0     0    -1     1     0     0;
    0     0     0     0     0     0     0     0     0     0     1     -1];
omega.R = sparse(omega.R);

%% DEFINE INITIAL NODAL POSITIONS/VELOCITIES

%Node XYZ positions
%current Rod Length - 0.6604
X.p = [
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

X.p(3:3:end) = X.p(3:3:end)-min(X.p(3:3:end))+0.1; %zero ground
omega.X0 = X.p;

%Node XYZ Velocities
X.pDOT = zeros(36,1);
%X.pDOT(3:3:end) = -15;

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

baseFloor = min(X.p(3:3:end))-0.1;
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

