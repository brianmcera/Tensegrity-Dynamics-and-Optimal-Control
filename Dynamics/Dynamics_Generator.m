function [nominalFcn,hFcns,jacobians,Gamma,genForces,constrForces,debugFcns] =...
    Dynamics_Generator(omega,constraints,generalForces)
%DYNAMICS_GENERATOR Generate the necessary dynamics equations for a given
%tensegrity topology. Receives input from tensegrity model codes in the
%'Models' folder.
%
%   INPUTS:
%   omega - a struct containing important tensegrity property subfields:
%       rods - 
%           stiffness - 
%           linVel - 
%       cables - 
%           stiffness - 
%           linVel - 
%           
%       mass - 
%       C - cable connectivity matrix
%       R - rod connectivity matrix
%       xyz - 
%       xyzDOT - 
%   constr - a struct relevant constraint information:
%       filename - matlab file in the 'Models/Constraints' folder
%       args - relevant input arguments to the {filename}.m function
%        plotInstr - 
%
%   OUTPUTS:
%   pDDOT - function handle for nonlinear pDDOT dynamics 
%   hFcns - function handles to helper functions useful for quickly
%   calculating the pDDOT dynamics (e.g., nodal separation calculation 
%   function)
%   jacobians - struct of Jacobians of f wrt p, pDOT, RL, and L
%   Gamma - handle to calculate vector of cable forces
%   genForces - handle to calculate vector of external general forces
%   constrForces - handle to calculate vector of constraint forces (e.g., 
%   rod constraints or stationary points)


%% Dynamic Constraints
%this section calculates forces due to the dynamic constraints defined in
%the {constraints} cell array. 
%{J}, the jacobian matrix of the constraints with respect to the nodal 
%positions, is used multiple times in the pDDOT calculation, and so its
%function handle is passed as a helper function.

%calculate costraint matrices (G, GDOT, J, JDOT) as functions of (X,U)
% G = @(X,U,hVars)[]; 
% GDOT = @(X,U,hVars)[];
% J = @(X,U,hVars)[];
% JDOT = @(X,U,hVars)[];

for i=1:length(constraints)
    if(i~=1)
        %evaluate constraint i, obtain anonymous functions and concatenate
        G = @(X,U,hVars)[G(X,U,hVars);feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'G')];
        GDOT = @(X,U,hVars)[GDOT(X,U,hVars);feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'GDOT')];
        J = @(X,U,hVars)[J(X,U,hVars);feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'dGdp')];
        JDOT = @(X,U,hVars)[JDOT(X,U,hVars);feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'dGDOTdp')];
        %{
        NEED TO ADD DG/GDOT WRT PDOT,RL,L
            %}
    else
        %evaluate constraint i, obtain anonymous functions and concatenate
        G = @(X,U,hVars) feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'G');
        GDOT = @(X,U,hVars) feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'GDOT');
        J = @(X,U,hVars) feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'dGdp');
        JDOT = @(X,U,hVars) feval(constraints{i}.filename,...
            X,U,hVars,omega,constraints{i}.args,'dGDOTdp');
        %{
        NEED TO ADD DG/GDOT WRT PDOT,RL,L
            %}
    end
end


%% Cable Forces in p-basis

%calculate cable forces (Gamma) as a function of (X,U,helperVariables)

C = omega.C;
R = omega.R;
n = size(omega.X0,1);

Gamma = @(X,U,hVars) zeros(n,1); %initialize cable forces as 0s
dGammadp = @(X,U,hVars) zeros(n,n);
dGammadRL = @(X,U,hVars) [];
Beta = 1e-4;

smax = @(x,b) (sqrt(x.^2+b^2)+x)/2; %smoothMax function
Alpha = @(X,U,z,i) sqrt(z{i}'*z{i}/2)-X.RL(i); %nodal separation difference

for i = 1:size(C,1)
    %cable {i} forces, represented in p-basis
    Gamma = @(X,U,hVars) Gamma(X,U,hVars)+...
        omega.cables.stiffness(i)*sqrt(2)*...
        smax(Alpha(X,U,hVars.z,i),Beta)*hVars.z{i}/norm(hVars.z{i});
    %Beta, above, is a scaling factor for the smoothMax (aka softPlus)
    
    %Jacobian of Gamma w.r.t. p state variable (derived by hand)
    dGammadp = @(X,U,hVars) dGammadp(X,U,hVars)+...
        omega.cables.stiffness(i)/sqrt(2)*(...
        diag(hVars.z{i}/norm(hVars.z{i}))*...
        repmat((1/sqrt(Alpha(X,U,hVars.z,i)^2+Beta^2)*...
        Alpha(X,U,hVars.z,i)+1)*...
        hVars.z{i}'/(sqrt(2)*norm(hVars.z{i}))*hVars.Chat{i},...
        n,1) + ...
        smax(Alpha(X,U,hVars.z,i),Beta)*...
        (eye(n)-...
        hVars.z{i}*hVars.z{i}'/norm(hVars.z{i})^2)/...
        norm(hVars.z{i})*hVars.Chat{i}...
        );
    
    %Jacobian of Gamma w.r.t. pDOT state variable
    dGammadpDOT = @(X,U,hVars) zeros(n,n); %cable damping ignored
    
    %Jacobian of Gamma w.r.t. RL input
    dGammadRL = @(X,U,hVars) [dGammadRL(X,U,hVars),...
        omega.cables.stiffness(i)/sqrt(2)*(hVars.z{i}/norm(hVars.z{i}))...
        *(-1/sqrt(Alpha(X,U,hVars.z,i)^2+Beta^2)*...
        Alpha(X,U,hVars.z,i)-1)];
    
    %Jacobian of Gamma w.r.t. L input
    dGammadL = @(X,U,hVars) zeros(n,size(R,1)); 
end


%% calculate General XYZ Forces
C = omega.C;
R = omega.R;
F = @(X,U,hVars) zeros(n,1);
dFdp = @(X,U,hVars) zeros(n,n);
dFdpDOT = @(X,U,hVars) zeros(n,n);
dFdRL = @(X,U,hVars) zeros(n,size(C,1));
dFdL = @(X,U,hVars) zeros(n,size(R,1));

for i=1:length(generalForces)
    %evaluate general force i, obtain anonymous functions and add together
    F = @(X,U,hVars) F(X,U,hVars)+ feval(generalForces{i}.filename,...
        X,U,omega,generalForces{i}.args,'genF'); 
    dFdp = @(X,U,hVars) dFdp(X,U,hVars)+ feval(generalForces{i}.filename,...
        X,U,omega,generalForces{i}.args,'dgenFdp'); 
    dFdpDOT = @(X,U,hVars) dFdpDOT(X,U,hVars)+ feval(generalForces{i}.filename,...
        X,U,omega,generalForces{i}.args,'dgenFdpDOT');
    dFdRL = @(X,U,hVars) dFdRL(X,U,hVars)+ feval(generalForces{i}.filename,...
        X,U,omega,generalForces{i}.args,'dgenFdRL');
    dFdL = @(X,U,hVars) dFdL(X,U,hVars)+ feval(generalForces{i}.filename,...
        X,U,omega,generalForces{i}.args,'dgenFdL');
end

%% prepare outputs

W = kron(diag(1./omega.M),eye(3)); %inverse mass matrix
Ks = 1e-4; %spring correcting feedback on G constraint forces
Kd = 1e-4; %damping correcting feedback on G constraint forces

genForces = F;

%nominal models
nominalFcn.pDDOT = @(X,U,hVars) W*(...
    -hVars.J'/(hVars.J*W*hVars.J')*...
    (JDOT(X,U,hVars)*X.pDOT-Ks*G(X,U,hVars)-Kd*GDOT(X,U,hVars))+... %constraint forces
    (eye(length(X.p))-hVars.J'/(hVars.J*W*hVars.J')*hVars.J*W)*...
    (Gamma(X,U,hVars)+F(X,U,hVars))...%cable and general external forces;
    );

nominalFcn.RLdot = @(X,U,hVars) U.RLdot; 

nominalFcn.Ldot = @(X,U,hVars) U.Ldot; 

%constraint forces
constrForces = @(X,U,hVars) (-hVars.J'/(hVars.J*W*hVars.J')*...
    (JDOT(X,U,hVars)*X.pDOT-Ks*G(X,U,hVars)-Kd*GDOT(X,U,hVars))+...
    (-hVars.J'/(hVars.J*W*hVars.J')*hVars.J*W)*...
    (Gamma(X,U,hVars)+F(X,U,hVars)));

%Jacobians
dpDDOTdp = @(X,U,hVars) W*(eye(length(X.p))-...
    hVars.J'/(hVars.J*W*hVars.J')*hVars.J*W)*...
    (dGammadp(X,U,hVars)+dFdp(X,U,hVars));

dpDDOTdpDOT = @(X,U,hVars) W*(...
    (-hVars.J'/(hVars.J*W*hVars.J')*JDOT(X,U,hVars)) + ...
    (eye(length(X.p))-hVars.J'/(hVars.J*W*hVars.J')*...
    hVars.J*W)*(dGammadpDOT(X,U,hVars)+dFdpDOT(X,U,hVars)));

dpDDOTdRL = @(X,U,hVars) W*(eye(length(X.p))-...
    hVars.J'/(hVars.J*W*hVars.J')*hVars.J*W)*...
    (dGammadRL(X,U,hVars)+dFdRL(X,U,hVars));

dpDDOTdL = @(X,U,hVars) W*(eye(length(X.p))-...
    hVars.J'/(hVars.J*W*hVars.J')*hVars.J*W)*...
    (dGammadL(X,U,hVars)+dFdL(X,U,hVars));


%aggregate into one struct of Jacobian function handles
%state jacobian function handles
jacobians.dpDDOTdp = dpDDOTdp;
jacobians.dpDDOTdpDOT = dpDDOTdpDOT;
jacobians.dpDDOTdRL = dpDDOTdRL;
jacobians.dpDDOTdL = dpDDOTdL;
jacobians.dpDDOTdRLdot = @(X,U,hVars) zeros(n,size(C,1));
jacobians.dpDDOTdLdot = @(X,U,hVars) zeros(n,size(R,1));

%RLdot input jacobian function handles <===================================== these are useless
jacobians.dRLdp = @(X,U,hVars) zeros(size(C,1),n);
jacobians.dRLdpDOT = @(X,U,hVars) zeros(size(C,1),n);
jacobians.dRLdRL = @(X,U,hVars) zeros(size(C,1),size(C,1));
jacobians.dRLdL = @(X,U,hVars) zeros(size(C,1),size(R,1));
jacobians.dRLdRLdot = @(X,U,hVars) eye(size(C,1));%omega.cableConstraintMatrix;
jacobians.dRLdLdot = @(X,U,hVars) zeros(size(C,1),size(R,1));
%Ldot input jacobian function handles
jacobians.dLdp = @(X,U,hVars) zeros(size(R,1),n);
jacobians.dLdpDOT = @(X,U,hVars) zeros(size(R,1),n);
jacobians.dLdRL = @(X,U,hVars) zeros(size(R,1),size(C,1));
jacobians.dLdL = @(X,U,hVars) zeros(size(R,1),size(R,1));
jacobians.dLdRLdot = @(X,U,hVars) zeros(size(R,1),size(C,1));
jacobians.dLdLdot = @(X,U,hVars) eye(size(R,1));%omega.rodConstraintMatrix;





%% Helper Functions
%these helper function handles are useful for calculating terms which come
%up repeatedly in the dynamics calculations (e.g., cable and rod nodal
%separation distances). By providing these helper functions, they can be
%calculated once per dynamics calculation and passed as input parameters,
%avoiding the need to repeat the same calculations multiple times.

%Note: Functions added here are likely identified by timing execution runs 
%and addressing function calls which bottleneck performance


R = omega.R;
%pre-calculate Kronecker products for anonymous functions
RkronMats = cell(size(R,1),1);
for i = 1:size(R,1)
    RkronMats{i} = -kron(R(i,:)'*R(i,:),eye(3));
end
for i = 1:size(R,1)
    %rod separation vectors
    hFcns.v{i} = @(X) RkronMats{i}*X.p; 
    %Rhat matrix multiplication repeatedly used in rod constraints (static)
    hFcns.RhatRhat{i} = (-kron(R(i,:)'*R(i,:),eye(3)))'*...
        (-kron(R(i,:)'*R(i,:),eye(3)));
end

C = omega.C;
CkronMats = cell(size(C,1),1);
for i = 1:size(C,1)
    CkronMats{i} = -kron(C(i,:)'*C(i,:),eye(3));
end
for i = 1:size(C,1)
    %cable separation vectors
    hFcns.z{i} = @(X) CkronMats{i}*X.p; 
    %Chat matrix repeatedly used in dGammadp (static)
    hFcns.Chat{i} = -kron(C(i,:)'*C(i,:),eye(3));
end

hFcns.J = J; %dGdp used multiple times in calculation of pDDOT


%% Create a struct of function handles passed 
%(not necessary for exexution but useful for debugging)
debugFcns.dGammadp = dGammadp;
debugFcns.dGammadRL = dGammadRL;
debugFcns.dFdp = dFdp;
debugFcns.dFdRL = dFdRL;

end
    

