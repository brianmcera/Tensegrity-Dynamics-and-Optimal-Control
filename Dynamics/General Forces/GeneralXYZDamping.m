function f = GeneralXYZDamping(X,U,omega,args,desFunc)
%FLOORFORCE Summary of this function goes here
%   Detailed explanation goes here

% %args fields/values
% damping = args.damping;
% 
% %apply general damping 
% genF = -damping*X.pDOT; 
% 
% %jacobian wrt pDOT
% dgenFdpDOT = diag(ones(length(X.pDOT),1)*(-damping));
% 
% 
% %node positions, cable lengths, and rod lengths do not affect XYZ damping
% dgenFdp = zeros(size(X.p,1),size(X.p,1),class(X.p));
% dgenFdRL = zeros(size(X.p,1),size(X.RL,1),class(X.p));
% dgenFdL = zeros(size(X.p,1),size(X.L,1),class(X.p));


%select which function to output for anonymous function handle
switch desFunc
    case 'genF'
        %args fields/values
        damping = args.damping;    
        %apply general damping
        genF = -damping*X.pDOT;
        f = genF;
    case 'dgenFdp'
        %node positions, cable lengths, and rod lengths do not affect XYZ damping
        dgenFdp = zeros(size(omega.X0,1),size(omega.X0,1),class(X.p));
        f = dgenFdp;
    case 'dgenFdpDOT'
        %args fields/values
        damping = args.damping;
        %jacobian wrt pDOT
        dgenFdpDOT = diag(ones(size(omega.X0,1),1)*(-damping));
        f = dgenFdpDOT;
    case 'dgenFdRL'
        %node positions, cable lengths, and rod lengths do not affect XYZ damping
        dgenFdRL = zeros(size(omega.X0,1),size(omega.C,1),class(X.p));
        f = dgenFdRL;
    case 'dgenFdL'
        %node positions, cable lengths, and rod lengths do not affect XYZ damping
        dgenFdL = zeros(size(omega.X0,1),size(omega.R,1),class(X.p));
        f = dgenFdL;
end

end

