function f = Gravity(X,U,omega,args,desFunc)
%FLOORFORCE Summary of this function goes here
%   Detailed explanation goes here

% %args fields/values
% gravity = args.gravity; 
% 
% %apply spring-damper floor forces
% %forces only apply to nodes in contact with ground
% %use differential smoothMax and smoothStep approximations
% genF = zeros(size(X.p),class(X.p));
% genF(3:3:end) = -omega.M.*gravity; 
% 
% %gravity is a conservative vector field
% dgenFdp = zeros(size(X.p,1),size(X.p,1),class(X.p));
% dgenFdpDOT = zeros(size(X.p,1),size(X.p,1),class(X.p));
% dgenFdRL = zeros(size(X.p,1),size(X.RL,1),class(X.p));
% dgenFdL = zeros(size(X.p,1),size(X.L,1),class(X.p));


%select which function to output for anonymous function handle
switch desFunc
    case 'genF'
        %args fields/values
        gravity = args.gravity;
        %apply spring-damper floor forces
        %forces only apply to nodes in contact with ground
        %use differential smoothMax and smoothStep approximations
        genF = zeros(size(omega.X0),class(X.p));
        genF(3:3:end) = -omega.M.*gravity;
        f = genF;
    case 'dgenFdp'
        %gravity is a conservative vector field
        dgenFdp = zeros(size(omega.X0,1),size(omega.X0,1),class(X.p));
        f = dgenFdp;
    case 'dgenFdpDOT'
        %gravity is a conservative vector field
        dgenFdpDOT = zeros(size(omega.X0,1),size(omega.X0,1),class(X.p));
        f = dgenFdpDOT;
    case 'dgenFdRL'
        dgenFdRL = zeros(size(omega.X0,1),size(omega.C,1),class(X.p));
        f = dgenFdRL;
    case 'dgenFdL'
        dgenFdL = zeros(size(omega.X0,1),size(omega.R,1),class(X.p));
        f = dgenFdL;
end

end

