function f = FloorForceHorizontal_ViscousFriction(X,U,omega,args,desFunc)
%FLOORFORCE Summary of this function goes here
%   Detailed explanation goes here

% %args fields/values
% baseFloor = args.baseFloor; %floor height Z-value
% damping = args.damping; %floor damping
% Beta = args.Beta; %parameter for softMax/softStep
% 
% %sStep (smooth approximation to 0-1 step)
% sStep = @(x,Beta) 1/2*((x.^2+Beta^2).^(-1/2).*x+1);
% 
% %apply spring-damper floor forces
% %forces only apply to nodes in contact with ground
% %use differential smoothMax and smoothStep approximations
% genF = zeros(size(X.p),class(X.p));
% 
% ztilde = baseFloor - X.p(3:3:end);
% 
% %x-friction
% genF(1:3:end) = ...
%     sStep(ztilde,Beta).*...
%     (-X.pDOT(1:3:end))*damping; 
% 
% %y-friction
% genF(2:3:end) = ...
%     sStep(ztilde,Beta).*...
%     (-X.pDOT(2:3:end))*damping; 
% 
% %dFx/dp
% dgenFxdp = 1/2*(-(ztilde.^2+(Beta)^2).^(-3/2).*(ztilde.^2)+...
%     (ztilde.^2+(Beta)^2).^(-1/2))*damping.*(-X.pDOT(1:3:end))*(-1);
% dgenFdp = diag(kron(dgenFxdp,[1 0 0]'));
% 
% %dFy/dp
% dgenFydp = 1/2*(-(ztilde.^2+(Beta)^2).^(-3/2).*(ztilde.^2)+...
%     (ztilde.^2+(Beta)^2).^(-1/2))*damping.*(-X.pDOT(2:3:end))*(-1);
% dgenFdp = dgenFdp + diag(kron(dgenFydp,[0 1 0]'));
% 
% %dFx/dpDOT (same for dFy/dpDOT)
% dgenFdpDOT = -1/2*((ztilde.^2+(Beta)^2).^(-1/2).*ztilde+1)*damping;
% dgenFdpDOT = diag(kron(dgenFdpDOT,[1 1 0]'));
% 
% %cable/rod inputs do not affect viscous XY force from floor
% dgenFdRL = zeros(size(X.p,1),size(X.RL,1),class(X.p));
% dgenFdL = zeros(size(X.p,1),size(X.L,1),class(X.p));


%select which function to output for anonymous function handle
switch desFunc
    case 'genF'
        %args fields/values
        baseFloor = args.baseFloor; %floor height Z-value
        damping = args.damping; %floor damping
        Beta = args.Beta; %parameter for softMax/softStep
        
        %sStep (smooth approximation to 0-1 step)
        sStep = @(x,Beta) 1/2*((x.^2+Beta^2).^(-1/2).*x+1);
        
        %apply spring-damper floor forces
        %forces only apply to nodes in contact with ground
        %use differential smoothMax and smoothStep approximations
        genF = zeros(size(omega.X0,1),1,class(X.p));
        ztilde = baseFloor - X.p(3:3:end); 
        %x-friction
        genF(1:3:end) = ...
            sStep(ztilde,Beta).*...
            (-X.pDOT(1:3:end))*damping;
        %y-friction
        genF(2:3:end) = ...
            sStep(ztilde,Beta).*...
            (-X.pDOT(2:3:end))*damping;
        
        f = genF;
    case 'dgenFdp'
        %args fields/values
        baseFloor = args.baseFloor; %floor height Z-value
        damping = args.damping; %floor damping
        Beta = args.Beta; %parameter for softMax/softStep
                
        ztilde = baseFloor - X.p(3:3:end);
        
        %dFx/dp
        dgenFxdp = 1/2*(-(ztilde.^2+(Beta)^2).^(-3/2).*(ztilde.^2)+...
            (ztilde.^2+(Beta)^2).^(-1/2))*damping.*(-X.pDOT(1:3:end))*(-1);
        dgenFdp = diag(kron(dgenFxdp,[1 0 0]'));     
        %dFy/dp
        dgenFydp = 1/2*(-(ztilde.^2+(Beta)^2).^(-3/2).*(ztilde.^2)+...
            (ztilde.^2+(Beta)^2).^(-1/2))*damping.*(-X.pDOT(2:3:end))*(-1);
        dgenFdp = dgenFdp + diag(kron(dgenFydp,[0 1 0]'));
        
        f = dgenFdp;
    case 'dgenFdpDOT'
        %args fields/values
        baseFloor = args.baseFloor; %floor height Z-value
        damping = args.damping; %floor damping
        Beta = args.Beta; %parameter for softMax/softStep
        
        ztilde = baseFloor - X.p(3:3:end);
                
        %dFx/dpDOT (same for dFy/dpDOT)
        dgenFdpDOT = -1/2*((ztilde.^2+(Beta)^2).^(-1/2).*ztilde+1)*damping;
        dgenFdpDOT = diag(kron(dgenFdpDOT,[1 1 0]'));
        
        f = dgenFdpDOT;
    case 'dgenFdRL'
        %cable/rod inputs do not affect viscous XY force from floor
        dgenFdRL = zeros(size(X.p,1),size(X.RL,1),class(X.p));
        
        f = dgenFdRL;
    case 'dgenFdL'
        dgenFdL = zeros(size(X.p,1),size(X.L,1),class(X.p));

        f = dgenFdL;
end

end

