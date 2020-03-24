function f = FloorForceVertical_XIncline(X,U,omega,args,desFunc)
%FLOORFORCE Summary of this function goes here
%   Detailed explanation goes here

% %args fields/values
% baseFloor = args.baseFloor; %floor height Z-value
% stiffness = args.stiffness; %floor stiffness
% damping = args.damping; %floor damping
% Beta = args.Beta; %parameter for softMax/softStep
% incline = args.incline; %incline in degrees in the +X direction
% 
% %smax (smooth approximation to max(x,0) and step)
% sMax = @(x,Beta) (sqrt(x.^2+Beta^2)+x)/2;
% sStep = @(x,Beta) 1/2*((x.^2+Beta^2).^(-1/2).*x+1);
% 
% %apply spring-damper floor forces
% %forces only apply to nodes in contact with ground
% %use differential smoothMax and smoothStep approximations
% genF = zeros(size(X.p),class(X.p));
% verticalLift = tand(incline)*X.p(1:3:end);
% genF(3:3:end) = ...
%     sMax((baseFloor + verticalLift - X.p(3:3:end))*stiffness,Beta) +...
%     sStep((baseFloor + verticalLift - X.p(3:3:end))*stiffness,Beta).*...
%     (-X.pDOT(3:3:end))*damping; 
% 
% x_tilde = baseFloor + verticalLift - X.p(3:3:end);
% x = x_tilde*stiffness;
% 
% dgenFdp = 1/2*((x.^2+(Beta)^2).^(-1/2).*x+1)*(-stiffness) +...
%     1/2*(-(x_tilde.^2+(Beta)^2).^(-3/2).*(x_tilde.^2)+...
%     (x_tilde.^2+(Beta)^2).^(-1/2))*damping.*(-X.pDOT(3:3:end))*(-1);
% dgenFdp = diag(kron(dgenFdp,[0 0 1]'));
% 
% dgenFdpDOT = -1/2*((x_tilde.^2+(Beta)^2).^(-1/2).*x_tilde+1)*damping;
% dgenFdpDOT = diag(kron(dgenFdpDOT,[0 0 1]'));
% 
% %cable/rod inputs do not affect vertical force from floor
% dgenFdRL = zeros(size(X.p,1),size(X.RL,1),class(X.p));
% dgenFdL = zeros(size(X.p,1),size(X.L,1),class(X.p));


%select which function to output for anonymous function handle
switch desFunc
    case 'genF'
        %args fields/values
        baseFloor = args.baseFloor; %floor height Z-value
        stiffness = args.stiffness; %floor stiffness
        damping = args.damping; %floor damping
        Beta = args.Beta; %parameter for softMax/softStep
        incline = args.incline; %incline in degrees in the +X direction
        
        %smax (smooth approximation to max(x,0) and step)
        sMax = @(x,Beta) (sqrt(x.^2+Beta^2)+x)/2;
        sStep = @(x,Beta) 1/2*((x.^2+Beta^2).^(-1/2).*x+1);
        
        %apply spring-damper floor forces
        %forces only apply to nodes in contact with ground
        %use differential smoothMax and smoothStep approximations
        genF = zeros(size(omega.X0,1),1,class(X.p));
        verticalLift = tand(incline)*X.p(1:3:end);
        %         genF(3:3:end) = ...
        %             sMax((baseFloor + verticalLift - X.p(3:3:end))*stiffness,Beta) +...
        %             sStep((baseFloor + verticalLift - X.p(3:3:end))*stiffness,Beta).*...
        %             (-X.pDOT(3:3:end))*damping;
        genF(3:3:end) = ...
            sMax((baseFloor + verticalLift - X.p(3:3:end))*stiffness,Beta) +...
            sStep((baseFloor + verticalLift - X.p(3:3:end)),Beta).*...
            (-X.pDOT(3:3:end))*damping;
        f = genF;
    case 'dgenFdp'
        %args fields/values
        baseFloor = args. baseFloor; %floor height Z-value
        stiffness = args.stiffness; %floor stiffness
        damping = args.damping; %floor damping
        Beta = args.Beta; %parameter for softMax/softStep
        incline = args.incline; %incline in degrees in the +X direction
                
        %apply spring-damper floor forces
        %forces only apply to nodes in contact with ground
        %use differential smoothMax and smoothStep approximations
        verticalLift = tand(incline)*X.p(1:3:end);

        x_tilde = baseFloor + verticalLift - X.p(3:3:end);
        x = x_tilde*stiffness;
        
        dgenFdp = 1/2*((x.^2+(Beta)^2).^(-1/2).*x+1)*(-stiffness) +...
            1/2*(-(x_tilde.^2+(Beta)^2).^(-3/2).*(x_tilde.^2)+...
            (x_tilde.^2+(Beta)^2).^(-1/2))*damping.*(-X.pDOT(3:3:end))*(-1);
        dgenFdp = diag(kron(dgenFdp,[0 0 1]'));

        f = dgenFdp;
    case 'dgenFdpDOT'
        %args fields/values
        baseFloor = args.baseFloor; %floor height Z-value
        damping = args.damping; %floor damping
        Beta = args.Beta; %parameter for softMax/softStep
        incline = args.incline; %incline in degrees in the +X direction
        
        %apply spring-damper floor forces
        %forces only apply to nodes in contact with ground
        %use differential smoothMax and smoothStep approximations
        verticalLift = tand(incline)*X.p(1:3:end);
        
        x_tilde = baseFloor + verticalLift - X.p(3:3:end);
                
        dgenFdpDOT = -1/2*((x_tilde.^2+(Beta)^2).^(-1/2).*x_tilde+1)*damping;
        dgenFdpDOT = diag(kron(dgenFdpDOT,[0 0 1]'));
        
        f = dgenFdpDOT;
    case 'dgenFdRL'
        %cable/rod inputs do not affect vertical force from floor
        dgenFdRL = zeros(size(omega.X0,1),size(omega.C,1),class(X.p));
        f = dgenFdRL;
    case 'dgenFdL'
        dgenFdL = zeros(size(omega.X0,1),size(omega.R,1),class(X.p));
        f = dgenFdL;
end

end

