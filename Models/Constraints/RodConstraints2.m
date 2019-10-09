function g = RodConstraints(X,U,hVars,omega,args,desFunc)
%RODCONSTRAINTS This function outputs the necessary matrices for
%calculating constraint forces in the rods. 
%   INPUTS:
%       X,U,hVars,omega,args,desFunc
%   OUTPUTS:
%       g - a function handle corresponding to a desired matrix depending
%       on the input argument 'desFunc'. This switch case allows for this
%       function to output all relevant matrices pertaining to this
%       specific constraint.

% 
% p = X.p;
% pDOT = X.pDOT;
% rodL = X.L;
% v = hVars.v;
% RhatRhat = hVars.RhatRhat;
% 
% numRods = length(rodL);
% 
% G = zeros(numRods,1,class(X.p));
% GDOT = zeros(numRods,1,class(X.p));
% dGdp = zeros(numRods,length(p),class(X.p));
% dGDOTdp = zeros(numRods,length(p),class(X.p));
% 
% for i = 1:numRods
%     RhatRhatbar = RhatRhat{i};
%     z = v{i};
%     G(i) = (1/2)*(z'*z) - rodL(i)^2;
%     GDOT(i) = p'*(RhatRhatbar)*pDOT;
%     dGdp(i,:) = p'*(RhatRhatbar);
%     dGDOTdp(i,:) = pDOT'*(RhatRhatbar);
% end


%select which function to output to anonymous function handle
switch desFunc
    case 'G'
        numRods = length(X.L);
        G = zeros(numRods,1,class(X.p));  
        for i = 1:numRods
            z = hVars.v{i};
            G(i) = (1/2)*(z'*z) - X.L(i)^2;
        end
        g = G;
    case 'GDOT'
        numRods = length(X.L);
        GDOT = zeros(numRods,1,class(X.p));
        for i = 1:numRods
            GDOT(i) = X.p'*(hVars.RhatRhat{i})*X.pDOT;
        end
        g = GDOT;
    case 'dGdp'      
        numRods = length(X.L);       
        dGdp = zeros(numRods,length(X.p),class(X.p));
        for i = 1:numRods
            dGdp(i,:) = X.p'*(hVars.RhatRhat{i});
        end
        g = dGdp;
    case 'dGDOTdp'  
        numRods = length(X.L);
        dGDOTdp = zeros(numRods,length(X.p),class(X.p));
        for i = 1:numRods
            dGDOTdp(i,:) = X.pDOT'*(hVars.RhatRhat{i});
        end
        g = dGDOTdp;
    case 'dGdpDOT'
        %todo
    case 'dGDOTdpDOT'
        %todo
    case 'dGdRL'
        %todo
    case 'dGDOTdRL'
        %todo
    case 'dGdL'
        %todo
    case 'dGDOTdL'
        %todo
end

end

