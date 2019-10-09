function g = StationaryConstraints(X,U,hVars,omega,args,desFunc)
%RODCONSTRAINTS This function outputs the necessary matrices for
%calculating constraint forces in the rods. 
%   INPUTS:
%       X,U,hVars,omega,args,desFunc
%   OUTPUTS:
%       g - a function handle corresponding to a desired matrix depending
%       on the input argument 'desFunc'. This switch case allows for this
%       function to output all relevant matrices pertaining to this
%       specific constraint.


p = X.p;
p0 = args.p0;
constrain = args.constrain;

constrainIdx = find(constrain);

numStationaryDir = numel(constrainIdx);

G = zeros(numStationaryDir,1);
GDOT = zeros(numStationaryDir,1);
dGdp = zeros(numStationaryDir,length(p));
dGDOTdp = zeros(numStationaryDir,length(p));

for i = 1:numStationaryDir
    G(i) = X.p(constrainIdx(i))-p0(constrainIdx(i));
    GDOT(i) = X.pDOT(constrainIdx(i));
    dGdp(i,constrainIdx(i)) = 1;
end

%select which function to output to anonymous function handle
switch desFunc
    case 'G'
        g = G;
    case 'GDOT'
        g = GDOT;
    case 'dGdp'
        g = dGdp;
    case 'dGDOTdp'
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

