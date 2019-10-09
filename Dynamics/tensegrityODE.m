function dXdt = tensegrityODE(t,Xstacked,U,nominalFnc,hFcns)
%DPDT Summary of this function goes here
%   Detailed explanation goes here

dXdt = zeros(size(Xstacked));
num_nodalPos = (size(Xstacked,1)-numel(U.RLdot)-numel(U.Ldot))/2;
dXdt(1:num_nodalPos) = ...
    Xstacked(num_nodalPos+1:num_nodalPos*2); %copy over pDOT
X.p = Xstacked(1:num_nodalPos);
X.pDOT = Xstacked(num_nodalPos+1:num_nodalPos*2);
X.RL = Xstacked(num_nodalPos*2+1:num_nodalPos*2+numel(U.RLdot));
X.L = Xstacked(end-numel(U.Ldot)+1:end);

hVars = [];

%pDDOT
dXdt(num_nodalPos+1:num_nodalPos*2) = nominalFnc.pDDOT(X.p,X.pDOT,X.RL,X.L);

%RLdot
dXdt(num_nodalPos*2+1:num_nodalPos*2+size(U.RLdot,1)) = nominalFnc.RLdot(X,U,hVars);

%Ldot
dXdt(end-size(U.Ldot,1)+1:end) = nominalFnc.Ldot(X,U,hVars);

end

