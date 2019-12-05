function [cost] = velocityCost(X,U,args,runtimeArgs)
%ROLLINGCOST Summary of this function goes here
%   Detailed explanation goes here

%static costs
RL_diff_weight = args.RL_diff_weight;
RL_actuation_weight = args.RL_actuation_weight;
L_diff_weight = args.L_diff_weight;
velocity_reward = args.velocity_reward;
omega = args.omega;
stepDiscount = args.stepDiscount; 

cost = 0;

N = size(X.p,2)-1;

desired_direction = runtimeArgs.desired_direction;

RL_diff = X.RL-repmat(omega.X.RL0,1,size(X.RL,2));
RL_actuation = U.RLdot;
L_diff = X.L-repmat(omega.X.L0,1,size(X.L,2));

for k=1:N+1
    %penalize cable difference from neutral length
    if(RL_diff_weight~=0) %check if possible to avoid computing norm
        cost = cost + (stepDiscount^(k-1))*RL_diff_weight*norm(RL_diff(:,k),1);
    end
    
    %penalize rod difference from neutral length
    if(L_diff_weight~=0) %check if possible to avoid computing norm
        cost = cost + (stepDiscount^(k-1))*L_diff_weight*norm(L_diff(:,k),1);
    end
    
    %penalize cable actuation
    if(RL_actuation_weight~=0) %check if possible to avoid computing norm
        cost = cost + (stepDiscount^(k-1))*RL_actuation_weight*norm(RL_actuation(:,k),1);
    end
    
    %reward nodal velocity
    totalMass = sum(omega.M);
%     for node=1:(length(X.p(:,1))/3)
%         cost = cost - (velocity_reward)*(stepDiscount^(k-1))*omega.M(node)/totalMass*...
%             (desired_direction'*X.pDOT(node*3-2:node*3,k));
%     end

    
    cost = cost - (velocity_reward)*(stepDiscount^(k-1))*...
        kron(omega.M'/totalMass,desired_direction')*X.pDOT(:,k);

end

end

