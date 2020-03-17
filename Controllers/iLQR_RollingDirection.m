classdef iLQR_RollingDirection < handle
    %QP_MPC_ROLLINGDIRECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        targetDestination
        horizon
        omega
        cableLinearVelocity
        cableMinLength
        cableMaxLength
        actuationMode %1-cables only,2-rods only,3-both
        dT
        Q
        R
        P
        nX_p
        nX_pDOT
        nX_RL
        nX_L
        nX
        nU_RLdot
        cableConstraintMatrix
        rodConstraintMatrix
        initial_step_guess
        
        %iLQR
        uGuess
        xTraj
        Amats
        Bmats
        Qmats
        Rmats
        Nmats
        Kgains
        
    end
    
    methods
        function obj = iLQR_RollingDirection(...
                X,omega,simTimeStep,horizon)
            %QP_MPC_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
            %controller timestep and horizon
            obj.dT = simTimeStep;
            obj.horizon = horizon;
            
            % ODE solver initial step guess
            obj.initial_step_guess = obj.dT/20;
            
            %initialize origin as default target of [0,0]
            obj.targetDestination = [0 0]';
            
            %actuator limits
            obj.omega = omega;
            obj.cableLinearVelocity = omega.cables.linear_velocity;
            obj.cableMinLength = omega.cables.minLength;
            obj.cableMaxLength = omega.cables.maxLength;
            
            obj.nX_p = numel(X.p);
            obj.nX_pDOT = numel(X.pDOT);
            obj.nX_RL = numel(X.RL);
            obj.nX_L = numel(X.L);
            obj.nX = obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L;
            obj.nU_RLdot = size(omega.cableConstraintMatrix,2); %constrained cable inputs
            
            obj.cableConstraintMatrix = omega.cableConstraintMatrix;
            obj.rodConstraintMatrix = omega.rodConstraintMatrix;
            
            %LQR cost function
            %state penalty
            obj.Q = zeros(obj.nX+1);
            
            %cable deviation penalty (careful, on if the state is the
            %deviation from RL0 or actual state itself)
            cableDeviationPenalty = 5e1;%5e1;
            obj.Q(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL,...
                obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
                cableDeviationPenalty*eye(obj.nX_RL); %1e0 works for 6-bar,
            
            %velocity input penalty
            cableVel_cost = 1e-1;
            obj.R = cableVel_cost*blkdiag(...
                omega.cableConstraintMatrix'*eye(obj.nX_RL)*omega.cableConstraintMatrix,...
                eye(obj.nX_L)); %~1e-2 works for dT=1e-3; 5e1 works for dT=5e-3
            obj.R(obj.nU_RLdot+1:end,obj.nU_RLdot+1:end) =...
                1e6*eye(obj.nX_L); %penalize rod actuation heavily
            
            %cell array of cost-to-go and feedback matrices
            obj.P = cell(horizon+1,1);
            
            %Jacobians
            obj.Amats = zeros(obj.nX+1,obj.nX+1,horizon+1);
            obj.Bmats = zeros(obj.nX+1,obj.nU_RLdot+obj.nX_L,horizon+1);
            obj.Qmats = zeros(obj.nX+1,obj.nX+1,horizon+1);
            obj.Rmats = zeros(obj.nU_RLdot+obj.nX_L,obj.nU_RLdot+obj.nX_L,horizon+1);
            obj.Nmats = zeros(obj.nX+1,obj.nU_RLdot+obj.nX_L,horizon+1);
            obj.Kgains = zeros(obj.nU_RLdot+obj.nX_L,obj.nX+1,horizon+1);
            
            %initialize iLQR input variables
            obj.uGuess = 1e-6*randn(size(omega.cableConstraintMatrix,2)+obj.nX_L,horizon);
            obj.uGuess(end-obj.nX_L+1:end,:)=zeros(obj.nX_L,size(obj.uGuess,2));
        end
        
        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                costOutput,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,nominalFnc,jacobianFcns,hFcns,costFcnHandle,...
                debugFcns,openLoopFlag)
            
            N = obj.horizon;
                        
            %calculate desired rolling direction according to target goal
            totalMass = sum(obj.omega.M);
            xCOM = sum(obj.omega.M.*Xhat.p(1:3:end))/totalMass;
            yCOM = sum(obj.omega.M.*Xhat.p(2:3:end))/totalMass;
            zCOM = sum(obj.omega.M.*Xhat.p(3:3:end))/totalMass;
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            
            zWeight = sum(obj.omega.M.*obj.omega.X.p0(3:3:end))/totalMass-zCOM;
%             zWeight = sum(obj.omega.M(end-1:end).*obj.omega.X.p0(end-3:3:end))/...
%                 sum(obj.omega.M(end-1:end))-...
%                 sum(obj.omega.M(end-1:end).*Xhat.p(end-3:3:end))/sum(obj.omega.M(end-1:end));

            desiredDirection = [desiredDirection;10*zWeight];
            desiredDirection = desiredDirection/norm(desiredDirection);
            %             desiredDirection(3) = desiredDirection(3)+0.5;
            %             desiredDirection = desiredDirection/norm(desiredDirection);
            
            %linear reward (velocity)
            velReward = 5e0;
            obj.Q(end,obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT) =...
                -velReward/2*desiredDirection(1)*obj.omega.M/totalMass;%.*(Xhat.p(3:3:end)-mean(Xhat.p(3:3:end)));
            obj.Q(end,obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT) =...
                -velReward/2*desiredDirection(2)*obj.omega.M/totalMass;%.*(Xhat.p(3:3:end)-mean(Xhat.p(3:3:end)));
            obj.Q(end,obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT) =...
                -velReward/2*desiredDirection(3)*obj.omega.M/totalMass;
            obj.Q(obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT,end) =...
                -velReward/2*desiredDirection(1)*obj.omega.M/totalMass;%.*(Xhat.p(3:3:end)-mean(Xhat.p(3:3:end)));
            obj.Q(obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT,end) =...
                -velReward/2*desiredDirection(2)*obj.omega.M/totalMass;%.*(Xhat.p(3:3:end)-mean(Xhat.p(3:3:end)));
            obj.Q(obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT,end) =...
                -velReward/2*desiredDirection(3)*obj.omega.M/totalMass;
                        
            
            %store important values for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            %reference state
            Xref = [zeros(obj.nX_p,1);
                zeros(obj.nX_pDOT,1);
                obj.omega.X.RL0; %neutral cable lengths
                zeros(obj.nX_L,1)];
            %concatenated state
            Xcombined = [Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L];
            [U_OL,X_OL,costOutput,hVars] = obj.iLQR(Xcombined,obj.Q,obj.R,Xref);
            
            %remap U input vector (potentially smaller due to constraints)
            %to full cable/rod actuation
            U_desired.RLdot = obj.cableConstraintMatrix*...
                U_OL(1:size(obj.cableConstraintMatrix,2),1);
            U_desired.Ldot = obj.rodConstraintMatrix*...
                U_OL(size(obj.cableConstraintMatrix,2)+1:end,1);
            
            if(any(abs(U_desired.RLdot)>obj.cableLinearVelocity))
                disp('WARNING: Desired Cable velocity exceeds motor capability')
                
                %                 disp('Clipping cable velocity input')
                %                 U_desired.RLdot = sign(U_desired.RLdot).*...
                %                     min(abs(U_desired.RLdot),...
                %                     obj.cableLinearVelocity);
                %
                %                 disp('Normalizing Velocity Input according to Motor Limits~~~~~~~~~~~~~~~~')
                %                 U_desired.RLdot = U_desired.RLdot/max(abs(U_desired.RLdot)).*obj.cableLinearVelocity*1.5;
            end
            
            OL_states.p = zeros(obj.nX_p,N+1);
            OL_states.pDOT = zeros(obj.nX_pDOT,N+1);
            OL_states.RL = zeros(obj.nX_RL,N+1);
            OL_states.L = zeros(obj.nX_L,N+1);
            OL_inputs.RLdot = zeros(obj.nX_RL,N+1);
            OL_inputs.Ldot = zeros(obj.nX_L,N+1);
            
            %parse into open-loop trajectory expected format
            if(openLoopFlag)
                for i=1:N+1
                    x_vec = X_OL(:,i);
                    %parse vector into struct
                    idx=0;
                    OL_states.p(:,i) = x_vec(idx+1:idx+obj.nX_p);
                    idx=idx+obj.nX_p;
                    OL_states.pDOT(:,i) = x_vec(idx+1:idx+obj.nX_pDOT);
                    idx=idx+obj.nX_pDOT;
                    OL_states.RL(:,i) = x_vec(idx+1:idx+obj.nX_RL);
                    idx=idx+obj.nX_RL;
                    OL_states.L(:,i) = x_vec(idx+1:idx+obj.nX_L);
                    if(i<=N)
                        u_vec = U_OL(:,i);
                        OL_inputs.RLdot(:,i) = obj.cableConstraintMatrix*...
                            u_vec(1:size(obj.cableConstraintMatrix,2));
                        OL_inputs.Ldot(:,i) = obj.rodConstraintMatrix*...
                            u_vec(size(obj.cableConstraintMatrix,2)+1:end);
                    end
                end
            else
                OL_states = [];
                OL_inputs = [];
            end
        end
        
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
        end
        
        
        function [U_optimal,OL_states,currCost,hVars] = iLQR(obj,X0,Q,R,Xref)
            N = obj.horizon;
            
            converged = 0; %initialize convergence flag to 'false'
            obj.xTraj = zeros(obj.nX,N+1);
            obj.xTraj(:,1) = X0; %initial state
            %             obj.uGuess(:,1:end-1) = obj.uGuess(:,2:end); %use last inputs as initial guess
            %             obj.uGuess = obj.uGuess+1e-3*randn(size(obj.uGuess));
            %             obj.uGuess(:,1) = obj.uGuess(:,1)+1e-3*randn(size(obj.uGuess,1),1);
            %             obj.uGuess(end-obj.nX_L+1:end,:)=zeros(obj.nX_L,size(obj.uGuess,2));
%             obj.uGuess = repmat(obj.uGuess(:,1),1,N);
            
            %forward simulate to get initial state trajectory
            for k = 1:N
                obj.xTraj(:,k+1) = stepForward(obj,obj.uGuess(:,k),obj.xTraj(:,k));
            end
            
%             obj.xTraj = stepForward(obj,obj.uGuess,obj.xTraj(:,1));
            
            %obj.xTraj0 = obj.xTraj;
            bestCost = 1e8;
            Alpha = 1;
            while(~converged)
                Alpha = Alpha*5;%  %increase stepsize from last linesearch
                uTempTop = obj.uGuess;
%                 obj.uGuess = obj.uGuess+1e-1*randn(size(obj.uGuess));
                for t = N:-1:1 %for each timestep horizon to 1
                    %1 is last, so that hVars at t=t0 can be passed as output to dynamics
                    
                    %parse x at timestep t
                    idx=0;
                    Xhat.p = obj.xTraj(idx+1:obj.nX_p,t);
                    idx = idx + obj.nX_p;
                    Xhat.pDOT = obj.xTraj(idx+1:idx+obj.nX_pDOT,t);
                    idx = idx + obj.nX_pDOT;
                    Xhat.RL = obj.xTraj(idx+1:idx+obj.nX_RL,t);
                    idx = idx + obj.nX_RL;
                    Xhat.L = obj.xTraj(idx+1:end,t);
                    %parse u at timestep t
                    Uhat.RLdot = obj.uGuess(1:obj.nU_RLdot,t);
                    Uhat.Ldot = obj.uGuess(obj.nU_RLdot+1:end,t);
                    
                    %helper variables
                    for i = 1:length(obj.omega.hFcns.z)
                        hVars.z{i} = obj.omega.hFcns.z{i}(Xhat);
                    end
                    for i = 1:length(obj.omega.hFcns.v)
                        hVars.v{i} = obj.omega.hFcns.v{i}(Xhat);
                    end
                    hVars.RhatRhat = obj.omega.hFcns.RhatRhat;
                    hVars.Chat = obj.omega.hFcns.Chat;
                    hVars.J = obj.omega.hFcns.J(Xhat,Uhat,hVars);
                    
%                     %Nominal values for linearization
%                     pDDOT_bar = obj.omega.nominalFcn.pDDOT(Xhat.p,Xhat.pDOT,Xhat.RL,Xhat.L);
%                     RLdot_bar =  obj.omega.nominalFcn.RLdot(Xhat,Uhat,hVars);
%                     Ldot_bar =  obj.omega.nominalFcn.Ldot(Xhat,Uhat,hVars);
                    
                    %dynamic linearization recalculated each controller iteration
                    dpDDOTdX = [obj.omega.jacobianFcns.dpDDOTdp(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dpDDOTdpDOT(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dpDDOTdRL(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dpDDOTdL(Xhat,Uhat,hVars)];
                    dRLdX = [obj.omega.jacobianFcns.dRLdp(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dRLdpDOT(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dRLdRL(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dRLdL(Xhat,Uhat,hVars)];
                    dRLdU= [obj.omega.jacobianFcns.dRLdRLdot(Xhat,Uhat,hVars)*obj.cableConstraintMatrix,...
                        zeros(size(obj.omega.C,1),size(obj.rodConstraintMatrix,2))];
                    dLdX = [obj.omega.jacobianFcns.dLdp(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dLdpDOT(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dLdRL(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dLdL(Xhat,Uhat,hVars)];
                    dLdU= [zeros(size(obj.omega.R,1),size(obj.cableConstraintMatrix,2)),...
                        obj.omega.jacobianFcns.dLdLdot(Xhat,Uhat,hVars)*obj.rodConstraintMatrix];
                    
                    
                    Abar = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L);
                    Abar(1:obj.nX_p,obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = ...
                        eye(obj.nX_pDOT);
                    Abar(obj.nX_p+1:end,:) = [dpDDOTdX;dRLdX;dLdX];
                    Bbar = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L,...
                        size(obj.cableConstraintMatrix,2)+size(obj.rodConstraintMatrix,2));
                    Bbar(end-(obj.nX_RL+obj.nX_L)+1:end,:) = ...
                        [dRLdU;dLdU];
                    %                     xDOTBar = [Xhat.pDOT;pDDOT_bar;obj.cableConstraintMatrix*RLdot_bar;Ldot_bar]; %nominal state velocities
                    %                     zHat = obj.dT*(xDOTBar-Abar*[Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L]...
                    %                         -[zeros(obj.nX_p+obj.nX_pDOT,obj.nX_RL+obj.nX_L);... %block of zeros in Bbar
                    %                         eye(obj.nX_RL),zeros(obj.nX_RL,obj.nX_L);... %Bbar: dRLdU
                    %                         zeros(obj.nX_L,obj.nX_RL),eye(obj.nX_L)]*... %Bbar: dLdU
                    %                         [Uhat.RLdot;Uhat.Ldot]);
                    zHat = zeros(obj.nX,1);
                    
                    %form large augmented matrices
                    obj.Amats(:,:,t) = [eye(obj.nX)+obj.dT*Abar,...
                        zHat; zeros(1,obj.nX),1];
                    obj.Bmats(:,:,t) = [obj.dT*Bbar;zeros(1,size(Bbar,2))];
                    decay = 0.95;
                    obj.Qmats(:,:,t) = Q;
                    obj.Qmats(:,end,t) = decay^t*obj.Qmats(:,end,t);
                    obj.Qmats(end,:,t) = decay^t*obj.Qmats(end,:,t);
                    obj.Qmats(:,end,t) = obj.Qmats(:,end,t) + ...
                        decay^t*Q*[(obj.xTraj(:,t)-Xref);1];
                    obj.Qmats(end,:,t) = obj.Qmats(end,:,t) + ...
                        (decay^t*Q*[(obj.xTraj(:,t)-Xref);1])';
                    obj.Qmats(end,end,t) = obj.Qmats(end,end,t) + ...
                        [(obj.xTraj(:,t)-Xref);1]'*decay^t*Q*[(obj.xTraj(:,t)-Xref);1]+...
                        obj.uGuess(:,t)'*decay^t*R*obj.uGuess(:,t);
                    obj.Rmats(:,:,t) = R*decay^t;
                    obj.Nmats(:,:,t) = [zeros(obj.nX,obj.nU_RLdot+obj.nX_L);
                        obj.uGuess(:,t)'*R*decay^t];
                    
                    %make Q positive definite
                    obj.Qmats(:,:,t) = obj.Qmats(:,:,t)+(abs(min(real(eig(obj.Qmats(:,:,t)))))+1e-6)*eye(size(obj.Qmats(:,:,t)));
                end
                
                %terminal cost
                SF = 1;
                obj.Qmats(:,:,N+1) = SF*Q;
                obj.Qmats(:,end,N+1) = obj.Qmats(:,end,N+1) + ...
                    SF*Q*[(obj.xTraj(:,N+1)-Xref);1];
                obj.Qmats(end,:,N+1) = obj.Qmats(end,:,N+1) + ...
                    SF*(Q*[(obj.xTraj(:,N+1)-Xref);1])';
                obj.Qmats(end,end,N+1) = obj.Qmats(end,end,N+1) + ...
                    [(obj.xTraj(:,N+1)-Xref);1]'*SF*Q*[(obj.xTraj(:,N+1)-Xref);1];
                
                %make Q positive definite
                obj.Qmats(:,:,N+1) = obj.Qmats(:,:,N+1)+(abs(min(real(eig(obj.Qmats(:,:,N+1)))))+1e-6)*eye(size(obj.Qmats(:,:,N+1)));
                
                %calculate LQR feedback gains
                %backward pass
                obj.P{N+1} = obj.Qmats(:,:,N+1); %terminal cost
                for k = N:-1:1
                    %stage cost
                    obj.P{k} = obj.Amats(:,:,k)'*obj.P{k+1}*obj.Amats(:,:,k)-...
                        (obj.Amats(:,:,k)'*obj.P{k+1}*obj.Bmats(:,:,k)+obj.Nmats(:,:,k))/...
                        (obj.Rmats(:,:,k)+...
                        obj.Bmats(:,:,k)'*obj.P{k+1}*obj.Bmats(:,:,k))*...
                        (obj.Bmats(:,:,k)'*obj.P{k+1}*obj.Amats(:,:,k)+obj.Nmats(:,:,k)')+...
                        obj.Qmats(:,:,k);
                    %optimal LQR feedback gain
                    obj.Kgains(:,:,k) = (obj.Rmats(:,:,k)+...
                        obj.Bmats(:,:,k)'*obj.P{k+1}*obj.Bmats(:,:,k))\...
                        (obj.Bmats(:,:,k)'*obj.P{k+1}*obj.Amats(:,:,k)+...
                        obj.Nmats(:,:,k)');
                    obj.Kgains(end-obj.nX_L+1:end,:,k) = zeros(obj.nX_L,obj.nX+1);
                end
                
                %calculate initial trajectory cost
                cost = calculateCost(obj,Xref,Q,R);
                
                %Iterate Variational Calculus
                while(true) %do continuously
                    deltaX = zeros(size(obj.xTraj)); %X0 is fixed
                    for t = 1:N
                        %check input constraint violation
                        checkConstraints = 1;
                        
                        velExceeded=[];
                        j=[];
                        obj.uGuess(:,t) = obj.uGuess(:,t)-Alpha*obj.Kgains(:,end,t);
                        if(checkConstraints)
                            %handle paired/similar cable constraints
                            velExceeded = abs(obj.cableConstraintMatrix*obj.uGuess(1:obj.nU_RLdot,t))>obj.omega.cables.linear_velocity;
                            cables = find(velExceeded);
                            [~,j] = find(obj.cableConstraintMatrix(velExceeded,:));
                            for k = 1:numel(j)
                                obj.uGuess(j(k),t) = sign(obj.uGuess(j(k),t))*...
                                    min(abs(obj.uGuess(j(k),t)),...
                                    obj.omega.cables.linear_velocity(cables(k)));
                            end
                        end
                        Khat = obj.Kgains(:,1:end-1,t);
                        Khat(j,:) = zeros(numel(j),obj.nX);
                        Khat(end-obj.nX_L+1:end,:) = zeros(obj.nX_L,obj.nX);
                        obj.uGuess(:,t) = obj.uGuess(:,t)-Khat*deltaX(:,t);
                        if(checkConstraints)
                            %handle paired/similar cable constraints
                            velExceeded = abs(obj.cableConstraintMatrix*obj.uGuess(1:obj.nU_RLdot,t))>obj.omega.cables.linear_velocity;
                            cables = find(velExceeded);
                            [~,j] = find(obj.cableConstraintMatrix(velExceeded,:));
                            for k = 1:numel(j)
                                obj.uGuess(j(k),t) = sign(obj.uGuess(j(k),t))*...
                                    min(abs(obj.uGuess(j(k),t)),...
                                    obj.omega.cables.linear_velocity(cables(k)));
                            end
                        end
                        temp = obj.xTraj(:,t+1);
                        obj.xTraj(:,t+1) = stepForward(obj,obj.uGuess(:,t),obj.xTraj(:,t));
                        deltaX(:,t+1) = obj.xTraj(:,t+1)-temp;
                        %                         deltaX = obj.Amats(1:end-1,1:end-1,t)*deltaX+...
                        %                             obj.Bmats(1:end-1,:,t)*(obj.uGuess(:,t)-uTemp(:,t));
                    end                    
%                     %forward simulate to get new state trajectory
%                     obj.xTraj(:,1) = X0;
%                     for k = 1:N
%                         obj.xTraj(:,k+1) = stepForward(obj,obj.uGuess(:,k),obj.xTraj(:,k));
%                     end

%                     obj.xTraj = stepForward(obj,obj.uGuess,obj.xTraj(:,1));
                    
                    %calculate new trajectory cost
                    currCost = calculateCost(obj,Xref,Q,R);
%                     z = [obj.xTraj(:,N+1)-Xref;1];
%                     currCost = z'*Q*z;
%                     for t = 1:N
%                         z = [obj.xTraj(:,t)-Xref;1];
%                         currCost = currCost + z'*Q*z +...
%                             obj.uGuess(:,t)'*R*obj.uGuess(:,t);
%                     end
                    
                    if(currCost<bestCost)
                        bestCost = currCost;
                    end
                    
                    if(1) %toggle on/off for debugging
                        costDiff = cost-currCost
                        figure(2)
                        cla
                        plot(obj.uGuess')
                        hold on
                        plot(uTempTop','--')
                        drawnow()
                        title(['Best Cost = ',num2str(bestCost,'%10.5e'),...
                            '; Current Cost = ',num2str(currCost,'%10.5e'),...
                            '; Alpha = ',num2str(Alpha)])
                    end
                                        
                    %                     if(currCost<=cost)
                    %                         %disp('Finished Line Search')
                    %                         %                         Alpha
                    %                         break
                    if(checkArmijo(obj,deltaX,obj.uGuess-uTempTop,Xref,Q,R,cost,currCost,0.5))
                        %disp('Finished Line Search')
                        %                         Alpha
                        break
                    elseif(Alpha<1e-2)
                        disp('Alpha small, moving on')
                        diff = cost-currCost
                        converged = 1;
                        break
                    else
                        %reset inputs, reduce line search stepsize
                        obj.uGuess = uTempTop;
                        %                         obj.xTraj = xTemp;
                        Alpha = Alpha/2;
                        %check if first two inputs have converged, if so - exit
%                         if(max(max(abs(obj.uGuess(:,1)-uTempTop(:,1))))<1e-4)
%                             converged = 1;
%                             disp('First few inputs converged, moving on...')
%                             break
%                         end
                    end
                end
                
                if(max(max(abs(obj.uGuess(:,1:2)-uTempTop(:,1:2))))<1e-4)
                    converged = 1;
                    disp('First few inputs converged, moving on...')
                end
                
                %                 %check if converged
                %                 plot(obj.uGuess')
                %                 plot(obj.xTraj(73:96,:)')
                %
                %                 uDiff = norm(obj.uGuess-uTempTop)
                                
            end
            
            U_optimal = obj.uGuess;
            OL_states = obj.xTraj;
        end
        
        
        function XOUT = stepForward(obj,u_vec,x_vec)
            %parse input vectors into structs
            idx=0;
            X.p = x_vec(idx+1:idx+obj.nX_p,:);
            idx=idx+obj.nX_p;
            X.pDOT = x_vec(idx+1:idx+obj.nX_pDOT,:);
            idx=idx+obj.nX_pDOT;
            X.RL = x_vec(idx+1:idx+obj.nX_RL,:);
            idx=idx+obj.nX_RL;
            X.L = x_vec(idx+1:idx+obj.nX_L,:);
            U.RLdot = obj.cableConstraintMatrix*...
                u_vec(1:size(obj.cableConstraintMatrix,2),:);
            U.Ldot = obj.rodConstraintMatrix*...
                u_vec(size(obj.cableConstraintMatrix,2)+1:end,:);
            
            N = size(U.RLdot,2);
            
            %cable update, subject to speed limits
            Uinput.RLdot = sign(U.RLdot).*...
                min(abs(U.RLdot),...
                obj.omega.cables.linear_velocity);
            %rod update, subject to speed limits
            Uinput.Ldot = sign(U.Ldot).*...
                min(abs(U.Ldot),...
                obj.omega.rods.linear_velocity);
            
            %cable min/max length limits
            Uinput.RLdot = max(Uinput.RLdot,...
                (obj.omega.cables.minLength-X.RL)/obj.dT);
            Uinput.RLdot = min(Uinput.RLdot,...
                (obj.omega.cables.maxLength-X.RL)/obj.dT);
            %rod min/max length limits
            Uinput.Ldot = max(Uinput.Ldot,...
                (obj.omega.rods.minLength-X.L)/obj.dT);
            Uinput.Ldot = min(Uinput.Ldot,...
                (obj.omega.rods.maxLength-X.L)/obj.dT);

            %forward simulate dynamics with ODE solver
            XIN = [X.p;X.pDOT;X.RL;X.L];
            options = odeset('RelTol',1e-2,'AbsTol',1e-2,...
                'InitialStep',obj.initial_step_guess);
            [t,XOUT] = ode113(@(t,Xstacked)...
                tensegrityODE(t,Xstacked, Uinput,obj.omega.nominalFcn,...
                obj.omega.hFcns), [0,obj.dT], XIN, options);
            obj.initial_step_guess = t(end)-t(end-4);
            XOUT = XOUT(end,:)';
        end
        
        function cost = calculateCost(obj,Xref,Q,R)
            N = obj.horizon;
            %calculate trajectory cost
            z = [obj.xTraj(:,N+1)-Xref;1];
            cost = z'*obj.Q*z;
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1];
                cost = cost + z'*Q*z +...
                    obj.uGuess(:,t)'*R*obj.uGuess(:,t);
            end
        end
        
        function satisfied = checkArmijo(obj,xDelta,uDelta,Xref,Q,R,prevCost,currCost,ratio)
            N = obj.horizon;
            %calculate trajectory cost
            J = prevCost;
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1]; %xTraj,uguess updated by this point
                J = J + ratio*(2*(z-[xDelta(:,t);0])'*Q*[xDelta(:,t);0] +...
                    2*(obj.uGuess(:,t)-uDelta(:,t))'*R*uDelta(:,t));
            end
            
            satisfied = currCost<J;
        end
        
    end
end

