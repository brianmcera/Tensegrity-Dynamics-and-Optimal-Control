classdef iLQR_RollingDirection_2 < handle
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
        nU_Ldot
        nU
        cableConstraintMatrix
        rodConstraintMatrix
        initial_step_guess
        
        %iLQR
        uGuess
        uDeltaGuess
        xTraj
        Amats
        Bmats
        QPrimemats
        RPrimemats
        NPrimemats
        Kgains
        inputChangePenalty;
        
        discount1
        discount2
        
        
    end
    
    methods
        function obj = iLQR_RollingDirection_2(...
                X,omega,simTimeStep,horizon)
            %QP_MPC_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
            %controller timestep and horizon
            obj.dT = simTimeStep;
            obj.horizon = horizon;
            
            % ode solver initial step
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
            obj.nU_Ldot = size(omega.rodConstraintMatrix,2); %constrained rod inputs
            obj.nU = obj.nU_RLdot + obj.nU_Ldot;
            
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
            cableVel_cost = 5e-4;
            obj.R = cableVel_cost*blkdiag(...
                omega.cableConstraintMatrix'*eye(obj.nX_RL)*...
                omega.cableConstraintMatrix,...
                eye(obj.nX_L)); 
            % note above: ~1e-2 works for dT=1e-3; 5e1 works for dT=5e-3
            rodVel_cost = 1e6;
            obj.R(obj.nU_RLdot+1:end,obj.nU_RLdot+1:end) =...
                rodVel_cost*eye(obj.nX_L); %penalize rod actuation heavily
            
            %cell array of cost-to-go and feedback matrices
            obj.P = cell(horizon+1,1);
            
            %Jacobians
            obj.Amats = zeros(obj.nX+1,obj.nX+1,horizon+1);
            obj.Bmats = zeros(obj.nX+1,obj.nU_RLdot+obj.nX_L,horizon+1);
            obj.QPrimemats = zeros(obj.nX+obj.nU+2,obj.nX+obj.nU+2,horizon+1);
            obj.RPrimemats = zeros(obj.nU_RLdot+obj.nU_Ldot,obj.nU_RLdot+obj.nU_Ldot,horizon+1);
            obj.NPrimemats = zeros(obj.nX+obj.nU+2,obj.nU_RLdot+obj.nX_L,horizon+1);
            obj.Kgains = zeros(obj.nU_RLdot+obj.nX_L,obj.nX+obj.nU+2,horizon+1);
            
            %initialize iLQR input variables
            obj.uGuess = zeros(size(omega.cableConstraintMatrix,2)+obj.nX_L,horizon+1);
            obj.uDeltaGuess = 0*randn(size(omega.cableConstraintMatrix,2)+obj.nU_Ldot,horizon);
            obj.uDeltaGuess(end-obj.nX_L+1:end,:)=zeros(obj.nX_L,size(obj.uDeltaGuess,2));
            
            % input penalty
            obj.inputChangePenalty = 5e0*eye(obj.nU);
            
            % discount factors
            obj.discount1 = 0.95; % discount factor on state penalty
            obj.discount2 = 0.95; % discount factor on input penalty
        end
        
        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                costOutput,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,nominalFnc,~,hFcns,costFcnHandle,...
                debugFcns,openLoopFlag)
            
            N = obj.horizon;
                        
            % calculate desired rolling direction according to target goal
            totalMass = sum(obj.omega.M);
            xCOM = sum(obj.omega.M.*Xhat.p(1:3:end))/totalMass;
            yCOM = sum(obj.omega.M.*Xhat.p(2:3:end))/totalMass;
            zCOM = sum(obj.omega.M.*Xhat.p(3:3:end))/totalMass;
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            
            zWeight = 1.1*sum(obj.omega.M.*obj.omega.X.p0(3:3:end))/...
                totalMass-zCOM; % 10% higher than neutral position
            zSF = 10; % how much to weigh z-deviation
            desiredDirection = [desiredDirection;zSF*zWeight];
            desiredDirection = desiredDirection/norm(desiredDirection);
            
            %linear penalty (velocity)
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
                U_OL(1:size(obj.cableConstraintMatrix,2),2);
            U_desired.Ldot = obj.rodConstraintMatrix*...
                U_OL(size(obj.cableConstraintMatrix,2)+1:end,2);
            %             U_desired.RLdot = obj.cableConstraintMatrix*...
            %                 mean(U_OL(1:size(obj.cableConstraintMatrix,2),:),2);
            %             U_desired.Ldot = obj.rodConstraintMatrix*...
            %                 mean(U_OL(size(obj.cableConstraintMatrix,2)+1:end,:),2);
            
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
        
        
        function [U_optimal,OL_states,currCost,hVars] = iLQR(...
                obj,X0,Q,R,Xref)
            N = obj.horizon;
            
            %advance initial input
%             obj.uGuess(:,1:end-1) = obj.uGuess(:,2:end);
            obj.uGuess = repmat(obj.uGuess(:,2),1,N+1);
            obj.uDeltaGuess = diff(obj.uGuess,1,2);
            
            converged = 0; %initialize convergence flag to 'false'
            obj.xTraj = zeros(obj.nX,N+1);
            obj.xTraj(:,1) = X0; %initial state
            
            %forward simulate to get initial state trajectory
            for k = 1:N
                obj.xTraj(:,k+1) = stepForward(obj,obj.uGuess(:,k),...
                    obj.xTraj(:,k));
            end
            
            bestCost = 1e8;
            Alpha = 0.5;
            while(~converged)
                Alpha = min(1,Alpha*1.5);%  %increase stepsize from last linesearch
%                 Alpha = 1;
                
                for t = N:-1:1 %for each timestep horizon to 1
                    %1 is last, so that hVars at t=t0 can be passed 
                    %as output to dynamics
                    
                    % parse x at timestep t
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
                    
                    % helper variables
                    for i = 1:length(obj.omega.hFcns.z)
                        hVars.z{i} = obj.omega.hFcns.z{i}(Xhat);
                    end
                    for i = 1:length(obj.omega.hFcns.v)
                        hVars.v{i} = obj.omega.hFcns.v{i}(Xhat);
                    end
                    hVars.RhatRhat = obj.omega.hFcns.RhatRhat;
                    hVars.Chat = obj.omega.hFcns.Chat;
                    hVars.J = obj.omega.hFcns.J(Xhat,Uhat,hVars);
                    
                    % dynamic linearization recalculated each iteration
%                     dpDDOTdX = [obj.omega.jacobianFcns.dpDDOTdp(Xhat,Uhat,hVars),...
%                         obj.omega.jacobianFcns.dpDDOTdpDOT(Xhat,Uhat,hVars),...
%                         obj.omega.jacobianFcns.dpDDOTdRL(Xhat,Uhat,hVars),...
%                         obj.omega.jacobianFcns.dpDDOTdL(Xhat,Uhat,hVars)];
                    dRLdX = [obj.omega.jacobianFcns.dRLdp(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dRLdpDOT(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dRLdRL(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dRLdL(Xhat,Uhat,hVars)];
                    dRLdU= [obj.omega.jacobianFcns.dRLdRLdot(Xhat,Uhat,hVars)*...
                        obj.cableConstraintMatrix,...
                        zeros(size(obj.omega.C,1),size(obj.rodConstraintMatrix,2))];
                    dLdX = [obj.omega.jacobianFcns.dLdp(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dLdpDOT(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dLdRL(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dLdL(Xhat,Uhat,hVars)];
                    dLdU= [zeros(size(obj.omega.R,1),size(obj.cableConstraintMatrix,2)),...
                        obj.omega.jacobianFcns.dLdLdot(Xhat,Uhat,hVars)*...
                        obj.rodConstraintMatrix];
                    
                    %TESTING JACOBIAN ACCURACY
                    [~,~,~,~,~,~,dpDDOTdX] = lsqnonlin(...
                        @(x0)pDDOTparser(obj,x0),...
                        [Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L],[],[],...
                        optimset('MaxIter',0,'Algorithm',...
                        'levenberg-marquardt',...
                        'Display','off'));
                    
                    %dynamic Jacobians
                    dfdx = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L);
                    dfdx(1:obj.nX_p,obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = ...
                        eye(obj.nX_pDOT);
                    dfdx(obj.nX_p+1:end,:) = [dpDDOTdX;dRLdX;dLdX];
                    dfdu = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L,...
                        size(obj.cableConstraintMatrix,2)+size(obj.rodConstraintMatrix,2));
                    dfdu(end-(obj.nX_RL+obj.nX_L)+1:end,:) = ...
                        [dRLdU;dLdU];
                    
                    % form large augmented matrices
                    obj.Amats(:,:,t) = [eye(obj.nX)+obj.dT*dfdx,...
                        zeros(obj.nX,1); zeros(1,obj.nX),1];
                    obj.Bmats(:,:,t) = [obj.dT*dfdu;zeros(1,size(dfdu,2))];
                    decay1 = obj.discount1;
                    decay2 = obj.discount2;
                    obj.RPrimemats(:,:,t) = decay2^(t-1)*obj.inputChangePenalty;
                    obj.QPrimemats(:,:,t) = [decay1^(t-1)*Q,zeros(obj.nX+1,obj.nU),...
                        decay1^(t-1)*Q*[obj.xTraj(:,t)-Xref;1];
                        zeros(obj.nU,obj.nX+1),decay2^(t-1)*obj.R,...
                        decay2^(t-1)*obj.R*obj.uGuess(:,t);...
                        (decay1^(t-1)*Q*[obj.xTraj(:,t)-Xref;1])',...
                        (decay2^(t-1)*obj.R*obj.uGuess(:,t))',...
                        [obj.xTraj(:,t)-Xref;1]'*decay1^(t-1)*Q*[obj.xTraj(:,t)-Xref;1]+...
                        obj.uGuess(:,t)'*decay2^(t-1)*obj.R*obj.uGuess(:,t)+...
                        obj.uDeltaGuess(:,t)'*decay2^(t-1)*obj.RPrimemats(:,:,t)*obj.uDeltaGuess(:,t)];
                    obj.NPrimemats(:,:,t) = [zeros(obj.nX+1+obj.nU,obj.nU);
                        obj.uDeltaGuess(:,t)'*obj.RPrimemats(:,:,t)];
                    
                    % make Q positive definite
                    obj.QPrimemats(:,:,t) = obj.QPrimemats(:,:,t)+(abs(min([real(eig(obj.QPrimemats(:,:,t)));0]))+1e-3)*eye(size(obj.QPrimemats(:,:,t)));
                end
                
                % terminal cost
                SF = 1;
                obj.QPrimemats(:,:,N+1) = SF*[decay1^(N)*Q,zeros(obj.nX+1,obj.nU),...
                        decay1^(N)*Q*[obj.xTraj(:,N+1)-Xref;1];
                        zeros(obj.nU,obj.nX+1),decay2^(N)*obj.R,...
                        decay2^(N)*obj.R*obj.uGuess(:,N+1);...
                        (decay1^(N)*Q*[obj.xTraj(:,N+1)-Xref;1])',...
                        (decay2^(N)*obj.R*obj.uGuess(:,N+1))',...
                        [obj.xTraj(:,N+1)-Xref;1]'*decay1^(N)*Q*[obj.xTraj(:,N+1)-Xref;1]+...
                        obj.uGuess(:,N+1)'*decay2^(N)*obj.R*obj.uGuess(:,N+1)];
                
                % make Q positive definite
                obj.QPrimemats(:,:,N+1) = obj.QPrimemats(:,:,N+1)+...
                    (abs(min([real(eig(obj.QPrimemats(:,:,t)));0]))+1e-3)*...
                    eye(size(obj.QPrimemats(:,:,N+1)));

                % calculate LQR feedback gains
                % backward pass
                obj.P{N+1} = obj.QPrimemats(:,:,N+1); %terminal cost
                for k = N:-1:1
                    %stage cost
                    Aprime = blkdiag([obj.Amats(:,:,k),obj.Bmats(:,:,k);
                        zeros(obj.nU,obj.nX+1),eye(obj.nU)],1);
                    Bprime = [obj.Bmats(:,:,k);eye(obj.nU);zeros(1,obj.nU)];
                    obj.P{k} = Aprime'*obj.P{k+1}*Aprime-...
                        (Aprime'*obj.P{k+1}*Bprime+obj.NPrimemats(:,:,k))/...
                        (obj.RPrimemats(:,:,k)+...
                        Bprime'*obj.P{k+1}*Bprime)*...
                        (Bprime'*obj.P{k+1}*Aprime+obj.NPrimemats(:,:,k)')+...
                        obj.QPrimemats(:,:,k);
                    %optimal LQR feedback gain
                    obj.Kgains(:,:,k) = (obj.RPrimemats(:,:,k)+...
                        Bprime'*obj.P{k+1}*Bprime)\...
                        (Bprime'*obj.P{k+1}*Aprime+...
                        obj.NPrimemats(:,:,k)');
                    
                    % constrain rod inputs
                    obj.Kgains(end-obj.nU_Ldot+1:end,:,k) = zeros(obj.nU_Ldot,obj.nX+obj.nU+2);
                end
                
                % calculate initial trajectory cost
                cost = calculateCurrentCost(obj,Xref,Q,R);
                
                uTempTop = obj.uGuess;
                uDeltaTempTop = obj.uDeltaGuess;
                xTemp = obj.xTraj;
                
                % BACKTRACKING LINESEARCH
                while(true) %do continuously
                    deltaX = zeros(size(obj.xTraj)); %X0 is fixed
                    deltaU = zeros(size(obj.uGuess));
                    for t = 1:N
                        % check input constraint violation
                        checkConstraints = 1;
                        
                        % initialize empty arrays of indices 
                        velExceeded=[];
                        j=[];
                        
                        % constant term update
                        obj.uDeltaGuess(:,t) = obj.uDeltaGuess(:,t)-Alpha*obj.Kgains(:,end,t);
                        
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
                        Khat(j,:) = zeros(numel(j),obj.nX+obj.nU+1);
                        Khat(end-obj.nU_Ldot+1:end,:) = zeros(obj.nU_Ldot,obj.nX+obj.nU+1);
                        
                        % feedback term update
                        obj.uDeltaGuess(:,t) = obj.uDeltaGuess(:,t)-Khat*[deltaX(:,t);0;deltaU(:,t)];
                        
                        if(checkConstraints)
                            % handle paired/similar cable constraints
                            velExceeded = abs(obj.cableConstraintMatrix*obj.uGuess(1:obj.nU_RLdot,t))>=obj.omega.cables.linear_velocity;
                            cables = find(velExceeded);
                            [~,j] = find(obj.cableConstraintMatrix(velExceeded,:));
                            for k = 1:numel(j)
                                obj.uGuess(j(k),t) = sign(obj.uGuess(j(k),t))*...
                                    min(abs(obj.uGuess(j(k),t)),...
                                    obj.omega.cables.linear_velocity(cables(k)));
                            end
                        end
                        
                        % update delta states
                        tempX = obj.xTraj(:,t+1);
                        obj.xTraj(:,t+1) = stepForward(obj,obj.uGuess(:,t),obj.xTraj(:,t));
                        deltaX(:,t+1) = obj.xTraj(:,t+1)-tempX;
                        tempU = obj.uGuess(:,t+1);
                        obj.uGuess(:,t+1) = obj.uGuess(:,t) +...
                            obj.uDeltaGuess(:,t); %update next uGuess at time t
                        deltaU(:,t+1) = obj.uGuess(:,t+1)-tempU;   
                    end          
                    
                    %check last input at t=N+1
                    if(checkConstraints)
                        % handle paired/similar cable constraints
                        velExceeded = abs(obj.cableConstraintMatrix*obj.uGuess(1:obj.nU_RLdot,N+1))>=obj.omega.cables.linear_velocity;
                        cables = find(velExceeded);
                        [~,j] = find(obj.cableConstraintMatrix(velExceeded,:));
                        for k = 1:numel(j)
                            obj.uGuess(j(k),N+1) = sign(obj.uGuess(j(k),N+1))*...
                                min(abs(obj.uGuess(j(k),N+1)),...
                                obj.omega.cables.linear_velocity(cables(k)));
                        end
                    end
                                            
                    %calculate new trajectory cost
                    currCost = calculateCurrentCost(obj,Xref,Q,R);
                    
                    % keep track of best cost
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
                        figure(3)
                        cla
                        plot(obj.uDeltaGuess')
                    end
                       
%                     udiff = norm(obj.uGuess-uTempTop)
%                     uDeltaDiff = norm(obj.uDeltaGuess-uDeltaTempTop)
                    if(currCost<cost)
                        %disp('Finished Line Search')
                        %                         Alpha
                        break
%                     if(checkArmijo(obj,deltaX,obj.uGuess-uTempTop,...
%                             obj.uDeltaGuess-uDeltaTempTop,Xref,Q,R,...
%                             cost,currCost,0.5))
%                         %disp('Finished Line Search')
%                         %                         Alpha
%                         break
%                     elseif(Alpha<1e-2)
%                         disp('Alpha small, moving on')
%                         costdiff = cost-currCost
%                         converged=1;
%                         break
                    else
                        %reset inputs, reduce line search stepsize
                        obj.uGuess = uTempTop;
                        obj.uDeltaGuess = uDeltaTempTop;
                        obj.xTraj = xTemp;
                        Alpha = Alpha/2;
                        if(Alpha<1e-2)
                            disp('Alpha small, moving on')
                            converged=1;
                            break
                        end
                    end
                end  
%                 if(max(max(abs(obj.uGuess(:,:)-uTempTop(:,:))))<1e-3)
%                     converged = 1;
%                     disp('First few inputs converged, moving on...')
%                 end   
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
            options = odeset('RelTol',1e-3,'AbsTol',1e-3,...
                'InitialStep',obj.initial_step_guess);
            XIN = [X.p;X.pDOT;X.RL;X.L];
            [t,XOUT] = ode23(@(t,Xstacked)...
                tensegrityODE(t, Xstacked ,Uinput, obj.omega.nominalFcn,...
                obj.omega.hFcns), [0 obj.dT], XIN, options);
            obj.initial_step_guess = t(end)-t(end-4);
            XOUT = XOUT(end,:)';
        end
        
        function cost = calculateCurrentCost(obj,Xref,Q,R)
            N = obj.horizon;
            %calculate trajectory cost
            z = [obj.xTraj(:,N+1)-Xref;1];
            cost = z'*obj.discount1^N*obj.Q*z +...
                obj.discount2^N*obj.uGuess(:,end)'*obj.R*obj.uGuess(:,end);
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1];
                cost = cost + z'*obj.discount1^(t-1)*Q*z +...
                    obj.uGuess(:,t)'*obj.discount2^(t-1)*R*obj.uGuess(:,t) +...
                    obj.uDeltaGuess(:,t)'*obj.discount2^(t-1)*obj.RPrimemats(:,:,t)*obj.uDeltaGuess(:,t);
            end
        end
        
        function satisfied = checkArmijo(obj,xDelta,uDelta,uChangeDelta,Xref,Q,R,prevCost,currCost,ratio)
            N = obj.horizon;
            %calculate trajectory cost
            J = prevCost;
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1]; %xTraj,uguess updated by this point
                check = ratio*(2*(z-[xDelta(:,t);0])'*obj.discount1^(t-1)*Q*[xDelta(:,t);0] +...
                    2*(obj.uGuess(:,t)-uDelta(:,t))'*obj.discount2^(t-1)*R*uDelta(:,t)+...
                    2*(obj.uDeltaGuess(:,t)-uChangeDelta(:,t))'*obj.discount2^(t-1)*obj.RPrimemats(:,:,t)*obj.uDeltaGuess(:,t));
                J = J + ratio*(2*(z-[xDelta(:,t);0])'*obj.discount1^(t-1)*Q*[xDelta(:,t);0] +...
                    2*(obj.uGuess(:,t)-uDelta(:,t))'*obj.discount2^(t-1)*R*uDelta(:,t)+...
                    2*(obj.uDeltaGuess(:,t)-uChangeDelta(:,t))'*obj.discount2^(t-1)*obj.RPrimemats(:,:,t)*obj.uDeltaGuess(:,t));
                if(check>0)
                    disp('error')
                end
            end
            
            satisfied = currCost<J;
        end
        
        function XOUT = discreteStepForward(obj,u_vec,x_vec,hVars)
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
            
            % helper variables
            for i = 1:length(obj.omega.hFcns.z)
                hVars.z{i} = obj.omega.hFcns.z{i}(X);
            end
            for i = 1:length(obj.omega.hFcns.v)
                hVars.v{i} = obj.omega.hFcns.v{i}(X);
            end
            hVars.RhatRhat = obj.omega.hFcns.RhatRhat;
            hVars.Chat = obj.omega.hFcns.Chat;
            hVars.J = obj.omega.hFcns.J(X,Uinput,hVars);
            
            XOUT = [X.p+obj.dT*X.pDOT;
                X.pDOT+obj.dT*obj.omega.nominalFcn.pDDOT(X.p,X.pDOT,X.RL,X.L);
                X.RL+obj.dT*obj.omega.nominalFcn.RLdot(X,Uinput,hVars);
                X.L+obj.dT*obj.omega.nominalFcn.Ldot(X,Uinput,hVars);
                ];
        end
        
        function pDDOT = pDDOTparser(obj,x0)
                        % parse x
                        idx=0;
                        Xhat.p = x0(idx+1:obj.nX_p);
                        idx = idx + obj.nX_p;
                        Xhat.pDOT = x0(idx+1:idx+obj.nX_pDOT);
                        idx = idx + obj.nX_pDOT;
                        Xhat.RL = x0(idx+1:idx+obj.nX_RL);
                        idx = idx + obj.nX_RL;
                        Xhat.L = x0(idx+1:end);
            
                        pDDOT = obj.omega.nominalFcn.pDDOT(Xhat.p,Xhat.pDOT,...
                            Xhat.RL,Xhat.L);
        end
    end
end

