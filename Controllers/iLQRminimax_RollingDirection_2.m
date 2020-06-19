classdef iLQRminimax_RollingDirection_2 < handle
    %QP_MPC_ROLLINGDIRECTION: This controller implements a minimax iterative
    %Dynamic Game, where optimal control inputs are obtained in the
    %presence of an adversarial disturbance
    %
    %   Inputs
    
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
        G
        P
        nX_p
        nX_pDOT
        nX_RL
        nX_L
        nX
        nU_RLdot
        nU_Ldot
        nU
        nW
        
        cableConstraintMatrix
        rodConstraintMatrix
        
        %iDG
        uGuess
        wGuess
        xTraj
        Amats
        Bmats
        Qmats
        Mmats
        Nmats
        Kgains
        wBounds
        
        %value function arrays for iDG
        V
        Vx
        Vxx     
        
        discount
    end
    
    methods
        function obj = iLQRminimax_RollingDirection_2(...
                X,omega,simTimeStep,horizon)
            %QP_MPC_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
            % controller timestep and horizon
            obj.dT = simTimeStep;
            obj.horizon = horizon;
            
            % initialize origin as default target of [0,0]
            obj.targetDestination = [0 0]';
            
            % actuator limits
            obj.omega = omega;
            obj.cableLinearVelocity = omega.cables.linear_velocity;
            obj.cableMinLength = omega.cables.minLength;
            obj.cableMaxLength = omega.cables.maxLength;
            
            % important dimensions
            obj.nX_p = numel(X.p);
            obj.nX_pDOT = numel(X.pDOT);
            obj.nX_RL = numel(X.RL);
            obj.nX_L = numel(X.L);
            obj.nX = obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L; %full state dimension
            obj.nU_RLdot = size(omega.cableConstraintMatrix,2); %constrained cable inputs
            obj.nU_Ldot = size(omega.rodConstraintMatrix,2); %constrained rod inputs
            obj.nU = obj.nU_RLdot + obj.nU_Ldot; %full input dimension
            obj.nW = obj.nX; %full disturbance dimension
            
            % input constraint mappings
            obj.cableConstraintMatrix = omega.cableConstraintMatrix;
            obj.rodConstraintMatrix = omega.rodConstraintMatrix;
            
            %state penalty matrix
            obj.Q = zeros(obj.nX+1);
            
            % cable deviation penalty matrix (the state x.RL should be the
            % deviation from XRef.RL or actual state itself)
            cableDeviationPenalty = 1e-1;
            obj.Q(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL,...
                obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
                cableDeviationPenalty*eye(obj.nX_RL); %1e0 works for 6-bar,
            
            % velocity input penalty
            cableVel_cost = 5e-4;
            obj.R = cableVel_cost*blkdiag(...
                omega.cableConstraintMatrix'*eye(obj.nX_RL)*omega.cableConstraintMatrix,...
                eye(obj.nX_L)); %~1e-2 works for dT=1e-3; 5e1 works for dT=5e-3
            rodVel_cost = 1e6;
            obj.R(obj.nU_RLdot+1:end,obj.nU_RLdot+1:end) =...
                rodVel_cost*eye(obj.nX_L); %penalize rod actuation heavily
            
            % Disturbance Penalty
            obj.G = 1e-2*eye(obj.nX);
                  
            % Array of Input Feedback Gains
            obj.Kgains = zeros(obj.nU+obj.nW,obj.nX+2,horizon+1);
                                    
            % Initialize iLQR input variables
            obj.uGuess = 1e-2*randn(size(omega.cableConstraintMatrix,2)+obj.nX_L,horizon);
%             obj.uGuess = 0*randn(size(omega.cableConstraintMatrix,2)+obj.nX_L,horizon);
            obj.uGuess(end-obj.nX_L+1:end,:)=zeros(obj.nX_L,size(obj.uGuess,2));
            
            % Initialize iLQR disturbance variables
            obj.wGuess = zeros(obj.nX,horizon);
            obj.wBounds = zeros(obj.nX,1);
            obj.wBounds(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
                1e-2*ones(obj.nX_RL,1); % RL disturbances
            obj.wBounds(1:obj.nX_p) = ...
                1e-2*ones(obj.nX_p,1); % position disturbances
            
            % Forward simulate to get initial state trajectory
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! might not need this
            obj.xTraj = zeros(obj.nX,horizon+1);
            for k = 1:horizon
                obj.xTraj(:,k+1) = stepForward(obj,obj.uGuess(:,k),obj.wGuess(:,k),obj.xTraj(:,k));
            end
            
            % Value function
            obj.V = zeros(1,horizon+1);
            obj.Vx = zeros(obj.nX+1,horizon+1);
            obj.Vxx = zeros(obj.nX+1,obj.nX+1,horizon+1);
            
            % Discount Factor for receding horizon control
            obj.discount = 1.0;
        end
        
        function setTargetDestination(obj,target)
            %SETTARGETDESTINATION This function changes the internal
            %property 'targetDestination' of the controller object
            
            obj.targetDestination = target;
        end
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                costOutput,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,nominalFnc,jacobianFcns,hFcns,costFcnHandle,...
                debugFcns,openLoopFlag)
                                    
            %calculate desired rolling direction according to target goal
            totalMass = sum(obj.omega.M);
            xCOM = sum(obj.omega.M.*Xhat.p(1:3:end))/totalMass;
            yCOM = sum(obj.omega.M.*Xhat.p(2:3:end))/totalMass;
            zCOM = sum(obj.omega.M.*Xhat.p(3:3:end))/totalMass;
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
                        
            %z-axis correction velocity
            zWeight = sum(obj.omega.M.*obj.omega.X.p0(3:3:end))/totalMass-zCOM;
            desiredDirection = [desiredDirection;0*zWeight]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changed to 0 here
            desiredDirection = desiredDirection/norm(desiredDirection);
            
            %linear penalty (mass-weighted velocity/momentum),
            velReward = 5e0;
            obj.Q(end,obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT) =...
                -velReward/2*desiredDirection(1)*obj.omega.M/totalMass;
            obj.Q(end,obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT) =...
                -velReward/2*desiredDirection(2)*obj.omega.M/totalMass;
            obj.Q(end,obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT) =...
                -velReward/2*desiredDirection(3)*obj.omega.M/totalMass;
            obj.Q(obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT,end) =...
                -velReward/2*desiredDirection(1)*obj.omega.M/totalMass;
            obj.Q(obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT,end) =...
                -velReward/2*desiredDirection(2)*obj.omega.M/totalMass;
            obj.Q(obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT,end) =...
                -velReward/2*desiredDirection(3)*obj.omega.M/totalMass;
                            
            %store rolling direction vector for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            %reference state for deviations
            Xref = [zeros(obj.nX_p,1);  % N/A, unweighted penalty
                zeros(obj.nX_pDOT,1);   % unweighted penalty
                obj.omega.X.RL0;        % N/A neutral cable lengths
                zeros(obj.nX_L,1)];     % N/A unweighted penalty
            
            %concatenated state
            Xcombined = [Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L];
            
            %call iterative Dynamic Game function
            [U_OL,X_OL,costOutput,hVars] = obj.iLQRminimax(...
                Xcombined,obj.Q,obj.R,obj.G,Xref);
            
            %remap U input vector (potentially smaller due to constraints)
            %to full cable/rod actuation
            U_desired.RLdot = obj.cableConstraintMatrix*...
                U_OL(1:size(obj.cableConstraintMatrix,2),1);
            U_desired.Ldot = obj.rodConstraintMatrix*...
                U_OL(size(obj.cableConstraintMatrix,2)+1:end,1);
            
            %cable velocity limits
            if(any(abs(U_desired.RLdot)>obj.cableLinearVelocity))
                disp(['WARNING: Desired Cable velocity ',...
                    'exceeds motor capability'])
                
                %                 disp('Clipping cable velocity input')
                %                 U_desired.RLdot = sign(U_desired.RLdot).*...
                %                     min(abs(U_desired.RLdot),...
                %                     obj.cableLinearVelocity);
                %
                %                 disp('Normalizing Velocity Input according to Motor Limits~~~~~~~~~~~~~~~~')
                %                 U_desired.RLdot = U_desired.RLdot/max(abs(U_desired.RLdot)).*obj.cableLinearVelocity*1.5;
            end
            
            %prepare open-loop outputs
            N = obj.horizon;
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
            else %send empty vectors, if not requested
                OL_states = [];
                OL_inputs = [];
            end
        end
        
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
        end
        
        
        function [U_optimal,OL_states,currCost,hVars] = iLQRminimax(...
                obj,X0,Q,R,G,Xref)
            % ILQRMINIMAX This function performs the bulk of the
            % calculations for the iterative dynamics games. Using an
            % iterative process, it calculates inputs/disturbances thath
            % have reached some sort of dynamic game equilibrium.
            
            %initial guess
            N = obj.horizon;
            %obj.uGuess = repmat(obj.uGuess(:,2),1,N+1);
            
            converged = 0; % initialize convergence flag to 'false'
            N = obj.horizon;
            obj.xTraj(:,1) = X0; % initial state
                        
            % forward simulate to get initial state trajectory
            for k = 1:N
                obj.xTraj(:,k+1) = stepForward(obj,obj.uGuess(:,k),...
                    obj.wGuess(:,k),obj.xTraj(:,k));
            end
            
            % initialize Q-values for linearized approximation
            Jx = zeros(obj.nX+1,N);
            Ju = zeros(obj.nU,N);
            Jw = zeros(obj.nW,N);
            Jxx = zeros(obj.nX+1,obj.nX+1,N);
            Jux = zeros(obj.nU,obj.nX+1,N);
            Jwx = zeros(obj.nW,obj.nX+1,N);
            Juu = zeros(obj.nU,obj.nU,N);
            Jww = zeros(obj.nW,obj.nW,N);
            Juw = zeros(obj.nU,obj.nW,N);
                                        
            bestCost = 1e8;
            Alpha = 1;
            
            % initialize value functions along current trajectory
            obj.V(end) = obj.discount^(N+1-1)*...
                [obj.xTraj(:,end)-Xref;1]'*Q*[obj.xTraj(:,end)-Xref;1];
            obj.Vx(:,end) = obj.discount^(N+1-1)*2*Q*[obj.xTraj(:,end)-Xref;1];
            obj.Vxx(:,:,end) = obj.discount^(N+1-1)*2*Q;
            
%             % populate value function initial guesses
%             for t = N:-1:1  
%                 obj.V(t) = obj.discount^(t-1)*(...
%                     [obj.xTraj(:,t)-Xref;1]'*Q*[obj.xTraj(:,t)-Xref;1]...
%                     + obj.uGuess(:,t)'*R*obj.uGuess(:,t)...
%                     - obj.wGuess(:,t)'*G*obj.wGuess(:,t)) + obj.V(t+1);
%             end
                
            %ITERATIVE LINEARIZATION OUTER-LOOP
            while(~converged) 
                Alpha = 1;%min(1,Alpha*2); %increase linesearch stepsize
                       
                %record 
                VTemp = obj.V;
                VxTemp = obj.Vx;
                VxxTemp = obj.Vxx;
                
                %                 Alpha = 1;
                for t = N:-1:1 % BACKWARDS PASS
                    % NOTE: t=1 is last for backwards path
                    % thus hVars at t=t0 can be passed as output to dynamics
                    
                    % parse x at timestep t
                    idx=0;
                    Xhat.p = obj.xTraj(idx+1:obj.nX_p,t);
                    idx = idx + obj.nX_p;
                    Xhat.pDOT = obj.xTraj(idx+1:idx+obj.nX_pDOT,t);
                    idx = idx + obj.nX_pDOT;
                    Xhat.RL = obj.xTraj(idx+1:idx+obj.nX_RL,t);
                    idx = idx + obj.nX_RL;
                    Xhat.L = obj.xTraj(idx+1:end,t);
                    % parse u at timestep t
                    Uhat.RLdot = obj.uGuess(1:obj.nU_RLdot,t);
                    Uhat.Ldot = obj.uGuess(obj.nU_RLdot+1:end,t);
                    
                    % helper variables for dynamics calculations
                    for i = 1:length(obj.omega.hFcns.z)
                        hVars.z{i} = obj.omega.hFcns.z{i}(Xhat);
                    end
                    for i = 1:length(obj.omega.hFcns.v)
                        hVars.v{i} = obj.omega.hFcns.v{i}(Xhat);
                    end
                    hVars.RhatRhat = obj.omega.hFcns.RhatRhat;
                    hVars.Chat = obj.omega.hFcns.Chat;
                    hVars.J = obj.omega.hFcns.J(Xhat,Uhat,hVars);
                    
                    % linearized dynamics recalculated at each timestep
                    dpDDOTdX = [obj.omega.jacobianFcns.dpDDOTdp(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dpDDOTdpDOT(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dpDDOTdRL(Xhat,Uhat,hVars),...
                        obj.omega.jacobianFcns.dpDDOTdL(Xhat,Uhat,hVars)];
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
                        obj.omega.jacobianFcns.dLdLdot(Xhat,Uhat,hVars)*obj.rodConstraintMatrix];
                                        
                    %dynamics Jacobian w.r.t. states (Fx)
                    Abar = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L);
                    Abar(1:obj.nX_p,obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = ...
                        eye(obj.nX_pDOT);
                    Abar(obj.nX_p+1:end,:) = [dpDDOTdX;dRLdX;dLdX];
                    Fx = [eye(obj.nX) + obj.dT*Abar,zeros(obj.nX,1);
                        zeros(1,obj.nX),1];
                    
                    %dynamics Jacobian w.r.t. input (Fu)
                    Bbar = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L,...
                        size(obj.cableConstraintMatrix,2)+size(obj.rodConstraintMatrix,2));
                    Bbar(end-(obj.nX_RL+obj.nX_L)+1:end,:) = ...
                        [dRLdU;dLdU];
                    Fu = [obj.dT*Bbar;zeros(1,obj.nU)];
                    
                    %dynamics Jacobian w.r.t. disturbance (Fw)
                    Fw = [obj.dT*eye(obj.nX);zeros(1,obj.nW)];
                    
                    %Calculate gradients/Jacobians of Q-function (J here)
                    Jx(:,t) = 2*obj.discount^(t-1)*Q*[obj.xTraj(:,t)-Xref;1] +...
                        Fx'*obj.Vx(:,t+1);
                    Ju(:,t) = 2*obj.discount^(t-1)*R*obj.uGuess(:,t) +...
                        Fu'*obj.Vx(:,t+1);
                    Jw(:,t) = -2*obj.discount^(t-1)*G*obj.wGuess(:,t) +...
                        Fw'*obj.Vx(:,t+1);
                    Jxx(:,:,t) = 2*obj.discount^(t-1)*Q + Fx'*obj.Vxx(:,:,t+1)*Fx;
                    Jux(:,:,t) = Fu'*obj.Vxx(:,:,t+1)*Fx;
                    Jwx(:,:,t) = Fw'*obj.Vxx(:,:,t+1)*Fx;
                    Juu(:,:,t) = 2*obj.discount^(t-1)*R + Fu'*obj.Vxx(:,:,t+1)*Fu;
                    Jww(:,:,t) = -2*obj.discount^(t-1)*G + Fw'*obj.Vxx(:,:,t+1)*Fw;
                    Juw(:,:,t) = Fu'*obj.Vxx(:,:,t+1)*Fw;
                    
                    %regularize Hessians
                    rho = 1e-3;
%                     Juu(:,:,t) = Juu(:,:,t) + (max(0,-min(real(eig(Juu(:,:,t)))))+rho)*eye(obj.nU);
%                     Jww(:,:,t) = Jww(:,:,t) + (min(0,-max(real(eig(Jww(:,:,t)))))-rho)*eye(obj.nW);
%                     Jxx(:,:,t) = Jxx(:,:,t) + (max(0,-min(real(eig(Jxx(:,:,t)))))+rho)*eye(obj.nX+1);
                    
                    %calculate optimal inputs and feedback gains
                    K_ut = inv(eye(obj.nU)-Juu(:,:,t)\Juw(:,:,t)/...
                        Jww(:,:,t)*Juw(:,:,t)')/Juu(:,:,t);
                    K_wt = inv(eye(obj.nW)-Jww(:,:,t)\Juw(:,:,t)'/...
                        Juu(:,:,t)*Juw(:,:,t))/Jww(:,:,t);
%                     K_ut = inv((eye(obj.nU)-Juu(:,:,t)\Juw(:,:,t)/...
%                         Jww(:,:,t)*Juw(:,:,t)')/Juu(:,:,t));
%                     K_wt = inv((eye(obj.nW)-Jww(:,:,t)\Juw(:,:,t)'/...
%                         Juu(:,:,t)*Juw(:,:,t))/Jww(:,:,t));
                    g_ut = K_ut*(Juw(:,:,t)/Jww(:,:,t)*Jw(:,t)-Ju(:,t));
                    g_wt = K_wt*(Juw(:,:,t)'/Juu(:,:,t)*Ju(:,t)-Jw(:,t));
                    G_ut = K_ut*(Juw(:,:,t)/Jww(:,:,t)*Jwx(:,:,t)-Jux(:,:,t));
                    G_wt = K_wt*(Juw(:,:,t)'/Juu(:,:,t)*Jux(:,:,t)-Jwx(:,:,t));

                    %optimal LQR feedback gain
                    obj.Kgains(1:obj.nU,:,t) = [G_ut,g_ut];
                    obj.Kgains(obj.nU+1:end,:,t) = [G_wt,g_wt];
                    
%                     %debugging for now
%                     obj.Kgains(obj.nU+1:end,:,t) = 0;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                    %zero out gain on rods (fix this later, should look at
                    %rod velocity constraints instead of hardcoding)
                    obj.Kgains(obj.nU_RLdot+1:obj.nU_RLdot+obj.nU_Ldot,:,t) = zeros(obj.nU_Ldot,obj.nX+2);
                end
 
                %calculate initial trajectory cost
                cost = calculateCost(obj,Xref,Q,R,G);
                                
                %save initial input/disturbance guesses before backtracking
                %linesearch
                uTemp = obj.uGuess;
                wTemp = obj.wGuess;
                xTemp = obj.xTraj;
                
                %BACKTRACKING LINESEARCH INNER-LOOP
                velExceeded=[];
                distExceeded=[];
                j=[];
                while(true) 
                    deltaX = zeros(size(obj.xTraj)); 
                    
                    % update value function quadratic approximation
                    obj.V(end) = obj.discount^(N+1-1)*...
                        [obj.xTraj(:,end)-Xref;1]'*Q*[obj.xTraj(:,end)-Xref;1];
                    obj.Vx(:,end) = obj.discount^(N+1-1)*2*Q*[obj.xTraj(:,end)-Xref;1];
                    obj.Vxx(:,:,end) = obj.discount^(N+1-1)*2*Q;
                    for t = N:-1:1
                        g_ut = Alpha*obj.Kgains(1:obj.nU,end,t);
                        G_ut = obj.Kgains(1:obj.nU,1:end-1,t);
                        %augment due to constraints
                        g_ut(j,:) = zeros(numel(j),1);
                        G_ut(j,:) = zeros(numel(j),size(G_ut,2));
                        g_wt = Alpha*obj.Kgains(obj.nU+1:end,end,t);
                        G_wt = obj.Kgains(obj.nU+1:end,1:end-1,t);
                        %augment due to constraints
                        g_wt(distExceeded,:) = zeros(sum(distExceeded),1);
                        G_wt(distExceeded,:) = zeros(sum(distExceeded),size(G_wt,2));
                      
                        obj.V(t) = obj.V(t) + ...
                            (g_ut'*Ju(:,t) + g_wt'*Jw(:,t)) + ...
                            g_ut'*Juw(:,:,t)*g_wt + 1/2*...
                            (g_ut'*Juu(:,:,t)*g_ut + g_wt'*Jww(:,:,t)*g_wt);
                        obj.Vx(:,t) = Jx(:,t) + G_ut'*Ju(:,t) + G_wt'*Jw(:,t)...
                            + G_ut'*Juu(:,:,t)*g_ut + Jux(:,:,t)'*g_ut +...
                            Jwx(:,:,t)'*g_wt + G_wt'*Jww(:,:,t)*g_wt +...
                            G_wt'*Juw(:,:,t)'*g_ut + G_ut'*Juw(:,:,t)*g_wt;
                        %                     obj.Vxx(:,:,t) = 1/2*(Jxx(:,:,t) +...
                        %                         G_ut'*Juu(:,:,t)*G_ut + G_wt'*Jww(:,:,t)*G_wt) +...
                        %                         G_ut'*Jux(:,:,t) + G_wt'*Jwx(:,:,t) + G_ut'*Juw(:,:,t)*G_wt;
                        obj.Vxx(:,:,t) = 1/2*(Jxx(:,:,t) +...
                            G_ut'*Juu(:,:,t)*G_ut + G_wt'*Jww(:,:,t)*G_wt +...
                            G_ut'*Jux(:,:,t) + G_wt'*Jwx(:,:,t) + G_ut'*Juw(:,:,t)*G_wt +...
                            (G_ut'*Jux(:,:,t) + G_wt'*Jwx(:,:,t) + G_ut'*Juw(:,:,t)*G_wt)');
%                         
%                         obj.V(t) = obj.V(t+1) - Ju(:,t)'/Juu(:,:,t)*Ju(:,t) -...
%                             Jw(:,t)'/Jww(:,:,t)*Jw(:,t);
%                         obj.Vx(:,t) = Jx(:,t) - Ju(:,t)'*Juu(:,:,t)*Jux(:,t) -...
%                             Jw(:,t)'*Jww(:,:,t)*Jwx(:,t);
%                         obj.Vxx(:,:,t) = Jxx(:,:,t) - ...
%                             Jwx(:,:,t)'/Jww(:,:,t)*Jwx(:,:,t) -...
%                             Jux(:,:,t)'/Juu(:,:,t)*Jux(:,:,t);
                    end
                    
                    for t = 1:N                                                
                        %add affine terms to input/disturbance
                        obj.uGuess(:,t) = obj.uGuess(:,t)+Alpha*obj.Kgains(1:size(obj.uGuess,1),end,t);
                        obj.wGuess(:,t) = obj.wGuess(:,t)+Alpha*obj.Kgains(size(obj.uGuess,1)+1:end,end,t);
                        
                        %check input constraint violation
                        checkConstraints = 1;
                        
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
                            %handle disturbance bounds
                            distExceeded = abs(obj.wGuess(:,t))>obj.wBounds;
                            obj.wGuess(distExceeded,t) = sign(...
                                obj.wGuess(distExceeded,t)).*...
                                obj.wBounds(distExceeded);
                        end
                        Khat = obj.Kgains(:,1:end-1,t);
                        Khat(j,:) = zeros(numel(j),obj.nX+1); 
                        obj.Kgains(j,:,t) = zeros(numel(j),obj.nX+2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Khat(end-obj.nX_L+1:end,:) = zeros(obj.nX_L,obj.nX+1); %rods can't move
                        Khat(end-obj.nX+find(distExceeded),:) = zeros(sum(distExceeded),obj.nX+1); %zero out disturbances at active constraints
                        
                        %feedback update for input/disturbance
                        obj.uGuess(:,t) = obj.uGuess(:,t)+Khat(1:size(obj.uGuess,1),:)*[deltaX(:,t);0];
                        obj.wGuess(:,t) = obj.wGuess(:,t)+Khat(size(obj.uGuess,1)+1:end,:)*[deltaX(:,t);0];
                        
                        %check constraints after feedback update
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
                            %handle disturbance bounds
                            distExceeded = abs(obj.wGuess(:,t))>obj.wBounds;
                            obj.wGuess(distExceeded,t) = sign(...
                                obj.wGuess(distExceeded,t)).*...
                                obj.wBounds(distExceeded);
                        end
                        
                        %calculate delta for next time step
                        temp = obj.xTraj(:,t+1);
                        obj.xTraj(:,t+1) = stepForward(...
                            obj,obj.uGuess(:,t),obj.wGuess(:,t),...
                            obj.xTraj(:,t));
                        deltaX(:,t+1) = obj.xTraj(:,t+1)-temp; 
                        
                    end                    
                    
                    % calculate new trajectory cost
                    currCost = calculateCost(obj,Xref,Q,R,G);
                    
                    % update best cost
                    bestCost = min(bestCost,currCost);
                    
                    if(1) %toggle on/off for debugging
                        costDiff = cost-currCost
                        figure(2)
%                         plot(obj.uGuess')
%                         hold on
                        subplot(1,2,1)
                        cla
                        plot(obj.uGuess','linewidth',2)
                        
                        subplot(1,2,2)
                        cla
                        plot(obj.wGuess')
%                         hold on
%                         plot(wTemp')
                        
                        drawnow()

                        title(['Best Cost = ',num2str(bestCost,'%10.5e'),...
                            '; Current Cost = ',num2str(currCost,'%10.5e'),...
                            '; Alpha = ',num2str(Alpha)])    
                    end
                      
                    Alpha
                    if(currCost < cost)
                        %disp('Finished Line Search')
                        %                         Alpha
%                         if(cost-currCost < 1e-3 )
%                             converged = true;
%                         end
                        break
%                     if(checkArmijo(obj,deltaX,obj.uGuess-uTemp,obj.wGuess-wTemp,Xref,Q,R,G,cost,currCost,0.5))
%                         %disp('Finished Line Search')
%                         %                         Alpha
%                         break
                    else
                        % reset inputs, reduce line search stepsize
                        obj.uGuess = uTemp;
                        obj.wGuess = wTemp;
                        obj.xTraj = xTemp;
%                         sum(abs(obj.V-VTemp))
%                         if(sum(abs(obj.V-VTemp))<1e-3)
%                             converged=1;
%                         end
                        Alpha = Alpha/2;
                        if(Alpha<1e-3)
                            disp('Alpha small, moving on')
                            converged = 1;
                            break
                        end
%                         obj.V = VTemp;
%                         obj.Vx = VxTemp;
%                         obj.Vxx = VxxTemp;
                        
%                         VTemp = obj.V;
%                         VxTemp = obj.Vx;
%                         VxxTemp = obj.Vxx;
                    end
                    
                end
                
%                 sum(abs(obj.V-VTemp))
%                 if(sum(abs(obj.V-VTemp))<1e-3)
%                     converged=1;
%                 end
                
                
                    
%                 if(max(max(abs(obj.uGuess(:,:)-uTemp(:,:))))<1e-4 &&...
%                         max(max(abs(obj.wGuess(:,:)-wTemp(:,:))))<1e-4)
%                     converged = 1;
%                     disp('First few inputs converged, moving on...')
%                 end                               
            end
            
            U_optimal = obj.uGuess;
            OL_states = obj.xTraj;
        end
        
        
        function XOUT = stepForward(obj,u_vec,w_vec,x_vec)
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

            %forward simulate dynamics with ODE solver
            XIN = [X.p;X.pDOT;X.RL;X.L];
            [~,XOUT] = ode45(@(t,Xstacked)...
                tensegrityODE(t,Xstacked,Uinput,obj.omega.nominalFcn,obj.omega.hFcns),[0 obj.dT],XIN);
            
            XOUT = XOUT(end,:)' + obj.dT*w_vec; %really simple disturbance model
        end
        
        function cost = calculateCost(obj,Xref,Q,R,G)
            N = obj.horizon;
            obj.discount;
            %calculate trajectory cost
            z = [obj.xTraj(:,N+1)-Xref;1];
            cost = (obj.discount^(N+1-1))*z'*obj.Q*z;
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1];
                cost = cost + (obj.discount^(t-1))*(z'*Q*z +...
                    obj.uGuess(:,t)'*R*obj.uGuess(:,t) -...
                    obj.wGuess(:,t)'*G*obj.wGuess(:,t));
            end
        end
        
        function cost = calculateCostNoDist(obj,Xref,Q,R)
            N = obj.horizon;
            %calculate trajectory cost
            z = [obj.xTraj(:,N+1)-Xref;1];
            cost = (obj.discount^(N+1-1))*z'*obj.Q*z;
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1];
                cost = cost + (obj.discount^(t-1))*(z'*Q*z +...
                    obj.uGuess(:,t)'*R*obj.uGuess(:,t));
            end
        end
        
        function satisfied = checkArmijo(obj,xDelta,uDelta,wDelta,Xref,...
                Q,R,G,prevCost,currCost,ratio)
            N = obj.horizon;
            %calculate trajectory cost
            J = prevCost; 
            for t = 1:N
                z = [obj.xTraj(:,t)-Xref;1]; %xTraj,uguess updated by this point
                J = J + ratio*(2*(z-[xDelta(:,t);0])'*Q*[xDelta(:,t);0] +...
                    2*(obj.uGuess(:,t)-uDelta(:,t))'*R*uDelta(:,t) -...
                    2*(obj.wGuess(:,t)-wDelta(:,t))'*G*wDelta(:,t));
            end  
            satisfied = currCost<J;
        end
        
    end
end

