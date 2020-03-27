classdef LQR_RollingDirection_2 < handle
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
    end
    
    methods
        function obj = LQR_RollingDirection_2(...
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
            % state and input dimensions
            obj.nX_p = numel(X.p);
            obj.nX_pDOT = numel(X.pDOT);
            obj.nX_RL = numel(X.RL);
            obj.nX_L = numel(X.L);
            obj.nX = obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L;
            obj.nU_RLdot = size(omega.cableConstraintMatrix,2); %constrained cable inputs
            obj.nU_RLdot = size(omega.cableConstraintMatrix,2); %constrained cable inputs
            % cable/rod input remapping
            obj.cableConstraintMatrix = omega.cableConstraintMatrix;
            obj.rodConstraintMatrix = omega.rodConstraintMatrix;       
            %== LQR cost function =========================================
            %state penalty
            obj.Q = zeros(obj.nX+1);
            %cable deviation penalty
            cableDeviationPenalty = 5e-1;
            obj.Q(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL,...
                obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
                cableDeviationPenalty*eye(obj.nX_RL); 
            %input penalty
            obj.R = 1e1*obj.dT*blkdiag(...
                omega.cableConstraintMatrix'*...
                eye(obj.nX_RL)*omega.cableConstraintMatrix,...
                omega.rodConstraintMatrix'*...
                eye(obj.nX_L)*omega.rodConstraintMatrix');
            obj.R(obj.nU_RLdot+1:end,obj.nU_RLdot+1:end) =...
                1e6*eye(obj.nX_L); %penalize rod actuation heavily
            %==============================================================
            %cell array of cost-to-go and feedback matrices
            obj.P = cell(horizon,1);
        end

        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                costOutput,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,nominalFnc,jacobianFcns,hFcns,costFcnHandle,...
                debugFcns,openLoopFlag)
                        
            N = obj.horizon;          
            Uhat.RLdot = zeros(size(Uhat.RLdot));
            Uhat.Ldot = zeros(size(Uhat.Ldot));
            
            %First calculate pDDOT and linearization about Xhat,Uhat
            %helper variables
            for i = 1:length(hFcns.z)
                hVars.z{i} = hFcns.z{i}(Xhat);
            end
            for i = 1:length(hFcns.v)
                hVars.v{i} = hFcns.v{i}(Xhat);
            end
            hVars.RhatRhat = hFcns.RhatRhat;
            hVars.Chat = hFcns.Chat;
            hVars.J = hFcns.J(Xhat,Uhat,hVars);
            %Nominal values for linearization
            pDDOT_bar = nominalFnc.pDDOT(Xhat.p,Xhat.pDOT,Xhat.RL,Xhat.L);
            RLdot_bar =  nominalFnc.RLdot(Xhat,Uhat,hVars);
            Ldot_bar =  nominalFnc.Ldot(Xhat,Uhat,hVars);
            %dynamic linearization recalculated each controller iteration
            %             dpDDOTdX = [jacobianFcns.dpDDOTdp(Xhat,Uhat,hVars),...
            %                 jacobianFcns.dpDDOTdpDOT(Xhat,Uhat,hVars),...
            %                 jacobianFcns.dpDDOTdRL(Xhat,Uhat,hVars),...
            %                 jacobianFcns.dpDDOTdL(Xhat,Uhat,hVars)];
            dpDDOTdU = [jacobianFcns.dpDDOTdRLdot(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdLdot(Xhat,Uhat,hVars)];
            dRLdX = [jacobianFcns.dRLdp(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdpDOT(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdL(Xhat,Uhat,hVars)];
            dRLdU= [jacobianFcns.dRLdRLdot(Xhat,Uhat,hVars)*...
                obj.cableConstraintMatrix,...
                zeros(size(obj.omega.C,1),...
                size(obj.rodConstraintMatrix,2))];
            dLdX = [jacobianFcns.dLdp(Xhat,Uhat,hVars),...
                jacobianFcns.dLdpDOT(Xhat,Uhat,hVars),...
                jacobianFcns.dLdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dLdL(Xhat,Uhat,hVars)];
            dLdU= [zeros(size(obj.omega.R,1),size(obj.cableConstraintMatrix,2)),...
                jacobianFcns.dLdLdot(Xhat,Uhat,hVars)*obj.rodConstraintMatrix];
            
            [~,~,~,~,~,~,dpDDOTdX] = lsqnonlin(...
                        @(x0)pDDOTparser(obj,x0),...
                        [Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L],[],[],...
                        optimset('MaxIter',0,'MaxFunEvals',0,'Algorithm',...
                        'levenberg-marquardt',...
                        'Display','off'));
            
            %calculate desired rolling direction according to target goal
            xCOM = mean(Xhat.p(1:3:end));
            yCOM = mean(Xhat.p(2:3:end));
            zCOM = mean(Xhat.p(3:3:end));
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            
            % vertical direction
            totalMass = sum(obj.omega.M);
            zWeight = 1.1*sum(obj.omega.M.*obj.omega.X.p0(3:3:end))/...
                totalMass-zCOM; % 10% higher than neutral position
            zSF = 10; % how much to weigh z-deviation
            desiredDirection = [desiredDirection;zSF*zWeight];
            desiredDirection = desiredDirection/norm(desiredDirection);
            
%             % maintain momentum
%             COM_vel = [sum(obj.omega.M.*Xhat.pDOT(1:3:end))/totalMass;
%                 sum(obj.omega.M.*Xhat.pDOT(2:3:end))/totalMass;
%                 0
%                 ];
%             COM_vel = COM_vel/norm(COM_vel);
%             % Combine desired direction and current robot momentum
%             Beta = 0.25;  % weight on current momentum
%             desiredDirection = (1-Beta)*desiredDirection + Beta*COM_vel;  
%             desiredDirection = desiredDirection/norm(desiredDirection);
            
            %linear penalty (velocity)
            velReward = 5e2;
%             obj.Q(end,obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT) =...
%                 -velReward/2*desiredDirection(1)*...
%                 obj.omega.M/sum(obj.omega.M);
%             obj.Q(end,obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT) =...
%                 -velReward/2*desiredDirection(2)*...
%                 obj.omega.M/sum(obj.omega.M);
%             obj.Q(end,obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT) =...
%                 -velReward/2*desiredDirection(3)*...
%                 obj.omega.M/sum(obj.omega.M);
%             obj.Q(obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT,end) =...
%                 -velReward/2*desiredDirection(1)*...
%                 obj.omega.M/sum(obj.omega.M);
%             obj.Q(obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT,end) =...
%                 -velReward/2*desiredDirection(2)*...
%                 obj.omega.M/sum(obj.omega.M);
%             obj.Q(obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT,end) =...
%                 -velReward/2*desiredDirection(3)*...
%                 obj.omega.M/sum(obj.omega.M);
            dir_diag = zeros(1,obj.nX_pDOT);
            dir_diag(1:3:end) = velReward*desiredDirection(1);
            dir_diag(2:3:end) = velReward*desiredDirection(2);
            dir_diag(3:3:end) = velReward*desiredDirection(3);
            obj.Q(obj.nX_p+1:obj.nX_p+obj.nX_pDOT,obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = ...
                diag(dir_diag);
            %store important values for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            %reference state
            Xref = [zeros(obj.nX_p,1);
                zeros(obj.nX_pDOT,1);
                obj.omega.X.RL0;        %neutral cable lengths
                obj.omega.X.L0          %neutral cable lengths      
                ];
                                        
            %setup LQR variables
            Abar = zeros(obj.nX);
            Abar(1:obj.nX_p,obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = ...
                eye(obj.nX_pDOT);  % pDOT
            Abar(obj.nX_p+1:end,:) = [dpDDOTdX;dRLdX;dLdX];
            Bbar = zeros(obj.nX, size(obj.cableConstraintMatrix,2) + ...
                size(obj.rodConstraintMatrix,2));
            Bbar(end-(obj.nX_RL+obj.nX_L)+1:end,:) = ...
                [dRLdU;dLdU];
            xDOTBar = [Xhat.pDOT;
                        pDDOT_bar;
                        RLdot_bar;
                        Ldot_bar]; % nominal state velocities
                    
            Xcombined = [Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L];
            z = obj.dT*(xDOTBar-Abar*(Xcombined - Xref)...
                -Bbar*[Uhat.RLdot;Uhat.Ldot]);
            Aaug = [eye(obj.nX)+obj.dT*Abar,z;
                        zeros(1,obj.nX),1];
            Baug = [obj.dT*Bbar;zeros(1,size(Bbar,2))];
                
            %backward pass
%             obj.Q = obj.Q + eye(size(obj.Q))*(min(eig(obj.Q))+1e-6);
            obj.P{obj.horizon} = obj.Q;
            for k = obj.horizon-1:-1:1
                obj.P{k} = Aaug'*obj.P{k+1}*Aaug-(Aaug'*obj.P{k+1}*Baug)/...
                    (obj.R+Baug'*obj.P{k+1}*Baug)*...
                    (Baug'*obj.P{k+1}*Aaug)+obj.Q;
            end

            % optimal feedback control
            K = (obj.R+Baug'*obj.P{2}*Baug)\...
                (Baug'*obj.P{2}*Aaug);
            U = -K*([Xcombined - Xref; 1]);
            
            % remap U input vector to full cable/rod actuation
            % potentially smaller dimensiondue to constraints)
            U_desired.RLdot = obj.cableConstraintMatrix*...
                U(1:size(obj.cableConstraintMatrix,2));
            U_desired.Ldot = obj.rodConstraintMatrix*...
                U(size(obj.cableConstraintMatrix,2)+1:end);
            
            if(any(abs(U_desired.RLdot)>obj.cableLinearVelocity))
                disp('WARNING: Desired Cable velocity exceeds motor capability')
                
%                 disp('Clipping cable velocity input')
%                 U_desired.RLdot = sign(U_desired.RLdot).*...
%                 min(abs(U_desired.RLdot),...
%                 obj.cableLinearVelocity);
            
%                 disp('Normalizing Velocity Input according to Motor Limits~~~~~~~~~~~~~~~~')
%                 U_desired.RLdot = U_desired.RLdot/max(abs(U_desired.RLdot)).*obj.cableLinearVelocity*1.5;
            end

                            
            OL_states = [];
            OL_inputs = [];
            costOutput =  [Xcombined-Xref;1]'*obj.P{1}*[Xcombined-Xref;1];

        end
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
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

