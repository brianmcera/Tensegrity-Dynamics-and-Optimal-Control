classdef LQR_RollingDirection_centralPayload < handle
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
        function obj = LQR_RollingDirection_centralPayload(...
                X,omega,simTimeStep,horizon)
            %QP_MPC_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
            %controller timestep and horizon
            obj.dT = simTimeStep;
            obj.horizon = horizon;
            
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
            %             %weigh desired velocities proportional to node mass
            %             obj.Q(obj.nX_p+1:obj.nX_p+obj.nX_pDOT,...
            %                 obj.nX_p+1:obj.nX_p+obj.nX_pDOT) =...
            %                 10*diag(kron(omega.M/sum(omega.M),[1 0 0]'));
            %             obj.Q(1:obj.nX_p,...
            %                 1:obj.nX_p) =...
            %                 5*diag(kron(omega.M/sum(omega.M),[1 0 0]'));
            %             obj.Q(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL,...
            %                 obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
            %                 0.1*eye(obj.nX_RL);
            %
            
            %this works for 6bar models without payload
            %             obj.Q(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL,...
            %                 obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
            %                 5e3*eye(obj.nX_RL); %cable deviation penalty
            
            %cable deviation penalty (careful, on if the state is the
            %deviation from RL0 or actual state itself)
            obj.Q(obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL,...
                obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
                1e0*eye(obj.nX_RL); %1e0 works for 6-bar,


            
            %velocity input penalty
            %             obj.R = 1e2*obj.dT*blkdiag(...
            %                 omega.cableConstraintMatrix'*eye(obj.nX_RL)*omega.cableConstraintMatrix,...
            %                 eye(obj.nX_L));
            obj.R = 5e-1*blkdiag(...
                omega.cableConstraintMatrix'*eye(obj.nX_RL)*omega.cableConstraintMatrix,...
                eye(obj.nX_L)); %~1e-2 works for dT=1e-3; 5e1 works for dT=5e-3
            obj.R(obj.nU_RLdot+1:end,obj.nU_RLdot+1:end) =...
                1e6*eye(obj.nX_L); %penalize rod actuation heavily
            
            %cell array of cost-to-go and feedback matrices
            obj.P = cell(horizon+1,1);
        end

        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                costOutput,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,nominalFnc,jacobianFcns,hFcns,costFcnHandle,...
                debugFcns,openLoopFlag)
                        
            N = obj.horizon;
                          
%             Uhat.RLdot = zeros(size(Uhat.RLdot));
%             Uhat.Ldot = zeros(size(Uhat.Ldot));
            
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
            dpDDOTdX = [jacobianFcns.dpDDOTdp(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdpDOT(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdL(Xhat,Uhat,hVars)];
            dpDDOTdU = [jacobianFcns.dpDDOTdRLdot(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdLdot(Xhat,Uhat,hVars)];
            dRLdX = [jacobianFcns.dRLdp(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdpDOT(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdL(Xhat,Uhat,hVars)];
            dRLdU= [jacobianFcns.dRLdRLdot(Xhat,Uhat,hVars)*obj.cableConstraintMatrix,...
                zeros(size(obj.omega.C,1),size(obj.rodConstraintMatrix,2))];
            dLdX = [jacobianFcns.dLdp(Xhat,Uhat,hVars),...
                jacobianFcns.dLdpDOT(Xhat,Uhat,hVars),...
                jacobianFcns.dLdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dLdL(Xhat,Uhat,hVars)];
            dLdU= [zeros(size(obj.omega.R,1),size(obj.cableConstraintMatrix,2)),...
                jacobianFcns.dLdLdot(Xhat,Uhat,hVars)*obj.rodConstraintMatrix];
            
            %calculate desired rolling direction according to target goal
            xCOM = mean(Xhat.p(1:3:end));
            yCOM = mean(Xhat.p(2:3:end));
            zCOM = mean(Xhat.p(3:3:end));
            centroid = [xCOM;yCOM;zCOM];
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            zWeight = mean(obj.omega.X.p0(3:3:end))-zCOM;
            desiredDirection = [desiredDirection;10*zWeight];
            desiredDirection = desiredDirection/norm(desiredDirection);
                                    
            %linear penalty (velocity)
            obj.Q(end,obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT) =...
                -5/2*desiredDirection(1)*obj.omega.M/sum(obj.omega.M);
            obj.Q(end,obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT) =...
                -5/2*desiredDirection(2)*obj.omega.M/sum(obj.omega.M);
            obj.Q(end,obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT) =...
                -5/2*desiredDirection(3)*obj.omega.M/sum(obj.omega.M);
            obj.Q(obj.nX_p+1:3:obj.nX_p+obj.nX_pDOT,end) =...
                -5/2*desiredDirection(1)*obj.omega.M/sum(obj.omega.M);
            obj.Q(obj.nX_p+2:3:obj.nX_p+obj.nX_pDOT,end) =...
                -5/2*desiredDirection(2)*obj.omega.M/sum(obj.omega.M);
            obj.Q(obj.nX_p+3:3:obj.nX_p+obj.nX_pDOT,end) =...
                -5/2*desiredDirection(3)*obj.omega.M/sum(obj.omega.M);
                            
            %store important values for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            %setup LQR variables
            Abar = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L);
            Abar(1:obj.nX_p,obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = ...
                eye(obj.nX_pDOT);
            Abar(obj.nX_p+1:end,:) = [dpDDOTdX;dRLdX;dLdX];
            Bbar = zeros(obj.nX_p+obj.nX_pDOT+obj.nX_RL+obj.nX_L,...
                size(obj.cableConstraintMatrix,2)+size(obj.rodConstraintMatrix,2));
            Bbar(end-(obj.nX_RL+obj.nX_L)+1:end,:) = ...
                [dRLdU;dLdU];
            xDOTBar = [Xhat.pDOT;pDDOT_bar;RLdot_bar;Ldot_bar]; %nominal state velocities
            
            %reference state
            Xref = [zeros(obj.nX_p,1);
                zeros(obj.nX_pDOT,1);
                obj.omega.X.RL0; %neutral cable lengths
                zeros(obj.nX_L,1)];
            
            z = obj.dT*(xDOTBar-Abar*[Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L]...
                -[zeros(obj.nX_p+obj.nX_pDOT,obj.nX_RL+obj.nX_L);... %block of zeros
                eye(obj.nX_RL),zeros(obj.nX_RL,obj.nX_L);...
                zeros(obj.nX_L,obj.nX_RL),eye(obj.nX_L)]*... %Bbar
                [Uhat.RLdot;Uhat.Ldot])+...
                obj.dT*Abar*Xref;
            
            Aaug = [eye(obj.nX)+obj.dT*Abar,z;zeros(1,obj.nX),1];
            Baug = [obj.dT*Bbar;zeros(1,size(Bbar,2))];
                
            %backward pass
            obj.P{obj.horizon+1} = obj.Q;
            for k = obj.horizon:-1:1
                obj.P{k} = Aaug'*obj.P{k+1}*Aaug-(Aaug'*obj.P{k+1}*Baug)/...
                    (obj.R+Baug'*obj.P{k+1}*Baug)*...
                    (Baug'*obj.P{k+1}*Aaug)+obj.Q;
            end

            %%optimal LQR feedback gain
            K = (obj.R+Baug'*obj.P{2}*Baug)\...
                (Baug'*obj.P{2}*Aaug);
            
            Xcombined = [Xhat.p;Xhat.pDOT;Xhat.RL;Xhat.L];
            U = -K*([Xcombined-Xref;1]);
            
            %remap U input vector (potentially smaller due to constraints)
            %to full cable/rod actuation
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
            
            OL_states.p = zeros(obj.nX_p,obj.horizon+1);
            OL_states.pDOT = zeros(obj.nX_pDOT,obj.horizon+1);
            OL_states.RL = zeros(obj.nX_RL,obj.horizon+1);
            OL_states.L = zeros(obj.nX_L,obj.horizon+1);
            OL_inputs.RLdot = zeros(obj.nX_RL,obj.horizon+1);
            OL_inputs.Ldot = zeros(obj.nX_L,obj.horizon+1);
            
            if(openLoopFlag)
                OL_states.p(:,1) = Xhat.p;
                OL_states.pDOT(:,1) = Xhat.pDOT;
                OL_states.RL(:,1) = Xhat.RL;
                OL_states.L(:,1) = Xhat.L;
                x_vec = [OL_states.p(:,1);OL_states.pDOT(:,1);
                    OL_states.RL(:,1);OL_states.L(:,1)];
                for i=1:obj.horizon
                    K = (obj.R+Baug'*obj.P{i+1}*Baug)\...
                        (Baug'*obj.P{i+1}*Aaug);
                    u_vec = -K*([x_vec-Xref;1]);
                    OL_inputs.RLdot(:,i) = obj.cableConstraintMatrix*...
                        u_vec(1:size(obj.cableConstraintMatrix,2));
                    OL_inputs.L(:,i).Ldot = obj.rodConstraintMatrix*...
                        u_vec(size(obj.cableConstraintMatrix,2)+1:end);
                    x_Diff = Aaug*[x_vec-Xref;1] + Baug*u_vec + [Xref;0];
                    x_vec = x_Diff(1:end-1);
                    
                    %parse vector into struct
                    idx=0;
                    OL_states.p(:,i+1) = x_vec(idx+1:idx+obj.nX_p);
                    idx=idx+obj.nX_p;
                    OL_states.pDOT(:,i+1) = x_vec(idx+1:idx+obj.nX_pDOT);
                    idx=idx+obj.nX_pDOT;
                    OL_states.RL(:,i+1) = x_vec(idx+1:idx+obj.nX_RL);
                    idx=idx+obj.nX_RL;
                    OL_states.L(:,i+1) = x_vec(idx+1:idx+obj.nX_L);
                end
                K = (obj.R+Baug'*obj.P{obj.horizon+1}*Baug)\...
                        (Baug'*obj.P{obj.horizon+1}*Aaug);
                u_vec = -K*([x_vec-Xref;1]);
                OL_inputs.RLdot(:,obj.horizon+1) = obj.cableConstraintMatrix*...
                    u_vec(1:size(obj.cableConstraintMatrix,2));
                OL_inputs.Ldot(:,obj.horizon+1) = obj.rodConstraintMatrix*...
                    u_vec(size(obj.cableConstraintMatrix,2)+1:end);
            else
                OL_states = [];
                OL_inputs = [];
            end

            costOutput =  [Xcombined;1]'*obj.P{1}*[Xcombined;1];

        end
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
        end
    end
end

