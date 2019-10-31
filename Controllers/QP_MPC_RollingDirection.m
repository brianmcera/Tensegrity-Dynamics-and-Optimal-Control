classdef QP_MPC_RollingDirection < handle
    %QP_MPC_ROLLINGDIRECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %yalmip variables
        p_sdpvar
        pDOT_sdpvar
        RL_sdpvar
        L_sdpvar
        RLdot_sdpvar
        Ldot_sdpvar
%         pDDOT_sdpvar
%         J1_sdpvar
%         J2_sdpvar
        
        targetDestination
        horizon
        omega
        cableLinearVelocity
        cableMinLength
        cableMaxLength
        rodMinLength
        rodMaxLength
        rodLinearVelocity
        actuationMode %1-cables only,2-rods only,3-both
        dT
        StaticConstraints
        optimizerObject
        options
        cost
        X
        U
        
    end
    
    methods
        function obj = QP_MPC_RollingDirection(...
                X,omega,simTimeStep,horizon)
            %QP_MPC_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
            %controller timestep and horizon
            obj.dT = simTimeStep;
            obj.horizon = horizon;
            
            %initialize origin as default target 
            obj.targetDestination = [0 0]';
            
            %set actuation mode to cables-only by default
            obj.actuationMode = 1;
            
            %actuator limits
            obj.omega = omega;
            obj.cableLinearVelocity = omega.cables.linear_velocity;
            obj.rodLinearVelocity = omega.rods.linear_velocity;
            obj.cableMinLength = omega.cables.minLength;
            obj.cableMaxLength = omega.cables.maxLength;
            obj.rodMaxLength = omega.rods.maxLength;
            obj.rodMinLength = omega.rods.minLength;
            
            %instantiate state SDPVARs for Yalmip
            obj.p_sdpvar = sdpvar(size(X.p,1),horizon+1,'full');
            obj.pDOT_sdpvar = sdpvar(size(X.pDOT,1),horizon+1,'full');
            obj.RL_sdpvar = sdpvar(size(X.RL,1),horizon+1,'full');
            obj.L_sdpvar = sdpvar(size(X.L,1),horizon+1,'full');
            %instantiate input SDPVARs for Yalmip
            obj.RLdot_sdpvar = sdpvar(size(X.RL,1),horizon+1,'full');
            obj.Ldot_sdpvar = sdpvar(size(X.L,1),horizon+1,'full');
                        
            %solver options
            obj.options = sdpsettings('showprogress',0,...
                'cachesolvers',1,'warning',1,'verbose',0,...
                'solver','gurobi','usex0',0,'debug',0);
            
            %instantiate SDPVARs
            obj.initializeSDPVARs(X)
            
        end
      
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
        end
        
        
        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                costOutput,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,nominalFnc,jacobianFcns,hFcns,costFcnHandle,...
                debugFcns,openLoopFlag)
                        
            N = obj.horizon;
            constr = obj.StaticConstraints; %copy static constraints
                        
            %initial conditions for MPC step
            constr = [constr, (obj.p_sdpvar(:,1) == Xhat.p):'initial p',...
                (obj.pDOT_sdpvar(:,1) == Xhat.pDOT):'initial pDOT',...
                (obj.RL_sdpvar(:,1) == Xhat.RL):'initial RL',...
                (obj.L_sdpvar(:,1) == Xhat.L):'initial L'];
            
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
            dRLdU= [jacobianFcns.dRLdRLdot(Xhat,Uhat,hVars),...
                jacobianFcns.dRLdLdot(Xhat,Uhat,hVars)];
            dLdX = [jacobianFcns.dLdp(Xhat,Uhat,hVars),...
                jacobianFcns.dLdpDOT(Xhat,Uhat,hVars),...
                jacobianFcns.dLdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dLdL(Xhat,Uhat,hVars)];
            dLdU= [jacobianFcns.dLdRLdot(Xhat,Uhat,hVars),...
                jacobianFcns.dLdLdot(Xhat,Uhat,hVars)];
            
            xCOM = mean(Xhat.p(1:3:end));
            yCOM = mean(Xhat.p(2:3:end));
            zCOM = mean(Xhat.p(3:3:end));
            centroid = [xCOM;yCOM;zCOM];
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            %zWeight = mean(obj.omega.X.p0(3:3:end))-zCOM;
            zWeight = mean(obj.omega.X.p0(end-3:3:end))-mean(Xhat.p(end-3:3:end)); %weight only on payload Z
            desiredDirection = [desiredDirection;10*zWeight];
            desiredDirection = desiredDirection/norm(desiredDirection);
            crossAxis = cross(desiredDirection,[0 0 1])';
            
            %trapezoidal dynamic constraint on p
            constr = [constr, (obj.p_sdpvar(:,2:N+1)==...
                obj.p_sdpvar(:,1:N)+obj.dT/2*...
                (obj.pDOT_sdpvar(:,1:N)+obj.pDOT_sdpvar(:,2:N+1))):...
                'Trapezoidal p-state Constraint'];         
            
            %trapezoidal approximation for pDOT
            Xvec_deviation = [...
                obj.p_sdpvar(:,1:N+1)-repmat(obj.p_sdpvar(:,1),1,N+1);...
                obj.pDOT_sdpvar(:,1:N+1)-repmat(obj.pDOT_sdpvar(:,1),1,N+1);
                obj.RL_sdpvar(:,1:N+1)-repmat(obj.RL_sdpvar(:,1),1,N+1);
                obj.L_sdpvar(:,1:N+1)-repmat(obj.L_sdpvar(:,1),1,N+1)];
            Uvec_deviation = [...
                obj.RLdot_sdpvar(:,1:N+1)-repmat(Uhat.RLdot,1,N+1);
                obj.Ldot_sdpvar(:,1:N+1)-repmat(Uhat.Ldot,1,N+1)];
            
            %dynamic constraint on pDOT SDPVARs
            constr = [constr, (obj.pDOT_sdpvar(:,2:N+1)==...
                obj.pDOT_sdpvar(:,1:N)+obj.dT/2*...
                (2*repmat(pDDOT_bar,1,N)+...
                dpDDOTdX*(Xvec_deviation(:,1:N)+Xvec_deviation(:,2:N+1))+...
                dpDDOTdU*(Uvec_deviation(:,1:N)+Uvec_deviation(:,2:N+1)))):...
                'Trapezoidal pDOT-state linearized constraint'];
            
            %dynamic constraint on RL SDPVARs
            constr = [constr, (obj.RL_sdpvar(:,2:N+1)==...
                obj.RL_sdpvar(:,1:N)+obj.dT/2*...
                (2*repmat(RLdot_bar,1,N)+...
                dRLdX*(Xvec_deviation(:,1:N)+Xvec_deviation(:,2:N+1))+...
                dRLdU*(Uvec_deviation(:,1:N)+Uvec_deviation(:,2:N+1)))):...
                'Trapezoidal RL-state linearized constraint'];
            
            %dynamic constraint on L SDPVARs
            constr = [constr, (obj.L_sdpvar(:,2:N+1)==...
                obj.L_sdpvar(:,1:N)+obj.dT/2*...
                (2*repmat(Ldot_bar,1,N)+...
                dLdX*(Xvec_deviation(:,1:N)+Xvec_deviation(:,2:N+1))+...
                dLdU*(Uvec_deviation(:,1:N)+Uvec_deviation(:,2:N+1)))):...
                'Trapezoidal L-state linearized constraint'];
                                
            %store important values for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            runtimeArgs.desired_direction = desiredDirection;
            runtimeArgs.crossAxis = crossAxis;
            runtimeArgs.p0 = Xhat.p;
            runtimeArgs.Centroid = centroid;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            obj.cost = costFcnHandle(obj.X,obj.U,runtimeArgs);
            
            diagnostic = optimize(constr,obj.cost,obj.options);
            constr = [];
            

            if (diagnostic.problem~=0)
                diagnostic
                disp('Error in Yalmip optimization')
                check(obj.StaticConstraints)
                error('Optimization failed in Yalmip solver')
            end
            
%             OL_states = struct();
            OL_states.p = value(obj.p_sdpvar);
            OL_states.pDOT = value(obj.pDOT_sdpvar);
            OL_states.RL = value(obj.RL_sdpvar);
            OL_states.L = value(obj.L_sdpvar);
            OL_inputs.RLdot = value(obj.RLdot_sdpvar);
            OL_inputs.Ldot = value(obj.Ldot_sdpvar);
            costOutput =  value(obj.cost);
            
            U_desired.RLdot = mean(OL_inputs.RLdot(:,1:2),2);
            U_desired.Ldot = mean(OL_inputs.Ldot(:,1:2),2);
                      
            if(isnan(U_desired.RLdot(1)))
                warning(['Controller is not executing properly. ',...
                    'SDPVARs are being returned as NaN. Check that',...
                    ' Gurobi and Yalmip are correctly installed'])
            end
            
            %re-instantiate SDPVARs
            %(performance decreased over time if they were persistent)
            %(performance is better when called at end of this loop)
            obj.initializeSDPVARs(Xhat)
        end
        
        function initializeSDPVARs(obj,Xhat)
            yalmip('clear')
            N = obj.horizon;
            
            %instantiate state SDPVARs for Yalmip
            obj.p_sdpvar = sdpvar(size(Xhat.p,1),N+1,'full');
            obj.pDOT_sdpvar = sdpvar(size(Xhat.pDOT,1),N+1,'full');
            obj.RL_sdpvar = sdpvar(size(Xhat.RL,1),N+1,'full');
            obj.L_sdpvar = sdpvar(size(Xhat.L,1),N+1,'full');
            %instantiate input SDPVARs for Yalmip
            obj.RLdot_sdpvar = sdpvar(size(Xhat.RL,1),N+1,'full');
            obj.Ldot_sdpvar = sdpvar(size(Xhat.L,1),N+1,'full');
            
            %reset constraint set
            constr = [];
            
            %actuation mode constraints 
            if obj.actuationMode == 1 %rods constrained
                constr = [constr,...
                    (obj.L_sdpvar(:,2:N+1) == obj.L_sdpvar(:,1:N)):...
                    'Constant Rod Lengths'];
            elseif obj.actuationMode == 2 %cables constrained
                constr = [constr,...
                    (obj.RL_sdpvar(:,2:N+1) == obj.RL_sdpvar(:,1:N)):...
                    'Constant Cable Lengths'];
            elseif obj.actuationMode == 3 %no additional constraints
                %do nothing
            end
            
            %passive cables
            passive_cables = obj.omega.cables.passive;
            if numel(passive_cables>0)
                constr = [constr,...
                    (obj.RL_sdpvar(passive_cables,2:N+1) ==...
                    obj.RL_sdpvar(passive_cables,1:N)):...
                    'Passive Cables'];
            end
            
            %actuator limits
            constr = [constr,...
                (obj.RL_sdpvar >= repmat(obj.cableMinLength,1,N+1)):...
                'Minimum Cable Length Limit'];
            constr = [constr,...
                (obj.RL_sdpvar <= repmat(obj.cableMaxLength,1,N+1)):...
                'Maximum Cable Length Limit'];
            constr = [constr,...
                (obj.L_sdpvar >= repmat(obj.rodMinLength,1,N+1)):...
                'Minimum Rod Length Limit'];
            constr = [constr,...
                (obj.L_sdpvar <= repmat(obj.rodMaxLength,1,N+1)):...
                'Maximum Rod Length Limit'];
            
            %speed limits
            constr = [constr,...
                (obj.RLdot_sdpvar(:,1:N+1)<=...
                repmat(obj.cableLinearVelocity,1,N+1)):...
                'Cable Speed Max Limit'];
            constr = [constr,...
                (obj.RLdot_sdpvar(:,1:N+1)>=...
                repmat(-obj.cableLinearVelocity,1,N+1)):...
                'Cable Speed Min Limit'];
            constr = [constr,...
                (obj.Ldot_sdpvar(:,1:N+1)<=...
                repmat(obj.rodLinearVelocity,1,N+1)):...
                'Rod Speed Max Limit'];
            constr = [constr,...
                (obj.Ldot_sdpvar(:,1:N+1)>=...
                repmat(-obj.rodLinearVelocity,1,N+1)):...
                'Rod Speed Min Limit'];
            
            %             %constant cable velocity over MPC horizon
            %             for k = 2:N
            %                 constr = [constr,...
            %                     (obj.RL_sdpvar(:,k)-obj.RL_sdpvar(:,k-1) == ...
            %                     obj.RL_sdpvar(:,k+1)-obj.RL_sdpvar(:,k)):...
            %                     'Constant Cable Velocity'];
            %             end

            %paired cables
            cable_pairs = obj.omega.cables.paired;
            if(numel(cable_pairs)~=0)
                for p=2:N+1
                    for i = 1:size(cable_pairs,1)
                        constr = [constr,...
                            (obj.RL_sdpvar(cable_pairs(i,1),p)+...
                            obj.RL_sdpvar(cable_pairs(i,2),p) == ...
                            obj.RL_sdpvar(cable_pairs(i,1),1)+...
                            obj.RL_sdpvar(cable_pairs(i,2),1)):...
                            'Cable Pairs'];
                    end
                end
            end
            
            %equal actuated cables
            cable_similar = obj.omega.cables.similar;
            if(numel(cable_similar)~=0)
                for p=1:N+1
                    for i = 1:size(cable_similar,1)
                        constr = [constr,...
                            (obj.RLdot_sdpvar(cable_similar(i,1),p) == ...
                            obj.RLdot_sdpvar(cable_similar(i,2),p)):...
                            'Cable Similar Actuation'];
                    end
                end
            end
            
            %save static constraints
            obj.StaticConstraints = constr;
            
            %initialize cost
            obj.X.p = obj.p_sdpvar;
            obj.X.pDOT = obj.pDOT_sdpvar;
            obj.X.RL = obj.RL_sdpvar;
            obj.X.L = obj.L_sdpvar;
            obj.U.RLdot = obj.RLdot_sdpvar;
            obj.U.Ldot = obj.Ldot_sdpvar;
        end
        
    end
end

