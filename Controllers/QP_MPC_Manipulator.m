classdef QP_MPC_RollingDirection < handle
    %QP_MPC_ROLLINGDIRECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %yalmip variables
        p_sdpvar
        pDOT_sdpvar
        RL_sdpvar
        L_sdpvar
        pDDOT_sdpvar
        J1_sdpvar
        J2_sdpvar
        
        targetDestination
        horizon
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
        
    end
    
    methods
        function obj = QP_MPC_RollingDirection(...
                X,U,omega,simTimeStep,horizon)
            %QP_MPC_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
            %instantiate SDPVARs for Yalmip
            obj.p_sdpvar = sdpvar(size(X.p,1),horizon+1,'full');
            obj.pDOT_sdpvar = sdpvar(size(X.pDOT,1),horizon+1,'full');
            obj.RL_sdpvar = sdpvar(size(U.RL,1),horizon+1,'full');
            obj.L_sdpvar = sdpvar(size(U.L,1),horizon+1,'full');
            obj.pDDOT_sdpvar = sdpvar(size(X.p,1),1,'full');
            obj.J1_sdpvar = sdpvar(size(X.p,1),size(X.p,1)+...
                size(X.pDOT,1),'full');
            obj.J2_sdpvar = sdpvar(size(X.p,1),size(U.RL,1)+...
                size(U.L,1),'full');
            
            %controller timestep and horizon
            obj.dT = simTimeStep;
            obj.horizon = horizon;
            
            %initialize origin as default target 
            obj.targetDestination = [0 0]';
            
            %set actuation mode to cables-only by default
            obj.actuationMode = 1;
            
            %actuator limits
            obj.cableLinearVelocity = omega.cables.linear_velocity;
            obj.rodLinearVelocity = omega.rods.linear_velocity;
            obj.cableMinLength = omega.cables.minLength;
            obj.cableMaxLength = omega.cables.maxLength;
            obj.rodMaxLength = omega.rods.maxLength;
            obj.rodMinLength = omega.rods.minLength;
            
            %% set up static constraints for Yalmip~~~~~~~~~~~~~~~~~~~~~~~~
            constr = [];
            N = horizon;
            
            %{ 
                    TODO: ROD ACTUATION NOT COMPLETED YET
            %}
  
            %actuation mode constraints 
            if obj.actuationMode == 1 %rods constrained
                constr = [constr,...
                    (obj.L_sdpvar == repmat(omega.U.L0,1,N+1)):...
                    'Constant Rod Lengths'];
            elseif obj.actuationMode == 2 %cables constrained
                constr = [constr,...
                    (obj.RL_sdpvar == repmat(omega.U.RL0,1,N+1)):...
                    'Constant Cable Lengths'];
            elseif obj.actuationMode == 3 %no additional constraints
                constr = constr;
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
                (obj.RL_sdpvar(:,2:N+1)-obj.RL_sdpvar(:,1:N)<=...
                repmat(obj.cableLinearVelocity*obj.dT,1,N)):...
                'Cable Speed Max Limit'];
            constr = [constr,...
                (obj.RL_sdpvar(:,2:N+1)-obj.RL_sdpvar(:,1:N)>=...
                repmat(-obj.cableLinearVelocity*obj.dT,1,N)):...
                'Cable Speed Min Limit'];
            constr = [constr,...
                (obj.L_sdpvar(:,2:N+1)-obj.L_sdpvar(:,1:N)<=...
                repmat(obj.rodLinearVelocity*obj.dT,1,N)):...
                'Rod Speed Max Limit'];
            constr = [constr,...
                (obj.L_sdpvar(:,2:N+1)-obj.L_sdpvar(:,1:N)>=...
                repmat(-obj.rodLinearVelocity*obj.dT,1,N)):...
                'Rod Speed Min Limit'];
            
            %constant cable velocity over MPC horizon
            for k = 2:N
                constr = [constr,...
                    (obj.RL_sdpvar(:,k)-obj.RL_sdpvar(:,k-1) == ...
                    obj.RL_sdpvar(:,k+1)-obj.RL_sdpvar(:,k)):...
                    'Constant Cable Velocity'];
            end

            %paired cables
            cable_pairs = omega.cables.paired;
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
            
            %save static constraints
            obj.StaticConstraints = constr;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
      
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
        end
        
        
        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
        
        
        function createOptimizerObject(obj,costFcnHandle)
            disp('~Creating Optimizer Object for MPC Controller~')
            
            N = obj.horizon;
            constr = obj.StaticConstraints; %copy static constraints
            
            %trapezoidal dynamic constraint on p
            constr = [constr, (obj.p_sdpvar(:,2:N+1)==...
                obj.p_sdpvar(:,1:N)+obj.dT/2*...
                (obj.pDOT_sdpvar(:,1:N)+obj.pDOT_sdpvar(:,2:N+1))):...
                'Trapezoidal p-state Constraint'];         
            
            %trapezoidal approximation for pDOT 
            X_deviation = [...
                obj.p_sdpvar(:,1:N)-repmat(obj.p_sdpvar(:,1),1,N);...
                obj.pDOT_sdpvar(:,1:N)-repmat(obj.pDOT_sdpvar(:,1),1,N)];
            Xp1_deviation = [...
                obj.p_sdpvar(:,2:N+1)-repmat(obj.p_sdpvar(:,1),1,N);...
                obj.pDOT_sdpvar(:,2:N+1)-repmat(obj.pDOT_sdpvar(:,1),1,N)];
            U_deviation = [...
                obj.RL_sdpvar(:,1:N)-repmat(obj.RL_sdpvar(:,1),1,N);
                obj.L_sdpvar(:,1:N)-repmat(obj.L_sdpvar(:,1),1,N)];
            Up1_deviation = [...
                obj.RL_sdpvar(:,2:N+1)-repmat(obj.RL_sdpvar(:,1),1,N);
                obj.L_sdpvar(:,2:N+1)-repmat(obj.L_sdpvar(:,1),1,N)];
            %dynamic constraint on pDOT SDPVARs
            constr = [constr, (obj.pDOT_sdpvar(:,2:N+1)==...
                obj.pDOT_sdpvar(:,1:N)+obj.dT/2*...
                (2*repmat(obj.pDDOT_sdpvar,1,N)+...
                obj.J1_sdpvar*(X_deviation+Xp1_deviation)+...
                obj.J2_sdpvar*(U_deviation+Up1_deviation))):...
                'Trapezoidal pDOT-state linearized constraint'];
            
            X.p = obj.p_sdpvar;
            X.pDOT = obj.pDOT_sdpvar;
            U.RL = obj.RL_sdpvar;
            U.L = obj.L_sdpvar;
            
            %temporary SDPVARs (dependent on user-defined  cost function)~~
            desDirection = sdpvar(3,1);
            costArgs.desired_direction = desDirection;
            crossAxis = sdpvar(3,1); %create temporary SDPVAR
            costArgs.crossAxis = crossAxis;
            p0 = sdpvar(size(X.p,1),1); %create temporary SDPVAR
            costArgs.p0 = p0;
            Centroid = sdpvar(3,1);
            costArgs.Centroid= Centroid;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            cost = costFcnHandle(X,U,costArgs);
            
            %input parameter arguments to optimizer object
            %last arguments in list may change depending on cost arguments
            %specific to the individual cost function, create SDPVARs above
            %and add them to the 'costArgs' struct. Also remember to add
            %the appropriate arguments to the function call in
            %'getOptimalInput' below
            parameters_in = {obj.p_sdpvar(:,1),...
                obj.pDOT_sdpvar(:,1),...
                obj.RL_sdpvar(:,1),...
                obj.L_sdpvar(:,1),...
                obj.pDDOT_sdpvar,...
                obj.J1_sdpvar,...
                obj.J2_sdpvar,...
                desDirection,...
                crossAxis,...
                Centroid}; 
            
            %optimizer object output
            solutions_out = {obj.p_sdpvar,...
                obj.pDOT_sdpvar,...
                obj.RL_sdpvar,...
                obj.L_sdpvar,...
                cost};
            
            options = sdpsettings('showprogress',1,...
                'cachesolvers',1,'warning',1,'verbose',0,...
                'solver','gurobi','usex0',1,'debug',0);

            obj.optimizerObject = optimizer(constr,cost,options,...
                parameters_in,solutions_out);
        end
        
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                cost,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,pDDOTFcn,jacobianFcns,hFcns,debugFcns)
            
            %First calculate pDDOT and linearization about Xhat,Uhat
            %helper variables
            for i = 1:length(hFcns.z)
                hVars.zbar{i} = hFcns.z{i}(Xhat);
            end
            for i = 1:length(hFcns.v)
                hVars.vbar{i} = hFcns.v{i}(Xhat);
            end
            hVars.RhatRhat = hFcns.RhatRhat;
            hVars.Chat = hFcns.Chat;
            hVars.Jbar = hFcns.J(Xhat,Uhat,hVars);
            %dynamic linearization recalculated each controller iteration
            pDDOT_bar = pDDOTFcn(Xhat,Uhat,hVars);
            JacobianMat1 = [jacobianFcns.dpDDOTdp(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdpDOT(Xhat,Uhat,hVars)];
            JacobianMat2 = [jacobianFcns.dpDDOTdRL(Xhat,Uhat,hVars),...
                jacobianFcns.dpDDOTdL(Xhat,Uhat,hVars)];
            
            xCOM = mean(Xhat.p(1:3:end));
            yCOM = mean(Xhat.p(2:3:end));
            zCOM = mean(Xhat.p(3:3:end));
            centroid = [xCOM;yCOM;zCOM];
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            zWeight = 1;
            desiredDirection = [desiredDirection;zWeight];
            crossAxis = cross(desiredDirection,[0 0 1])';
            
            %store important values for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            inputs = {Xhat.p,Xhat.pDOT,Uhat.RL,Uhat.L,pDDOT_bar,...
                JacobianMat1,JacobianMat2,desiredDirection,crossAxis,...
                centroid};
            [solutions,diagnostic,~,~,obj.optimizerObject] =...
                obj.optimizerObject(inputs);
            if (diagnostic~=0)
                disp('Error in Yalmip optimization')
                check(obj.StaticConstraints)
                if(range(range(JacobianMat1))>=1e4 ||...
                        range(range(JacobianMat2))>=1e4)
                    disp('Large matrix elements in Jacobian calculations')
                    disp('Double check dynamics functions')
                    disp('Consider decreasing "beta" in the model.m file')
                end
                error('Optimization failed in Yalmip solver')
            end
             
            OL_states.p = solutions{1};
            OL_states.pDOT = solutions{2};
            OL_inputs.RL = solutions{3};
            OL_inputs.L = solutions{4};
            cost = solutions{5};
            
            U_desired.RL = solutions{3}(:,2);
            U_desired.L = solutions{4}(:,2);
                      
            if(isnan(U_desired.RL(1)))
                warning(['Controller is not executing properly. ',...
                    'SDPVARs are being returned as NaN. Check that',...
                    ' Gurobi and Yalmip are correctly installed'])
            end
        end
        
    end
end

