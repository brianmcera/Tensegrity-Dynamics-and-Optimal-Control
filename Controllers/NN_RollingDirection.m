classdef NN_RollingDirection < handle
    %QP_MPC_ROLLINGDIRECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        targetDestination
        horizon
        cableLinearVelocity
        cableMinLength
        cableMaxLength
        pairedCables
        rodMinLength
        rodMaxLength
        rodLinearVelocity
        actuationMode %1-cables only,2-rods only,3-both
        dT
        StaticConstraints
        optimizerObject
        NN
        
    end
    
    methods
        function obj = NN_RollingDirection(...
                X,U,omega,simTimeStep,horizon,NNfilename)
            %NN_ROLLINGDIRECTION Construct an instance of this class
            %   Detailed explanation goes here
            
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
            obj.pairedCables = omega.cables.paired;
            
            %initialize neural network
            setNNnet(NNfilename);
                    
        end
      
        
        function setActuationMode(obj,desMode)
            obj.actuationMode = desMode;
        end
        
        
        function setNNnet(obj,NNfilename)
            %instantiate NN
            load(['.\Results\Learning\TrialNeuralNetworks\',NNfilename],...
                'net2')
            obj.NN = net2;
        end
        
        
        function setTargetDestination(obj,target)
            obj.targetDestination = target;
        end
 
        
        function [U_desired,OL_states,OL_inputs,hVars,...
                cost,controllerOutputArgs] = getOptimalInput(...
                obj,Xhat,Uhat,pDDOTFcn,jacobianFcns,hFcns,debugFcns)
            
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
            
            xCOM = mean(Xhat.p(1:3:end));
            yCOM = mean(Xhat.p(2:3:end));
            zCOM = mean(Xhat.p(3:3:end));
            centroid = [xCOM;yCOM;zCOM];
            desiredDirection = obj.targetDestination-[xCOM;yCOM];
            desiredDirection = desiredDirection/norm(desiredDirection);
            zWeight = 0;
            desiredDirection = [desiredDirection;zWeight];
            
            %store important values for recordkeeping
            controllerOutputArgs.desDirection = desiredDirection;
            
            Xnodes = Xhat.p(1:3:end)-centroid(1);
            Ynodes = Xhat.p(2:3:end)-centroid(2);
            Znodes = Xhat.p(3:3:end)-centroid(3);
            
            rot_angle = atan2(desiredDirection(2),desiredDirection(1));
            rot_mat = [cos(-rot_angle) -sin(-rot_angle);
                sin(-rot_angle) cos(-rot_angle)];
            XYrotatedNodes = rot_mat*[Xnodes';Ynodes'];
            XrotatedNodes = XYrotatedNodes(1,:)';
            YrotatedNodes = XYrotatedNodes(2,:)';
            
            XYDOTrotatedNodes = rot_mat*[Xhat.pDOT(1:3:end)';
                Xhat.pDOT(2:3:end)'];
            XDOTrotatedNodes = XYDOTrotatedNodes(1,:)';
            YDOTrotatedNodes = XYDOTrotatedNodes(2,:)';
            
            input = [
                XrotatedNodes;
                YrotatedNodes;
                Znodes;
                %Xhat.pDOT;
                %                 XDOTrotatedNodes;
                %                 YDOTrotatedNodes;
                %                 Xhat.pDOT(3:3:end);
                Uhat.RL;
                %Uhat.L;
                %desiredDirection
                %[1 0 0]';
                ];
             
            OL_states.p = nan(size(Xhat.p));
            OL_states.pDOT = nan(size(Xhat.pDOT));
            OL_inputs.RL = nan(size(Uhat.RL));
            OL_inputs.L = nan(size(Uhat.L));
            cost = nan(1);
            
            U_desired.RL = predict(obj.NN,input);
            if(isrow(U_desired.RL))
                U_desired.RL = U_desired.RL';
            end
            U_desired.L = Uhat.L;
            
            %paired cable restriction
            for pairs = 1:size(obj.pairedCables,1)
                RL1change = U_desired.RL(obj.pairedCables(pairs,1))-...
                    Uhat.RL(obj.pairedCables(pairs,1));
                RL2change = U_desired.RL(obj.pairedCables(pairs,2))-...
                    Uhat.RL(obj.pairedCables(pairs,2));
                if(RL1change*RL2change>(-0.1))
                    %if cables both extend/both retract, do nothing
                    U_desired.RL(obj.pairedCables(pairs,1)) =...
                        Uhat.RL(obj.pairedCables(pairs,1));
                    U_desired.RL(obj.pairedCables(pairs,2)) =...
                        Uhat.RL(obj.pairedCables(pairs,2));
                else
                    meanChange = mean([abs(RL1change),abs(RL2change)]);
                    timeScalingFactor = 1;
                    U_desired.RL(obj.pairedCables(pairs,1))=...
                        Uhat.RL(obj.pairedCables(pairs,1))+...
                        sign(RL1change)*meanChange*timeScalingFactor;
                    U_desired.RL(obj.pairedCables(pairs,2))=...
                        Uhat.RL(obj.pairedCables(pairs,2))+...
                        sign(RL2change)*meanChange*timeScalingFactor;
                end
            end
                    
            
                      
        end
        
    end
end

