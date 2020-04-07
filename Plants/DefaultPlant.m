classdef DefaultPlant < handle
    %DEFAULTPLANT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        Current_p 
        Current_pDOT
        Current_RL
        Current_L
        cableLinearVelocity
        rodLinearVelocity
        cableMaxLength
        cableMinLength
        rodMinLength
        rodMaxLength
        simulationTime
        dT
    end
    
    methods
        function obj = DefaultPlant(x0,modelInfo,dT)
            %DEFAULTPLANT Construct an instance of this class
            %   Detailed explanation goes here
            obj.simulationTime = 0;
            obj.Current_p = x0.p;
            obj.Current_pDOT = x0.pDOT;
            obj.Current_RL = x0.RL;
            obj.Current_L = x0.L;
            obj.dT = dT;
            obj.cableLinearVelocity = modelInfo.cables.linear_velocity;
            obj.rodLinearVelocity = modelInfo.rods.linear_velocity;
            obj.cableMaxLength = modelInfo.cables.maxLength;
            obj.cableMinLength = modelInfo.cables.minLength;
            obj.rodMaxLength = modelInfo.rods.maxLength;
            obj.rodMinLength = modelInfo.rods.minLength;
        end
        
        
        function stepForward(obj,U,nominalFnc,hFcns)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            %cable update, subject to speed limits
            Uinput.RLdot = sign(U.RLdot).*...
                min(abs(U.RLdot),...
                obj.cableLinearVelocity);
            %rod update, subject to speed limits
            Uinput.Ldot = sign(U.Ldot).*...
                min(abs(U.Ldot),...
                obj.rodLinearVelocity);
            
            %cable min/max length limits
            Uinput.RLdot = max(Uinput.RLdot,(obj.cableMinLength-obj.Current_RL)/obj.dT);
            Uinput.RLdot = min(Uinput.RLdot,(obj.cableMaxLength-obj.Current_RL)/obj.dT);
            %rod min/max length limits
            Uinput.Ldot = max(Uinput.Ldot,(obj.rodMinLength-obj.Current_L)/obj.dT);
            Uinput.Ldot = min(Uinput.Ldot,(obj.rodMaxLength-obj.Current_L)/obj.dT);
            
            %forward simulate dynamics with ODE solver
            XIN = [obj.Current_p;obj.Current_pDOT;...
                obj.Current_RL;obj.Current_L];
            [~,XOUT] = ode45(@(t,Xstacked)...
                tensegrityODE(t,Xstacked,Uinput,nominalFnc,hFcns),[0 obj.dT],XIN);
            
            obj.Current_p = XOUT(end,...
                1:...
                length(obj.Current_p))';
            obj.Current_pDOT = XOUT(end,...
                length(obj.Current_p)+1:...
                length(obj.Current_p)*2)';
            obj.Current_RL = XOUT(end,...
                length(obj.Current_p)*2+1:...
                length(obj.Current_p)*2+length(obj.Current_RL))';
            obj.Current_L = XOUT(end,...
                end-length(obj.Current_L)+1:...
                end)'; 
            
            obj.simulationTime = obj.simulationTime + obj.dT;
        end
        
        
        function Y = outputSensors(obj)
            %full-state information
            Y.p = obj.Current_p;
            Y.pDOT = obj.Current_pDOT;
            Y.RL = obj.Current_RL;
            Y.L = obj.Current_L;
        end
                
        
        function KE = getKineticEnergy(obj)
            numNodes = numel(obj.Current_p)/3;
            KE = 1/2*sum(kron(eye(numNodes),[1 1 1])*...
                obj.Current_pDOT.^2);
        end
        
        
        function resetClock(obj)
            obj.simulationTime = 0;
        end
        
    end
end

