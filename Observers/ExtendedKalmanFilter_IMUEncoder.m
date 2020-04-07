classdef ExtendedKalmanFilter_IMUEncoder < handle
    %DEFAULTOBSERVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nominalFcn
        Xhat_m
        Xhat_p
        
%         Xhat_m_p
%         Xhat_m_pDOT
%         Xhat_m_RL
%         Xhat_m_L
%         Xhat_p_p
%         Xhat_p_pDOT
%         Xhat_p_RL
%         Xhat_p_L
        
        hFcns
        jacobianFcns
        Pp % Prior state variance
        Pm % Posterior Measurement state variance
%         Q % Process Variance
%         R % Measurement Variance
%         H % Linearized Measurement Model about state; y = h(x) = H*x
%         M %Linearized Measurement Model about Noise
%         L %Linearized Process Model about Noise
        
        dT 
        
        modelInfo
        
        nX_p
        nX_pDOT
        nX_RL
        nX_L
        nX
        
        %idx
        totalSimSteps
        nY
        
        v % Process noise
        w % Measurement noise
        
        % Step forward
        Current_RL
        Current_p
        
    end
    
    methods
        function obj = ExtendedKalmanFilter_IMUEncoder(nominalFcn,hFcns,jacobianFcns,x0,modelInfo,dT)
            % pDDOTFcn calculates double derivative, don't need it
            % hFcns helper functions. 
            % Add a state estimate that is kept track in KF
            % add jacobian functions, pDDOT dynamics function, etc. into
            % object properties
            
            %dynamics function handles used for the forward simulation
            obj.nominalFcn = nominalFcn;
            obj.hFcns = hFcns;
            obj.jacobianFcns = jacobianFcns;
            
            obj.dT = dT;
            obj.modelInfo = modelInfo;
            
            % Generate Noise
            % Init noise with randn
            rng('shuffle')
            
            %%%%%%%%%%
            % EKF SETUP
            %%%%%%%%%%
            
            %obj.idx = 1;
            obj.nX_p    = length(x0.p);
            obj.nX_pDOT = length(x0.pDOT); 
            obj.nX_RL	= length(x0.RL); 
            obj.nX_L    = length(x0.L); 
            obj.nX      = obj.nX_p + obj.nX_pDOT + obj.nX_RL + obj.nX_L;

            %size of observation (depends on hardware sensor configuration)
            %6xOrientation (X,Y,Z); 24xCableEncoders = 42-dim vector
            obj.nY = 3*size(modelInfo.R,1) + size(modelInfo.C,1) + obj.nX_pDOT; % Adjust obj.H accordingly. 
            
%             % Init recordkeeping variables
%             obj.Xhat_m_p    = zeros(obj.nX_p, obj.totalSimSteps);
%             obj.Xhat_m_pDOT = zeros(obj.nX_pDOT, obj.totalSimSteps);
%             obj.Xhat_p_p    = zeros(obj.nX_p, obj.totalSimSteps);
%             obj.Xhat_p_pDOT = zeros(obj.nX_pDOT, obj.totalSimSteps);
%             obj.Xhat_p_RL   = zeros(obj.nX_RL, obj.totalSimSteps);
%             obj.Xhat_m_RL   = zeros(obj.nX_RL, obj.totalSimSteps);
%             obj.Xhat_p_L    = zeros(obj.nX_RL, obj.totalSimSteps);
%             obj.Xhat_m_L    = zeros(obj.nX_RL, obj.totalSimSteps);
%             obj.Pm          = zeros(obj.nX, obj.nX, obj.totalSimSteps);
%             obj.Pp          = zeros(obj.nX, obj.nX, obj.totalSimSteps);
            
            % Init values at k = 0
            obj.Xhat_m.p   	= x0.p;
            obj.Xhat_m.pDOT = x0.pDOT;
            obj.Xhat_m.RL   = x0.RL;
            obj.Xhat_m.L  	= x0.L;
            obj.Xhat_p.p   	= x0.p;
            obj.Xhat_p.pDOT = x0.pDOT;
            obj.Xhat_p.RL   = x0.RL;
            obj.Xhat_p.L  	= x0.L;

            obj.Pm         	= (1e-6)*eye(obj.nX);
            obj.Pp       	= obj.Pm;

            %%%%%%%%%%%% 
            %%%%%%%%%%%%
            
        end 
        
        function [Xhat_m,z,PmVar,PpVar] = estimateState(obj,Y,U,Q,R)
            
            % ESTIMATESTATE This function takes in the full-state
            % observations from 'defaultPlant' and parses the struct for the
            % necessary Xhat,Uhat output for the controller
            % find hvars the same way the controller does 
            % hvars is used to find jacobians
            % Might not need to use gamma, genforces, constrforces for kf 
            
            %Recording down current observations as a property may
            %be useful for observers that require updating based on past
            %estimates (e.g., Kalman Filter)
            
            
            % Set Noise variance
            S_v = 0.001*eye(obj.nX); % Desired variance
            A_v = sqrtm(S_v); 
            S_w = 1e-6*eye(obj.nY); % Desired variance
            A_w = sqrtm(S_w);
            % Generate Random Noise
            obj.v = randn(obj.nX,1);
            obj.w = randn(obj.nY,1);
            obj.v = A_v*obj.v;
            obj.w = A_w*obj.w;
            
            %USER-DEFINED:parse full-state vector to simulate sensors =====
            %==============================================================
            
            z = zeros(obj.nY,1);
            for i = 1:size(obj.modelInfo.R,1)
                %simulated IMU sensor per rod
                IMUvec = -kron(obj.modelInfo.R(i,:),eye(3))*Y.p; 
                z(i*3-2:i*3) = IMUvec/norm(IMUvec); 
            end
            %simulated Encoder sensor per cable
            EncoderVec = Y.RL; %incorporate noise here
            z(size(obj.modelInfo.R,1)*3+1:size(obj.modelInfo.R,1)*3+obj.nX_RL) = EncoderVec;
            
            z(size(obj.modelInfo.R,1)*3+obj.nX_RL+1:size(obj.modelInfo.R,1)*3+obj.nX_RL+obj.nX_pDOT) = Y.pDOT;
            
            %incorporate Random Noise here

            %==============================================================
            %==============================================================
            
            
            %%%%%%%%%%%%%%%%%%%%%
            % EKF Implementation 
            %%%%%%%%%%%%%%%%%%%%%
            
            %First calculate pDDOT and linearization about Xhat,Uhat
            % Jacobian is calculated about Xhat_m(k-1)
            for i = 1:length(obj.hFcns.z)
                hVars.z{i} = obj.hFcns.z{i}(obj.Xhat_m);
            end
            for i = 1:length(obj.hFcns.v)
                hVars.v{i} = obj.hFcns.v{i}(obj.Xhat_m);
            end
            hVars.RhatRhat = obj.hFcns.RhatRhat;
            hVars.Chat = obj.hFcns.Chat;
            hVars.J = obj.hFcns.J(obj.Xhat_m, U, hVars);
            % Find jacbians 
            jacobians.dpDDOTdp =  obj.jacobianFcns.dpDDOTdp(obj.Xhat_m,U,hVars);
            jacobians.dpDDOTdpDOT = obj.jacobianFcns.dpDDOTdpDOT(obj.Xhat_m,U,hVars);
            jacobians.dpDDOTdRL = obj.jacobianFcns.dpDDOTdRL(obj.Xhat_m,U,hVars);
            jacobians.dpDDOTdL = obj.jacobianFcns.dpDDOTdL(obj.Xhat_m,U,hVars);
            
            %Linearized Process Matrices
            A = [zeros(obj.nX_p,obj.nX_p), eye(obj.nX_p,obj.nX_pDOT),zeros(obj.nX_p,obj.nX_RL),zeros(obj.nX_p,obj.nX_L);
                jacobians.dpDDOTdp, jacobians.dpDDOTdpDOT,jacobians.dpDDOTdRL,jacobians.dpDDOTdL;
                zeros(obj.nX_RL,obj.nX);
                zeros(obj.nX_L,obj.nX)]; 
            L = eye(obj.nX);
            
            % State prediction
            [xout] = stepForward(obj,U,obj.Xhat_m);
            obj.Xhat_p.p = xout.p;
            obj.Xhat_p.pDOT = xout.pDOT;
            obj.Xhat_p.RL = xout.RL;
            obj.Xhat_p.L = xout.L;
            
            % Prior Variance update
            obj.Pp = A*obj.Pm*A' + L*Q*L';
            
            
            % MEASUREMENT UPDATE
            
            %Linearized (about current estimate)Measurement Matrices
            H = zeros(obj.nY,obj.nX);
            for i = 1:size(obj.modelInfo.R,1)
                kronMat = -kron(obj.modelInfo.R(i,:),eye(3));
                s = kronMat*obj.Xhat_p.p;
                H(i*3-2:i*3,1:obj.nX_p) = ...
                    1/norm(s)*(eye(3)-s*s'/(s'*s))*kronMat;     
            end
            H(size(obj.modelInfo.R,1)*3+1:size(obj.modelInfo.R,1)*3+obj.nX_RL,...
                obj.nX_p+obj.nX_pDOT+1:obj.nX_p+obj.nX_pDOT+obj.nX_RL) = ...
                eye(obj.nX_RL);
            H(size(obj.modelInfo.R,1)*3+obj.nX_RL+1:size(obj.modelInfo.R,1)*3+obj.nX_RL+obj.nX_pDOT,...
                obj.nX_p+1:obj.nX_p+obj.nX_pDOT) = eye(obj.nX_pDOT);
            M = eye(obj.nY);
            
            %USER-DEFINED:expected readings based on prior estimates ======
            %==============================================================
            zHat = zeros(obj.nY,1);
            for i = 1:size(obj.modelInfo.R,1)
                IMUvec = -kron(obj.modelInfo.R(i,:),eye(3))*obj.Xhat_p.p; 
                zHat(i*3-2:i*3) = IMUvec/norm(IMUvec); 
            end
            EncoderVec = obj.Xhat_p.RL; %incorporate noise here
            zHat(size(obj.modelInfo.R,1)*3+1:size(obj.modelInfo.R,1)*3+obj.nX_RL) = EncoderVec;
            zHat(size(obj.modelInfo.R,1)*3+obj.nX_RL+1:size(obj.modelInfo.R,1)*3+obj.nX_RL+obj.nX_pDOT) = Y.pDOT;
            %==============================================================
            %==============================================================
            
            K = obj.Pp*H' / (H*obj.Pp*H' + M*R*M');
            Xhat_p_Vec = [obj.Xhat_p.p;obj.Xhat_p.pDOT;obj.Xhat_p.RL;obj.Xhat_p.L];
            Xhat_m_Vec = Xhat_p_Vec + K*(z - zHat);
            %parse Xhat_m vector into struct
            idx = 1;
            obj.Xhat_m.p = Xhat_m_Vec(idx:idx-1+numel(obj.Xhat_m.p));
            idx = idx + numel(obj.Xhat_m.p);
            obj.Xhat_m.pDOT = Xhat_m_Vec(idx:idx-1+numel(obj.Xhat_m.pDOT));
            idx = idx + numel(obj.Xhat_m.pDOT);
            obj.Xhat_m.RL = Xhat_m_Vec(idx:idx-1+numel(obj.Xhat_m.RL));
            idx = idx + numel(obj.Xhat_m.RL);
            obj.Xhat_m.L = Xhat_m_Vec(idx:end);
            
            
            % Measurement Variance Update (Joseph Form)
            obj.Pm= (eye(size(K,1))-K*H)*...
                obj.Pp*(eye(size(K,1))-K*H)'+...
                K*R*K';
            
            %check known box constraints
            obj.Xhat_m.RL = min(obj.Xhat_m.RL,obj.modelInfo.cables.maxLength);
            obj.Xhat_m.RL = max(obj.Xhat_m.RL,obj.modelInfo.cables.minLength);
            
            %maintain known rod length
            constrainedP = zeros(size(obj.Xhat_m.p));
            for i = 1:size(obj.modelInfo.R,1)
                midPoint = kron(abs(obj.modelInfo.R(i,:))'*...
                    abs(obj.modelInfo.R(i,:)),eye(3))*obj.Xhat_m.p/2;
                dirVector = kron(obj.modelInfo.R(i,:)'*...
                    obj.modelInfo.R(i,:),eye(3))*obj.Xhat_m.p;
                origRodLength = norm(kron(obj.modelInfo.R(i,:),eye(3))*...
                    obj.Xhat_m.p);
                dirVector = obj.Xhat_m.L(i)/sqrt(2)*dirVector/norm(dirVector);
                constrainedP = constrainedP + dirVector + midPoint;
            end
            obj.Xhat_m.p = constrainedP;
            
%             %set nodes on ground
%             obj.Xhat_m.p(3:3:end) = obj.Xhat_m.p(3:3:end) -...
%                 min(obj.Xhat_m.p(3:3:end)) + (-0.40);
            
            %center node XY positions
            obj.Xhat_m.p(1:3:end) = obj.Xhat_m.p(1:3:end)-...
                mean(obj.Xhat_m.p(1:3:end))+...
                mean(Y.p(1:3:end));
            obj.Xhat_m.p(2:3:end) = obj.Xhat_m.p(2:3:end)-...
                mean(obj.Xhat_m.p(2:3:end))+...
                mean(Y.p(2:3:end));
            
            %prepare output
            min(obj.Xhat_p.p)
            min(obj.Xhat_m.p)
            Xhat_m = obj.Xhat_m;
            PmVar = diag(obj.Pm);
            PpVar = diag(obj.Pp);
            %CAN ADD STATE LIMITS HERE

            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
           
        end
        
        function [variance] = getVariance(obj)
            variance = obj.Pm;
        end
        
        function xout = stepForward(obj,U,xHat_m)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            %cable update, subject to speed limits
            Uinput.RLdot = sign(U.RLdot).*...
                min(abs(U.RLdot),...
                obj.modelInfo.cables.linear_velocity);
            %rod update, subject to speed limits
            Uinput.Ldot = sign(U.Ldot).*...
                min(abs(U.Ldot),...
                obj.modelInfo.rods.linear_velocity);
            
            %cable min/max length limits
            Uinput.RLdot = max(Uinput.RLdot,...
                (obj.modelInfo.cables.minLength-xHat_m.RL)/obj.dT);
            Uinput.RLdot = min(Uinput.RLdot,...
                (obj.modelInfo.cables.maxLength-xHat_m.RL)/obj.dT);
            %rod min/max length limits
            Uinput.Ldot = max(Uinput.Ldot,...
                (obj.modelInfo.rods.minLength-xHat_m.L)/obj.dT);
            Uinput.Ldot = min(Uinput.Ldot,...
                (obj.modelInfo.rods.maxLength-xHat_m.L)/obj.dT);
            
            %forward simulate dynamics with ODE solver
            XIN = [xHat_m.p;xHat_m.pDOT;xHat_m.RL;xHat_m.L];
            [~,XOUT] = ode45(@(t,Xstacked)...
                tensegrityODE(t,Xstacked,Uinput,obj.nominalFcn,obj.hFcns),[0 obj.dT],XIN);     
            
            xout.p = XOUT(end,1:length(obj.Xhat_p.p))';
            idx = length(obj.Xhat_p.p);
            xout.pDOT = XOUT(end,idx+1:idx+length(obj.Xhat_p.pDOT))';
            idx = idx + length(obj.Xhat_p.pDOT);
            xout.RL = XOUT(end,idx+1:idx+length(obj.Xhat_p.RL))';
            idx = idx + length(obj.Xhat_p.RL);
            xout.L = XOUT(end,idx+1:idx+length(obj.Xhat_p.L))';
        end
    end
end

