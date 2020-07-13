classdef UnscentedKalmanFilter_IMUEncoderNodeVel_OLD2 < handle
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
        
        randomizedCableBias 
        
        v % Process noise
        w % Measurement noise
        
        % Step forward
        Current_RL
        Current_p
        
    end
    
    methods
        function obj = UnscentedKalmanFilter_IMUEncoderNodeVel_OLD2(nominalFcn,hFcns,jacobianFcns,x0,modelInfo,dT)
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
            obj.randomizedCableBias = .05*randn(24,1); %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%
            % UKF SETUP
            %%%%%%%%%%
            
            %obj.idx = 1;
            obj.nX_p    = length(x0.p);
            obj.nX_pDOT = length(x0.pDOT); 
            obj.nX_RL	= length(x0.RL); 
            obj.nX_L    = length(x0.L); 
            obj.nX      = obj.nX_p + obj.nX_pDOT + obj.nX_RL + obj.nX_L;

            %size of observation (depends on hardware sensor configuration)
            %18xOrientation (X,Y,Z); 24xCableEncoders; 36xNodeVelocity = 42-dim vector
            obj.nY = 3*size(modelInfo.R,1) + obj.nX_RL + obj.nX_pDOT + obj.nX_L; % Adjust obj.H accordingly. 
            
%             % Init recordkeeping variables
%             obj.Xhat_m_p    = zeros(obj.nX_p, obj.totalSimSteps);
%             obj.Xhat_m_pDOT = zeros(obj.nX_pDOT, obj.totalSimSteps);
            
            % Init values at k = 0
            randPercent = .001;
            obj.Xhat_m.p   	= modelInfo.X0+0.10*max(modelInfo.X0)*randn(size(modelInfo.X0));
            obj.Xhat_m.pDOT = x0.pDOT+randPercent*max(x0.pDOT)*zeros(size(x0.pDOT));
            obj.Xhat_m.RL   = x0.RL+randPercent*max(x0.RL)*randn(size(x0.RL));
            obj.Xhat_m.L  	= x0.L+randPercent*max(x0.L)*randn(size(x0.L));
            obj.Xhat_m.p(3:3:end) = obj.Xhat_m.p(3:3:end) -...
                min(obj.Xhat_m.p(3:3:end));
            obj.Xhat_p.p   	= obj.Xhat_m.p;
            obj.Xhat_p.pDOT = obj.Xhat_m.pDOT;
            obj.Xhat_p.RL   = obj.Xhat_m.RL;
            obj.Xhat_p.L  	= obj.Xhat_m.L;
            
%             obj.Pm         	= diag(abs(randPercent^2*[...
%                 max(modelInfo.X0)^2*ones(size(modelInfo.X0));...
%                 max(x0.pDOT)^2*ones(size(x0.pDOT));...
%                 max(x0.RL)^2*ones(size(x0.RL));...
%                 max(x0.L)^2*ones(size(x0.L))]));
%             obj.Pp       	= obj.Pm;

            obj.Pm = (5e-2)^2*eye(obj.nX);
            obj.Pm(1:obj.nX_p,1:obj.nX_p) = (5e-2)^2*eye(obj.nX_p);
            obj.Pp = obj.Pm;

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
            S_v = Q; % Desired variance
            A_v = sqrtm(S_v); 
            S_w = R; % Desired variance
            A_w = sqrtm(S_w);
            % Generate Random Noise
            obj.v = randn(obj.nX,1);
            obj.w = randn(obj.nY,1);
%             obj.v = A_v*obj.v;
%             obj.w = A_w*obj.w;
            
            %USER-DEFINED:parse full-state vector to simulate sensors =====
            %==============================================================
            z = zeros(obj.nY,1);
            %simulated IMU sensor per rod
            for i = 1:size(obj.modelInfo.R,1)
                IMUvec = -kron(obj.modelInfo.R(i,:),eye(3))*Y.p; 
                z(i*3-2:i*3) = IMUvec/norm(IMUvec); 
            end
            %simulated Encoder sensor per cable
            EncoderVec = Y.RL;% + obj.randomizedCableBias; %incorporate noise here <================ADDED BIAS HERE
            z(size(obj.modelInfo.R,1)*3+1:...
                size(obj.modelInfo.R,1)*3+obj.nX_RL) = EncoderVec;
            %simulated nodal velocities
            z(size(obj.modelInfo.R,1)*3+obj.nX_RL+1:...
                size(obj.modelInfo.R,1)*3+obj.nX_RL+obj.nX_pDOT) = Y.pDOT; 
            %simulated rod length measurements
            z(end-obj.nX_L+1:...
                end) = Y.L;
            
            %incorporate Random Noise here
            z = z + A_w*randn(size(z));
            
            %==============================================================
            %==============================================================
            
            
            %%%%%%%%%%%%%%%%%%%%%%
            % UKF Implementation %
            %%%%%%%%%%%%%%%%%%%%%%
            
            Xhat_m_Vec = [obj.Xhat_m.p;obj.Xhat_m.pDOT;obj.Xhat_m.RL;obj.Xhat_m.L];
            n = numel(Xhat_m_Vec);
            sigma_x_m = zeros(n,n*2);
            Lambda = n;
            matRoot = sqrtm(Lambda*obj.Pm);
            for i = 1:n
                sigma_x_m(:,i) = Xhat_m_Vec + matRoot(:,i);
                sigma_x_m(:,n+i) = Xhat_m_Vec - matRoot(:,i);
            end
            
%             %constrain node XYZ for sampled points
%             for sample = 1:size(sigma_x_m,2)
%                 constrainedP = zeros(obj.nX_p,1);
%                 for i = 1:size(obj.modelInfo.R,1)
%                     nodeP = sigma_x_m(1:obj.nX_p,sample);
%                     midPoint = kron(abs(obj.modelInfo.R(i,:))'*...
%                         abs(obj.modelInfo.R(i,:)),eye(3))*nodeP/2;
%                     dirVector = kron(obj.modelInfo.R(i,:)'*...
%                         obj.modelInfo.R(i,:),eye(3))*nodeP;
%                     origRodLength = norm(kron(obj.modelInfo.R(i,:),eye(3))*...
%                         nodeP);
%                     if(origRodLength < 0.4)
%                         warning('0 distance node behavior')
%                     end
%                     dirVector = sigma_x_m(end-obj.nX_L+i,sample)/sqrt(2)*dirVector/norm(dirVector);
%                     constrainedP = constrainedP + dirVector + midPoint;
% %                     distCheck = norm(kron(obj.modelInfo.R(i,:),eye(3))*...
% %                         constrainedP) - sigma_x_m(end-obj.nX_L+i,sample);
% %                     if(distCheck>1e-2)
% %                         distCheck
% %                         pause()
% %                     end
%                 end
%                 sigma_x_m(1:obj.nX_p,sample) = constrainedP;
%             end

            sigma_x_p = zeros(n,n*2);
            % State prediction
            parfor i = 1:size(sigma_x_p,2)
                currSigma = struct();
                %parse Xhat_m vector into struct
                idx = 1;
                currSigma.p = sigma_x_m(idx:idx-1+numel(obj.Xhat_m.p),i);
                idx = idx + numel(obj.Xhat_m.p);
                currSigma.pDOT = sigma_x_m(idx:idx-1+numel(obj.Xhat_m.pDOT),i);
                idx = idx + numel(obj.Xhat_m.pDOT);
                currSigma.RL = sigma_x_m(idx:idx-1+numel(obj.Xhat_m.RL),i);
                idx = idx + numel(obj.Xhat_m.RL);
                currSigma.L = sigma_x_m(idx:end,i);

                XOUT = stepForward(obj,U,currSigma);
                sigma_x_p(:,i) = XOUT(end,:)';
                
%                 if(min(sigma_x_p(3:3:36,i))<-0.6)
%                     disp('Forward Dynamics Error')
%                 end
                
            end
            
%             if(min(min(sigma_x_p(3:3:36,:))) < -0.6)
%                 pause()
%             end
            Xhat_p_Vec = mean(sigma_x_p,2);
            %parse Xhat_m vector into struct
            idx = 1;
            obj.Xhat_p.p = Xhat_p_Vec(idx:idx-1+numel(obj.Xhat_p.p));
            idx = idx + numel(obj.Xhat_m.p);
            obj.Xhat_p.pDOT = Xhat_p_Vec(idx:idx-1+numel(obj.Xhat_p.pDOT));
            idx = idx + numel(obj.Xhat_m.pDOT);
            obj.Xhat_p.RL = Xhat_p_Vec(idx:idx-1+numel(obj.Xhat_p.RL));
            idx = idx + numel(obj.Xhat_p.RL);
            obj.Xhat_p.L = Xhat_p_Vec(idx:end);
                        
            %calculate sample variance
            temp_Pp = Q;
            for i = 1:2*n
                temp_Pp = temp_Pp + 1/(2*n)*(sigma_x_p(:,i)-Xhat_p_Vec)*...
                    (sigma_x_p(:,i)-Xhat_p_Vec)';
            end
            obj.Pp = temp_Pp;
                        
            % MEASUREMENT UPDATE
            disp('Measurement Update')
            sigma_z = zeros(obj.nY,size(sigma_x_p,2));
            
            %% ===USER-DEFINED:expected readings based on prior estimates ==
            %==============================================================
            for sample = 1:size(sigma_z,2)
                currentX = struct();
                currentX.p = sigma_x_p(1:obj.nX_p,sample);
                currentX.pDOT = sigma_x_p(obj.nX_p+1:obj.nX_p+obj.nX_pDOT,sample);
                currentX.L = sigma_x_p(end-obj.nX_L+1:end,sample);
                for i = 1:size(obj.modelInfo.R,1)
                    IMUvec = -kron(obj.modelInfo.R(i,:),eye(3))*currentX.p;
                    sigma_z(i*3-2:i*3,sample) = IMUvec/norm(IMUvec);
                end
                currentX.RL = sigma_x_p(size(obj.modelInfo.R,2)*3*2+1:...
                    size(obj.modelInfo.R,2)*3*2+size(obj.modelInfo.C,1),sample);
                EncoderVec = currentX.RL; %incorporate noise here
                sigma_z(size(obj.modelInfo.R,1)*3+1:...
                    size(obj.modelInfo.R,1)*3+obj.nX_RL,sample) = EncoderVec;
                sigma_z(size(obj.modelInfo.R,1)*3+obj.nX_RL+1:...
                    size(obj.modelInfo.R,1)*3+obj.nX_RL+obj.nX_pDOT,...
                    sample) = currentX.pDOT;
                sigma_z(end-obj.nX_L+1:...
                    end,sample) = currentX.L;            
            end
            %==============================================================
            %==============================================================
            zHat = mean(sigma_z,2);
            Pzz = R;
            Pxz = zeros(obj.nX,obj.nY);
            for i = 1:2*n
                Pzz = Pzz + 1/(2*n)*(sigma_z(:,i)-zHat)*...
                    (sigma_z(:,i)-zHat)';
                Pxz = Pxz + 1/(2*n)*(sigma_x_p(:,i)-Xhat_p_Vec)*...
                    (sigma_z(:,i)-zHat)';
            end
 
            K = Pxz/Pzz;
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
            
            if(min(obj.Xhat_m.p(3:3:36)<-0.6))
                    disp('Measurement Update Error')
            end
                
            % Measurement Variance Update (Joseph Form)
            obj.Pm = obj.Pp-K*Pzz*K';
            
            %check known box constraints
            obj.Xhat_m.RL = min(obj.Xhat_m.RL,obj.modelInfo.cables.maxLength);
            obj.Xhat_m.RL = max(obj.Xhat_m.RL,obj.modelInfo.cables.minLength);
            
            %consistent node positions with rod length
            constrainedP = zeros(size(obj.Xhat_m.p));
            for i = 1:size(obj.modelInfo.R,1)
                midPoint = kron(abs(obj.modelInfo.R(i,:))'*...
                    abs(obj.modelInfo.R(i,:)),eye(3))*obj.Xhat_m.p/2;
                dirVector = kron(obj.modelInfo.R(i,:)'*...
                    obj.modelInfo.R(i,:),eye(3))*obj.Xhat_m.p;
                origRodLength = norm(kron(obj.modelInfo.R(i,:),eye(3))*...
                    obj.Xhat_m.p);
                if(origRodLength < 0.2)
                    warning('0 distance nodes')
                end
                dirVector = obj.Xhat_m.L(i)/sqrt(2)*dirVector/norm(dirVector);
                constrainedP = constrainedP + dirVector + midPoint;
            end
            obj.Xhat_m.p = constrainedP;
            
            %node on floor assumption
            if min(obj.Xhat_m.p(3:3:end) > 0)
                obj.Xhat_m.p(3:3:end) = obj.Xhat_m.p(3:3:end) - ...
                    min(obj.Xhat_m.p(3:3:end))-.005; %baseFloor hardcode
            end
            
%             %rodLength assumption
%             obj.Xhat_m.L = Y.L;
            
            %center node XY positions
            obj.Xhat_m.p(1:3:end) = obj.Xhat_m.p(1:3:end)-...
                mean(obj.Xhat_m.p(1:3:end))+...
                mean(Y.p(1:3:end));
            obj.Xhat_m.p(2:3:end) = obj.Xhat_m.p(2:3:end)-...
                mean(obj.Xhat_m.p(2:3:end))+...
                mean(Y.p(2:3:end));
                  
            %prepare output
%             min(obj.Xhat_p.p)
%             min(obj.Xhat_m.p)
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
        
        function XOUT = stepForward(obj,U,xHat_m)
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
            
            %{
                    **TODO: HANDLE TENSION CONSTRAINTS
            %}
            
            %forward simulate dynamics with ODE solver
            XIN = [xHat_m.p;xHat_m.pDOT;xHat_m.RL;xHat_m.L];
            [~,XOUT] = ode45(@(t,Xstacked)...
                tensegrityODE(t,Xstacked,Uinput,obj.nominalFcn,obj.hFcns),[0 obj.dT],XIN);     

        end
    end
end

