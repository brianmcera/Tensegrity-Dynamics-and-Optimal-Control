function TensegritySimulation()
%TENSEGRITYSIMULATION Summary of this function goes here
%   Detailed explanation goes here

%self notes:
%changed the model file to reflect 350 stiffness, 25 pretension
%change to slow motors (5 cm/s)
%1-8-19: changed LSTM_controller to use min paired-cable actuation

close all
clear all
clc
rng('shuffle') 
 
%% User-defined Parameters ////////////////////////////////////////////////

%model filename and custom input arguments
modelName = 'I_six_bar_model';
%plant filename and custom input arguments
plantName = 'DefaultPlant';
%observer filename and custom input arguments
observerName = 'DefaultObserver';%'UnscentedKalmanFilter_IMUEncoderNodeVel_OLD2';
%controller filename and custom input arguments
controllerName = 'iLQRminimax_RollingDirection_2';
%cost filename and custom input arguments
costName = 'velocityCost'; %specify costFunction parameters below

simTimeStep = 3e-3; %simulation timestep
totalSimSteps = 5000; %total number of simulation timesteps
controllerHorizon = 20;%5; %MPC horizon for controller
actuationMode = 1; %1-cables, 2-rods, 3-both 
%todo: rod actuation does not work right now

%dynamics generation parameters
saveDynamicsFile = false;
OptimizeFlag = false;
optimizedDynamics_Filename = 'optimizedDynamics_SixBar_19_08_28_15_53_23';

%save data
doLog =  true;
savePeriod = 10; %how often to save data (e.g., every 50 timesteps)

%forward simulation / initial conditions
show_initialization_graphics = false; %true/false
KineticEnergyDamping = false; %true/false
perturbCablesPercentage = 0.00; %double between 0 and 1
perturbRodsPercentage = 0.00; %double between 0 and 1
XYrandomRotate = false; %true/false
YZrandomRotate = false; %true/false
cameraViewInitial = [90 0]; %[az,el] or {1,2,3}

%simulation loop
show_simulation_graphics = false;
openLoopFlag = true; %toggle open-loop calculation in loop (takes longer)
cameraViewLoop = [0,15]; %syntax: either A)[az,el] or B)integer,{1,2,3}


%//////////////////////////////////////////////////////////////////////////

%% Generate Model and Dynamics
timeStamp = datestr(now,'yy_mm_dd_HH_MM_SS');
[omega,X] = feval(modelName);
constraints = omega.constraints;
generalForces = omega.generalForces;
structurePlot(X,omega,constraints,cameraViewInitial,1,0,1);
pause(1e-3)
[nominalFcn,hFcns,jacobianFcns,Gamma,genForces,constrForces,debugFcns] = ...
    Dynamics_Generator(omega,constraints,generalForces);

%build symbolic matlabFunction which optimizes computation time
%longer initialization, faster iteration loops
if(isfile(['./Dynamics/Generated Dynamics/',optimizedDynamics_Filename,'.m']))
    disp(['Existing dynamics function found: <'...
        './Dynamics/Generated Dynamics/',optimizedDynamics_Filename,'.m>'])
    disp('Using pre-generated dynamics file...')
    nominalFcn.pDDOT = @(x1,x2,x3,x4)feval(optimizedDynamics_Filename,...
        x1,x2,x3,x4);
else
    disp('Existing dynamics file NOT found, generating dynamics functions...')
    %symbolic State structure
    Xsym = struct();
    fields = fieldnames(X);
    
    %parse through all struct fields
    for fieldIdx = 1:numel(fields)
        Xsym.(fields{fieldIdx}) = ...
            sym(fields{fieldIdx},size(X.(fields{fieldIdx})),'real');
    end
    
    %symbolic Input structure
    Usym = struct();
    Usym.RLdot = sym('RLdot',size(omega.cableConstraintMatrix,2),'real');
    Usym.Ldot = sym('Ldot',size(X.L),'real');

    %symbolic Helper Variables structure
    hVarSym = struct();
    fields = fieldnames(hFcns);
    for fieldIdx = 1:numel(fields)-1
        if(~isa(hFcns.(fields{fieldIdx}){fieldIdx},'function_handle'))
            %doesn't take Xsym as an argument (element is not function handle)
            for i = 1:numel(hFcns.(fields{fieldIdx}))
                hVarSym.(fields{fieldIdx}){i} = ...
                    hFcns.(fields{fieldIdx}){i};
            end
        else
            %takes Xsym as an argument (element is function handle)
            for i = 1:numel(hFcns.(fields{fieldIdx}))
                hVarSym.(fields{fieldIdx}){i} = ...
                    hFcns.(fields{fieldIdx}){i}(Xsym);
            end
        end
    end
    hVarSym.J = hFcns.J(Xsym,Usym,hVarSym);
    disp('Generating Dynamic Function Handles...')
    if(OptimizeFlag)
        disp('Optimizing Dynamics This will take >15 minutes (optimizeFlag variable is set to "true")...')
    end
    if(saveDynamicsFile)
        nominalFcn.pDDOT = matlabFunction(nominalFcn.pDDOT(Xsym,Usym,hVarSym),...
            'File',['Dynamics/Generated Dynamics/',optimizedDynamics_Filename,...
            '_',num2str(timeStamp)],...
            'Optimize',OptimizeFlag,'Vars',{Xsym.p,Xsym.pDOT,Xsym.RL,Xsym.L});
    else
        nominalFcn.pDDOT = matlabFunction(nominalFcn.pDDOT(Xsym,Usym,hVarSym),...
            'Optimize',OptimizeFlag,'Vars',{Xsym.p,Xsym.pDOT,Xsym.RL,Xsym.L});
    end
end

%store dynamics into 'omega' struct
omega.nominalFcn = nominalFcn;
omega.hFcns = hFcns;
omega.jacobianFcns = jacobianFcns;


%% Set up desired initial conditions
%calculate initial cable/rod lengths for forward simulation

%helper function handles
z = @(p,V) -kron(V'*V,eye(3))*p; %tension direction in p-basis
separationDist = @(p,V) sqrt(z(p,V)'*z(p,V)/2); %euclidean separation dist.

%initialize cable restlengths
C = omega.C;
for cable = 1:size(C,1) %apply pretension
    X_desired.RL(cable) = separationDist(X.p,C(cable,:))-...
        omega.cables.pretension(cable)/omega.cables.stiffness(cable);
end
X_desired.RL = X_desired.RL'; %column vector

%initialize rod lengths
R = omega.R;
for rod = 1:size(R,1)
    X_desired.L(rod) = separationDist(X.p,R(rod,:));
end
X_desired.L = X_desired.L'; %column vector

%store neutral inputs of equally pretensioned robot
omega.X.RL0 = X_desired.RL;
omega.X.L0 = X_desired.L;

%perturb cable lengths around nominal pretensioned lengths
if(perturbCablesPercentage~=0)
    disp(['~~"perturbCablesPercentage" parameter not equal to nominal value of "0":',...
        ' Perturbing nominal pretensioned cable restlengths~~'])
    X_desired.RL = (eye(numel(X_desired.RL))+...
        diag(perturbCablesPercentage*rand(numel(X_desired.RL),1)-...
        perturbCablesPercentage/2))*X_desired.RL;
    %handle constrained cables
    if(~isempty(omega.cables.passive))
        X_desired.RL(omega.cables.passive) = omega.X.RL0(...
            omega.cables.passive);
    end
    %handle paired cables
    if(~isempty(omega.cables.paired))
        X_desired.RL(omega.cables.paired(:,2)) = ...
            -X_desired.RL(omega.cables.paired(:,1))+...
            omega.X.RL0(omega.cables.paired(:,1))+...
            omega.X.RL0(omega.cables.paired(:,2));
    end
    %handle similarly actuated cables
    if(~isempty(omega.cables.similar))
        X_desired.RL(omega.cables.similar(:,2)) = ...
            X_desired.RL(omega.cables.similar(:,1));
    end
end
if(perturbRodsPercentage~=0)
    X_desired.L = (eye(numel(X_desired.L))+...
        diag(perturbCablesPercentage*rand(numel(X_desired.L),1)-...
        perturbCablesPercentage/2))*X_desired.L;
    %constrained rods
    if(~isempty(omega.rods.constrained))
        X_desired.L(omega.rods.constrained) = omega.X.L0(...
            omega.rods.constrained);
    end
end

%rotate robot about vertical Z-axis
if(XYrandomRotate==true)
    disp(['~~"XYrandomRotate" parameter set to TRUE: ',...
        'Randomly rotating nodal positions about vertical Z-axis~~'])
    %center XY
    centeredX = X.p(1:3:end)-mean(X.p(1:3:end));
    centeredY = X.p(2:3:end)-mean(X.p(2:3:end));
    randRotAngle = rand*(2*pi);
    rotatedXY = [cos(randRotAngle),-sin(randRotAngle);...
        sin(randRotAngle),cos(randRotAngle)]*[centeredX';centeredY'];
    X.p(1:3:end) = rotatedXY(1,:); %rotated X values
    X.p(2:3:end) = rotatedXY(2,:); %rotated Y values
    structurePlot(X,omega,constraints,cameraViewInitial,1,0,1);
end

%rotate robot about X-axis
if(YZrandomRotate==true)
    disp('~~"YZrandomRotate" parameter set to TRUE: Randomly rotating nodal positions about horizontal X-axis~~')
    %center YZ
    minZ = min(X.p(3:3:end));
    centeredY = X.p(2:3:end)-mean(X.p(2:3:end));
    centeredZ = X.p(3:3:end)-mean(X.p(3:3:end));
    randRotAngle = rand*(2*pi);
    rotatedYZ = [cos(randRotAngle),-sin(randRotAngle);...
        sin(randRotAngle),cos(randRotAngle)]*[centeredY';centeredZ'];
    X.p(2:3:end) = rotatedYZ(1,:); %rotated X values
    X.p(3:3:end) = rotatedYZ(2,:); %rotated Y values
    X.p(3:3:end) = X.p(3:3:end)-min(X.p(3:3:end)) + minZ;
    structurePlot(X,omega,constraints,cameraViewInitial,1,0,0);
end

pause(1e-3) %give time to update and redraw structurePlot


%% Instantiate Cost function handle & Plant,Controller,Observer objects
%Handle and Objects created here after model instantiation
% Note: 'omega' parameter here reflects updated cable/rod lengths after any
%modifications to nominal model file (e.g., cable perturbations, model
%rotations, etc.

%% cost function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This function encapsulates all of the necessary parameters which can be 
%tuned to produce a user-define desired behavior from the controller.
%cost ar

%user-defined arguments to cost function
costArgs.RL_diff_weight = 5;
costArgs.RL_actuation_weight = 0;
costArgs.L_diff_weight = 0;
costArgs.velocity_reward = 5;
costArgs.omega = omega; 
costArgs.stepDiscount = 0.95;
costArgs.N = controllerHorizon;
costFcnHandle = @(X,U,runtimeArgs)feval(costName,X,U,costArgs,runtimeArgs);
disp('Generating Cost Function Handle...')

%% plant object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This object is the simulated object which represents the actual robot.
%Specifically, this object keeps track of the robot state throughout the
%simulation. This is also where any actuator input noise and sensor output
%noise should be incorporated.

plant = feval(plantName,X,omega,simTimeStep);

%% Forward Simulate to Equilibrium (with Kinetic Energy Damping)
forwardSimIter = 1;
kineticEnergy_prev = 0;
cablesStillMoving = 1; 

if(any(X_desired.RL<omega.cables.minLength))
    error(['Desired restlength is less than minimum cable ',...
        'length limits. Increase stiffness or decrease desired pretension'])
end

disp('Starting forward simulation to equilibrium...')
if(KineticEnergyDamping)
    disp('Kinetic Damping, Current KE:')
end
while((plant.getKineticEnergy()>1e-4)||... %non-negligible kinetic energy
        plant.getKineticEnergy()==0 ||... %KE just got reset
        cablesStillMoving) %cables haven't yet reached goal rest lengths
    
    %calculate appropriate cable length / rod length
    U.RLdot = sign(X_desired.RL - plant.Current_RL).*...
        min(omega.cables.linear_velocity,...
        abs(X_desired.RL - plant.Current_RL)/simTimeStep); 
    U.Ldot = zeros(size(plant.Current_L));
    
    %advance simulation one timestep
    plant.stepForward(U,nominalFcn,hFcns); 
    X.p = plant.Current_p;
    X.pDOT = plant.Current_pDOT;
    
    %plot update
    if(show_initialization_graphics==true)
        structurePlot(X,omega,constraints,cameraViewInitial,1,0,1);
        title(['t = ',num2str(plant.simulationTime)])
        drawnow()
    end
    
    %kinetic energy damping
    if(KineticEnergyDamping)
        if(plant.getKineticEnergy() < kineticEnergy_prev &&...
                forwardSimIter>3)
            plant.Current_pDOT = zeros(size(plant.Current_pDOT));
            kineticEnergy_prev = 0;
            forwardSimIter = 0;
            disp('KE Peak, wiping velocities')
        else
            kineticEnergy_prev = plant.getKineticEnergy();
            fprintf('%3.2e..',kineticEnergy_prev);
        end
    end
    
    %cables reach desired initial (perturbed) state
    if(norm(plant.Current_RL-X_desired.RL)<1e-6)
        cablesStillMoving = 0;
    end
    
    forwardSimIter = forwardSimIter+1;
end
plant.resetClock(); %reset simulation time for recordkeeping
omega.X.p0 = plant.Current_p; %record equilibrium state for reference


%% controller object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This object handles optimal control of the robot. A variety of approaches
%and algorithms can be plug-and-play inserted here. Examples include linear
%time-varying MPC, iLQR, neural network, etc. 

controllerTargets = {[30,0]'};
targetIdx = 1;
controller = feval(controllerName,X,omega,simTimeStep,controllerHorizon);
controller.setActuationMode(actuationMode);
controller.setTargetDestination(controllerTargets{targetIdx});

%% observer object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This object is the observer/estimator, which takes in the full information
%from the Plant object, simulates partially observable sensor data, and 
%outputs the best estimate of the current state.
%Additionally, any information about uncertainty/variance should be handled
%by this object.

observer = feval(observerName,nominalFcn,hFcns,jacobianFcns,...
    X,omega,simTimeStep);

%% Instantiate Record Arrays for Data Storage
%**include 'record' in filename for automated storage**
simulationParameters_record.omega = omega;
simulationParameters_record.model = modelName;
simulationParameters_record.plant = plantName;
simulationParameters_record.controller = controllerName;
simulationParameters_record.observer = observerName;
simulationParameters_record.costName = costName;
simulationParameters_record.costArgs = costArgs;
simulationParameters_record.timestep = simTimeStep;
simulationParameters_record.controllerHorizon = controllerHorizon;

%dimension of concatenated state vector
Xdim = numel(X.p)+numel(X.pDOT)+numel(X.RL)+numel(X.L);
%dimension of sensor observation vector
Ydim = 84;

X_record.p = zeros(size(X.p,1),totalSimSteps);
X_record.pDOT = zeros(size(X.pDOT,1),totalSimSteps);
X_record.RL = zeros(size(X.RL,1),totalSimSteps);
X_record.L = zeros(size(X.L,1),totalSimSteps);
U_record.RLdot = zeros(size(X.RL,1),totalSimSteps);
U_record.Ldot = zeros(size(X.L,1),totalSimSteps);
Xhat_record.p = zeros(size(X.p,1),totalSimSteps);
Xhat_record.pDOT = zeros(size(X.pDOT,1),totalSimSteps);
Xhat_record.RL = zeros(size(X.RL,1),totalSimSteps);
Xhat_record.L = zeros(size(X.L,1),totalSimSteps);
Z_record = zeros(Ydim,totalSimSteps);
Pm_record = zeros(Xdim,totalSimSteps);
Pp_record = zeros(Xdim,totalSimSteps);
X_openLoop_record = cell(1,totalSimSteps);
U_openLoop_record = cell(1,totalSimSteps);
Gamma_record = zeros(size(X.p,1),totalSimSteps);
GenForce_record = zeros(size(X.p,1),totalSimSteps);
ConstrForce_record = zeros(size(X.p,1),totalSimSteps);
cost_record = zeros(1,totalSimSteps);
controllerOutputArgs_record = cell(1,totalSimSteps);


%% Simulation Loop

LoopStart = tic;
resetCables = 0;

%initial zero input for timestep t=0;
U.RLdot = zeros(size(X.RL,1),1);
U.Ldot = zeros(size(X.L,1),1);

for iteration = 1:totalSimSteps
    disp('')
    disp(['Simulation Loop Timestamp: ',timeStamp])
    disp(['Iteration: ',num2str(iteration)])
    disp(['Total Elapsed Time: ',num2str(toc(LoopStart))])
    
    %% Simulate Plant (takes U, outputs Y)
    %handle cable actuation, motor dynamics, input saturation, input/output
    %noise
    plantElapsed = tic;
    plant.stepForward(U,nominalFcn,hFcns);
    Y = plant.outputSensors();    
    %plot results
    if(show_simulation_graphics)
        X.p = plant.Current_p;
        X.pDOT = plant.Current_pDOT;
        structurePlot(X,omega,constraints,cameraViewLoop,1,0,1);
        title(['t = ',num2str(plant.simulationTime)])
    end
    disp(['Simulate Plant Elapsed Time: ',num2str(toc(plantElapsed))])
    
    %% Observer (takes Y, outputs Xhat)
    %examples:full-state information, Luenberg Observer, Kalman Filter,
    %Extended Kalman Filter, Unscented Kalman Filter
    observerElapsed = tic;
    Qvar = (1e-2)^2*eye(size(X.p,1)+size(X.pDOT,1)+size(X.RL,1)+size(X.L,1));
    Rvar = diag([(1e-2)^2*ones(18,1);
        (10e-2)^2*ones(24,1);
        (1e-2)^2*ones(36,1);
        (1e-6)^2*ones(6,1)]);%observer vector dimension = 84
    %     Rvar = (1e-6)^2*eye(78);

    [Xhat,Z,Pm,Pp] = observer.estimateState(Y,U,Qvar,Rvar);
    disp(['Observer Elapsed Time: ',num2str(toc(observerElapsed))])
    
    %norm(Xhat.p-Y.p)
    %     structurePlot(Xhat,omega,constraints,cameraViewInitial,1,0,0)
    pause(1e-4)
    
    %% Controller (takes Xhat outputs U)
    %examples: constrained QP MPC, iLQR, LQR
    %input arguments - Xhat,Uhat,pDDOT,jacobians,hFcns
    controllerElapsed = tic;   
    
    [U,OL_states,OL_inputs,hVars,cost,controllerOutputs] =...
        controller.getOptimalInput(...
        Xhat,U,nominalFcn,jacobianFcns,hFcns,costFcnHandle,...
        debugFcns,openLoopFlag);
    
    disp(['Controller Elapsed Time: ',num2str(toc(controllerElapsed))])
    
    %% inner-loop checks
    xCOM = mean(plant.Current_p(1:3:end));
    yCOM = mean(plant.Current_p(2:3:end));
    zCOM = mean(plant.Current_p(3:3:end));
    centroid = [xCOM;yCOM;zCOM];
    if(norm(centroid(1:2)-controllerTargets{targetIdx})<0.5)
        targetIdx = targetIdx+1;
        %resetCables = 1; %comment out if using MPC
        if(targetIdx <= numel(controllerTargets))
            controller.setTargetDestination(controllerTargets{targetIdx});
        else
            return
        end
    end
    if(resetCables == 1 && norm(omega.U.RL0-plant.Current_RL)<=1e-3)
        resetCables = 0;
    end
    
    %% Save Data
    if(doLog)
        X_record.p(:,iteration) = plant.Current_p;
        X_record.pDOT(:,iteration) = plant.Current_pDOT;
        X_record.RL(:,iteration) = plant.Current_RL;
        X_record.L(:,iteration) = plant.Current_L;
        Xhat_record.p(:,iteration) = Xhat.p;
        Xhat_record.pDOT(:,iteration) = Xhat.pDOT;
        Xhat_record.RL(:,iteration) = Xhat.RL;
        Xhat_record.L(:,iteration) = Xhat.L;
        Pm_record(:,iteration) = Pm;
        Pp_record(:,iteration) = Pp;
        U_record.RLdot(:,iteration) = U.RLdot;
        U_record.Ldot(:,iteration) = U.Ldot;
        X_openLoop_record{iteration} = OL_states; %save memory
        U_openLoop_record{iteration} = OL_inputs; %save memory
        Xbar.p = plant.Current_p;
        Xbar.pDOT = plant.Current_pDOT;
        Xbar.RL = plant.Current_RL;
        Xbar.L = plant.Current_L;
        Ubar = U;
        Gamma_record(:,iteration) = Gamma(Xbar,Ubar,hVars);
        GenForce_record(:,iteration) = genForces(Xbar,Ubar,hVars);
        ConstrForce_record(:,iteration) = constrForces(Xbar,Ubar,hVars);
        cost_record(iteration) = cost;
        controllerOutputArgs_record{iteration} = controllerOutputs;
        if(mod(iteration,savePeriod)==0)
            workspaceVars = who;
            try
                save(['./Results/LatestResults_MPC_',num2str(timeStamp),'_',...
                    modelName,'_iLQR'],...
                    workspaceVars{not(~contains(workspaceVars, 'record'))})
            catch
                disp('File Temporarily Unavailable')
            end
        end
    end    
end

