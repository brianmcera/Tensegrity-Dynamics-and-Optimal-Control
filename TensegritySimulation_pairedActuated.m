function TensegritySimulation_pairedActuated()
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
modelName = 'I_six_bar_model_pairedActuated';
%plant filename and custom input arguments
plantName = 'DefaultPlant';
%observer filename and custom input arguments
observerName = 'DefaultObserver';
%controller filename and custom input arguments
controllerName = 'QP_MPC_RollingDirection2';
%cost filename and custom input arguments
costName = 'velocityCost'; %specify costFunction parameters below

simTimeStep = 1e-2; %simulation timestep
totalSimSteps = 750; %total number of simulation timesteps
controllerHorizon = 10; %MPC horizon for controller
actuationMode = 1; %1-cables, 2-rods, 3-both 
%todo: rod actuation does not work right now

%save data
doLog = true;
savePeriod = 50; %how often to save data (e.g., every 50 timesteps)

%forward simulation / initial conditions
show_initialization_graphics = false; %true/false
perturbCablesPercentage = 0.0; %between 0 and 1
perturbRodsPercentage = 0.00; %between 0 and 1
XYrandomRotate = false; %true/false
YZrandomRotate = false; %true/false
KineticEnergyDamping = 1; %true/false
cameraViewInitial = [90 0]; %[az,el] or {1,2,3}

%simulation loop
show_simulation_graphics = false;
cameraViewLoop = 3; %syntax: either A)[az,el] or B)integer,{1,2,3}



%//////////////////////////////////////////////////////////////////////////
%% Generate Model and Dynamics
[omega,X,U] = feval(modelName);
constraints = omega.constraints;
generalForces = omega.generalForces;
structurePlot(X,omega,constraints,cameraViewInitial,1);
pause(1e-3)
[pDDOTFcn,hFcns,jacobianFcns,Gamma,genForces,constrForces,debugFcns] = ...
    Dynamics_Generator(omega,constraints,generalForces);

%% Set up desired initial conditions
%calculate initial cable/rod lengths for forward simulation
%helper function handles
z = @(p,V) -kron(V'*V,eye(3))*p; %tension direction in p-basis
separationDist = @(p,V) sqrt(z(p,V)'*z(p,V)/2);

%initialize cable restlengths
C = omega.C;
for cable = 1:size(C,1)
    U_desired.RL(cable) = separationDist(X.p,C(cable,:))-...
        omega.cables.pretension(cable)/omega.cables.stiffness(cable);
end
U_desired.RL = U_desired.RL'; %column vector

%initialize rod lengths
R = omega.R;
for rod = 1:size(R,1)
    U_desired.L(rod) = separationDist(X.p,R(rod,:));
end
U_desired.L = U_desired.L'; %column vector

%store neutral inputs of pretensioned robot
omega.U.RL0 = U_desired.RL;
omega.U.L0 = U_desired.L;

%perturb cable lengths around nominal pretensioned lengths
if(perturbCablesPercentage~=0)
    disp('~~"perturbCablesPercentage" parameter not equal to nominal "0": Perturbing nominal pretensioned cable restlengths~~')
    U_desired.RL = (eye(numel(U_desired.RL))+...
        diag(perturbCablesPercentage*rand(numel(U_desired.RL),1)-...
        perturbCablesPercentage/2))*U_desired.RL;
    %constrained cables
    if(~isempty(omega.cables.passive))
        U_desired.RL(omega.cables.passive) = omega.U.RL0(...
            omega.cables.passive);
    end
    %paired cables
    if(~isempty(omega.cables.paired))
        U_desired.RL(omega.cables.paired(:,2)) = ...
            -U_desired.RL(omega.cables.paired(:,1))+...
            omega.U.RL0(omega.cables.paired(:,1))+...
            omega.U.RL0(omega.cables.paired(:,2));
    end
end
if(perturbRodsPercentage~=0)
    U_desired.L = (eye(numel(U_desired.L))+...
        diag(perturbCablesPercentage*rand(numel(U_desired.L),1)-...
        perturbCablesPercentage/2))*U_desired.L;
    %constrained rods
    if(~isempty(omega.rods.constrained))
        U_desired.L(omega.rods.constrained) = omega.U.L0(...
            omega.rods.constrained);
    end
end

%rotate robot about vertical Z-axis
if(XYrandomRotate)
    disp('~~"XYrandomRotate" parameter set to TRUE: Randomly rotating nodal positions about vertical Z-axis~~')
    %center XY
    centeredX = X.p(1:3:end)-mean(X.p(1:3:end));
    centeredY = X.p(2:3:end)-mean(X.p(2:3:end));
    randRotAngle = rand*(2*pi);
    rotatedXY = [cos(randRotAngle),-sin(randRotAngle);...
        sin(randRotAngle),cos(randRotAngle)]*[centeredX';centeredY'];
    X.p(1:3:end) = rotatedXY(1,:); %rotated X values
    X.p(2:3:end) = rotatedXY(2,:); %rotated Y values
    structurePlot(X,omega,constraints,cameraViewInitial,1);
end

%rotate robot about X-axis
if(YZrandomRotate)
    disp('~~"YZrandomRotate" parameter set to TRUE: Randomly rotating nodal positions about horizontal X-axis~~')
    %center YZ
    centeredY = X.p(2:3:end)-mean(X.p(2:3:end));
    centeredZ = X.p(3:3:end)-mean(X.p(3:3:end));
    randRotAngle = rand*(2*pi);
    rotatedYZ = [cos(randRotAngle),-sin(randRotAngle);...
        sin(randRotAngle),cos(randRotAngle)]*[centeredY';centeredZ'];
    X.p(2:3:end) = rotatedYZ(1,:); %rotated X values
    X.p(3:3:end) = rotatedYZ(2,:); %rotated Y values
    structurePlot(X,omega,constraints,cameraViewInitial,1);
end
 
pause(1e-3) %give time to update structurePlot

%% Instantiate Cost function handle & Plant,Controller,Observer objects
%Handle and Objects created here after model instantiation
%'omega' parameter here reflects updated cable/rod lengths after any
%modifications to nominal model file (e.g., cable perturbations, model
%rotations, etc.

%cost function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%user-defined arguments to cost function
costArgs.RL_diff_weight = 3;
costArgs.RL_actuation_weight = 3;%0
costArgs.L_diff_weight = 5;
costArgs.velocity_reward = 1;
costArgs.omega = omega; %latest omega after calculating RL0/L0 above
costArgs.stepDiscount = 0.95;
costFcnHandle = @(X,U,runtimeArgs)feval(costName,X,U,costArgs,runtimeArgs);

%plant object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plant = feval(plantName,X,U,omega,simTimeStep);

%observer object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
observer = feval(observerName);

%controller object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% controllerTargets = {[-3*sqrt(2),-3]',[3*sqrt(2),-3]',[0,0]'};
controllerTargets = {10*[cos(rand(1)*2*pi),sin(rand(1)*2*pi)]'};%{[3,0]',[3,-3]',[0,-3]',[0,0]'};
targetIdx = 1;
controller = feval(controllerName,X,U,omega,simTimeStep,controllerHorizon);
controller.setActuationMode(actuationMode);
controller.setTargetDestination(controllerTargets{targetIdx});
%controller.createOptimizerObject(costFcnHandle); %only for YALMIP MPC, optimizer object
% controller.setLSTMnet('Latest_FFNN_19-01-25_22_12_32'); %only for LSTM neural network

%% Forward Simulate to Equilibrium (with Kinetic Energy Damping)
forwardSimIter = 1;
kineticEnergy_prev = 0;
cablesStillMoving = 1; %debug for now
while((plant.getKineticEnergy()>1e-4)|| plant.getKineticEnergy()==0 ||...
        cablesStillMoving)
    
    %advance simulation one step
    plant.stepForward(U_desired,pDDOTFcn,hFcns);
    X.p = plant.Current_p;
    X.pDOT = plant.Current_pDOT;
    
    %plot
    if(show_initialization_graphics)
        structurePlot(X,omega,constraints,cameraViewInitial,1);
        title(['t = ',num2str(plant.simulationTime)])
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
            disp(['Kinetic Damping, Current KE: ',...
                num2str(kineticEnergy_prev)])
        end
    end
    
    %cables reach desired initial (perturbed) state
    if(all(plant.Current_RL == U_desired.RL))
        cablesStillMoving = 0;
    end
    
    forwardSimIter = forwardSimIter+1;
end
plant.resetClock();
Y = plant.outputSensors();


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

X_record.p = zeros(size(X.p,1),totalSimSteps);
X_record.pDOT = zeros(size(X.pDOT,1),totalSimSteps);
U_record.RL = zeros(size(U.RL,1),totalSimSteps);
U_record.L = zeros(size(U.L,1),totalSimSteps);
Xhat_record.p = zeros(size(X.p,1),totalSimSteps);
Xhat_record.pDOT = zeros(size(X.pDOT,1),totalSimSteps);
Uhat_record.RL = zeros(size(U.RL,1),totalSimSteps);
Uhat_record.L = zeros(size(U.L,1),totalSimSteps);
X_openLoop_record = cell(1,totalSimSteps);
U_openLoop_record = cell(1,totalSimSteps);
Gamma_record = zeros(size(X.p,1),totalSimSteps);
GenForce_record = zeros(size(X.p,1),totalSimSteps);
ConstrForce_record = zeros(size(X.p,1),totalSimSteps);
cost_record = zeros(1,totalSimSteps);
controllerOutputArgs_record = cell(1,totalSimSteps);


%% Simulation Loop

timeStamp = datestr(now,'yy-mm-dd_HH_MM_SS');
LoopStart = tic;
resetCables = 0;
for iteration = 1:totalSimSteps
    disp('')
    disp(['Simulation Loop Timestamp: ',timeStamp])
    disp(['Iteration: ',num2str(iteration)])
    disp(['Total Elapsed Time: ',num2str(toc(LoopStart))])
    
    %% Observer (takes Y, outputs Xhat,Uhat)
    %examples:full-state information, Luenberg Observer, Kalman Filter,
    %Extended Kalman Filter, Unscented Kalman Filter
    observerElapsed = tic;
    [Xhat,Uhat] = observer.estimateState(Y);
    disp(['Observer Elapsed Time: ',num2str(toc(observerElapsed))])
    
    %% Controller (takes Xhat,Uhat outputs U)
    %examples: constrained QP MPC, iLQR, LQR
    %input arguments - Xhat,Uhat,pDDOT,jacobians,hFcns
    controllerElapsed = tic;        
    [U_desired,OL_states,OL_inputs,hVars,cost,controllerOutputs] =...
        controller.getOptimalInput(...
        Xhat,Uhat,pDDOTFcn,jacobianFcns,hFcns,costFcnHandle,debugFcns);
    
    %reset cable state <-- incorporate this into the controller#################################
    if(resetCables==1)
        controllerOutputs.resetCables = 1;
        U_desired.RL = omega.U.RL0;
    else
        controllerOutputs.resetCables = 0;
    end
    disp(['Controller Elapsed Time: ',num2str(toc(controllerElapsed))])
    
    %% Simulate Plant (takes U, outputs Y)
    %handle cable actuation, motor dynamics, input saturation, input/output
    %noise
    plantElapsed = tic;
    plant.stepForward(U_desired,pDDOTFcn,hFcns);
    Y = plant.outputSensors();    
    %plot results
    if(show_simulation_graphics)
        X.p = plant.Current_p;
        X.pDOT = plant.Current_pDOT;
        structurePlot(X,omega,constraints,cameraViewLoop,1);
        title(['t = ',num2str(plant.simulationTime)])
    end
    disp(['Simulate Plant Elapsed Time: ',num2str(toc(plantElapsed))])
    
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
        U_record.RL(:,iteration) = plant.Current_RL;
        U_record.L(:,iteration) = plant.Current_L;
        Xhat_record.p(:,iteration) = Xhat.p;
        Xhat_record.pDOT(:,iteration) = Xhat.pDOT;
        Uhat_record.RL(:,iteration) = Uhat.RL;
        Uhat_record.L(:,iteration) = Uhat.L;
        X_openLoop_record{iteration} = OL_states; %save memory
        U_openLoop_record{iteration} = OL_inputs;
        Xbar.p = plant.Current_p;
        Xbar.pDOT = plant.Current_pDOT;
        Ubar.RL = plant.Current_RL;
        Ubar.L = plant.Current_L;
        Gamma_record(:,iteration) = Gamma(Xbar,Ubar,hVars);
        GenForce_record(:,iteration) = genForces(Xbar,Ubar,hVars);
        ConstrForce_record(:,iteration) = constrForces(Xbar,Ubar,hVars);
        cost_record(iteration) = cost;
        controllerOutputArgs_record{iteration} = controllerOutputs;
        if(mod(iteration,savePeriod)==0)
            workspaceVars = who;
            try
                save(['./Results/LatestResults_MPC_fastMotors_OL',num2str(timeStamp),'_',...
                    modelName,'_pairedActuated'],...
                    workspaceVars{not(~contains(workspaceVars, 'record'))})
            catch
                disp('File Temporarily Unavailable')
            end
        end
    end    
end

