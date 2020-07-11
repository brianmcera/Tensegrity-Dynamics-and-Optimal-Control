%  TENSEGRITYSIMULATION:  A top-level function to simulate general
%  tensegrity dynamics, subject to custom controllers, observers, and plant
%  dynamics. Custom configuration parameters can be modified in first code
%  section titled USER-DEFINED PARAMETERS.
%
%  Usage:      >> TensegritySimulation
%
%  Arguments:  N/A
%
%  Returns:    N/A
%
%  Author:     Brian Cera
%  Date:       Aug 2019
%  Modified:   Per Brian Cera 10/15/2019

close all
clear all
clc
rng('shuffle') 
 
%% User-defined Parameters ////////////////////////////////////////////////

% model filename
model_name = 'six_bar_model';
% plant filename
plant_name = 'DefaultPlant';
% observer filename
observer_name = 'DefaultObserver';
% controller filename
controller_name = 'iLQRminimax_RollingDirection_inputpenalty';
% cost filename
cost_name = 'velocityCost'; % NOTE:specify costFunction parameters below

sim_time_step = 5e-3;       % simulation timestep
total_sim_steps = 5000;     % total number of simulation timesteps
controller_horizon = 10;%5; % MPC horizon for controller
actuation_mode = 1;         % 1-cables, 2-rods, 3-both 

% dynamics generation parameters
save_dynamics_file = false;
optimize_flag = false;
optimized_dynamics_filename = 'optimizedDynamics_SixBar_20_03_26_10_19_07';

% save data
log_toggle =  true; % toggle flag to save simulation data to external file
save_period = 50;   % how often to save data (e.g., every 50 timesteps)

% forward simulation / initial conditions
show_initialization_graphics = false;   % true/false
kinetic_energy_damping = true;          % true/false
perturb_cables_percent = 0.00;          % double between 0 and 1
perturb_rods_percent = 0.00;            % double between 0 and 1
xy_random_rotate = false;               % true/false
yz_random_rotate = false;               % true/false
camera_view_initial = [90, 0];          % format: [az,el] or {1,2,3}

%simulation loop
show_simulation_graphics = false;   % toggle visualization in loop
open_loop_flag = true;              % toggle open-loop calculation in loop 
camera_view_loop = [0,15];          % either A)[az,el] or B)integer,{1,2,3}


%//////////////////////////////////////////////////////////////////////////

%% Generate Model and Dynamics
time_stamp = datestr(now,'yy_mm_dd_HH_MM_SS');
[omega,X] = feval(model_name);
constraints = omega.constraints;
general_forces = omega.generalForces;
structurePlot(X,omega,constraints,camera_view_initial,1,0,1,1,1);
pause(1e-3)
[nominal_fcn,hFcns,jacobian_fcns,Gamma,gen_forces,constr_forces,debug_fcns] = ...
    Dynamics_Generator(omega,constraints,general_forces);

% build symbolic matlabFunction which optimizes computation time
% longer initialization, faster iteration loops
if(isfile(['./Dynamics/Generated Dynamics/',optimized_dynamics_filename,'.m']))
    disp(['Existing dynamics function found: <'...
        './Dynamics/Generated Dynamics/',optimized_dynamics_filename,'.m>'])
    disp('Using pre-generated dynamics file...')
    nominal_fcn.pDDOT = @(x1,x2,x3,x4)feval(optimized_dynamics_filename,...
        x1,x2,x3,x4);
else
    disp('Existing dynamics file NOT found, generating dynamics functions...')
    % symbolic State structure
    Xsym = struct();
    fields = fieldnames(X);
    
    % parse through all struct fields
    for fieldIdx = 1:numel(fields)
        Xsym.(fields{fieldIdx}) = ...
            sym(fields{fieldIdx},size(X.(fields{fieldIdx})),'real');
    end
    
    % symbolic Input structure
    Usym = struct();
    Usym.RLdot = sym('RLdot',size(omega.cableConstraintMatrix,2),'real');
    Usym.Ldot = sym('Ldot',size(omega.rodConstraintMatrix,2),'real');

    % symbolic Helper Variables structure
    hVarSym = struct();
    fields = fieldnames(hFcns);
    for fieldIdx = 1:numel(fields)-1
        if(~isa(hFcns.(fields{fieldIdx}){fieldIdx},'function_handle'))
            % doesn't take Xsym as an argument (element is not function handle)
            for i = 1:numel(hFcns.(fields{fieldIdx}))
                hVarSym.(fields{fieldIdx}){i} = ...
                    hFcns.(fields{fieldIdx}){i};
            end
        else
            % takes Xsym as an argument (element is function handle)
            for i = 1:numel(hFcns.(fields{fieldIdx}))
                hVarSym.(fields{fieldIdx}){i} = ...
                    hFcns.(fields{fieldIdx}){i}(Xsym);
            end
        end
    end
    hVarSym.J = hFcns.J(Xsym,Usym,hVarSym);
    disp('Generating Dynamic Function Handles...')
    if(optimize_flag)
        disp(['Optimizing Dynamics This may take >15 minutes ',...
            '(optimizeFlag variable is set to "true")...'])
    end
    if(save_dynamics_file)
        if(optimize_flag)
            pDDOT_simplified = nominal_fcn.pDDOT(Xsym,Usym,hVarSym);
            for i = 1:numel(pDDOT_simplified)
                disp(['Simplifying element ',num2str(i),'/',num2str(numel(pDDOT_simplified))])
                simp_time = tic;
                pDDOT_simplified(i) = simplify(pDDOT_simplified(i));
                disp(['Elapsed Time: ',num2str(toc(simp_time))])
            end
        end
        nominal_fcn.pDDOT = matlabFunction(pDDOT_simplified,...
            'File',['Dynamics/Generated Dynamics/',optimized_dynamics_filename,...
            '_',num2str(time_stamp)],...
            'Optimize',optimize_flag,'Vars',{Xsym.p,Xsym.pDOT,Xsym.RL,Xsym.L});
    else
        nominal_fcn.pDDOT = matlabFunction(nominal_fcn.pDDOT(Xsym,Usym,hVarSym),...
            'Optimize',optimize_flag,'Vars',{Xsym.p,Xsym.pDOT,Xsym.RL,Xsym.L});
    end
end

%store dynamics into 'omega' struct
omega.nominalFcn = nominal_fcn;
omega.hFcns = hFcns;
omega.jacobianFcns = jacobian_fcns;


%% Set up desired initial conditions
% calculate initial cable/rod lengths for forward simulation

% helper function handles
z = @(p,V) -kron(V'*V,eye(3))*p; % tension direction in p-basis
separationDist = @(p,V) sqrt(z(p,V)'*z(p,V)/2); % euclidean separation dist.

% initialize cable restlengths
C = omega.C;
for cable = 1:size(C,1) % apply pretension
    X_desired.RL(cable) = separationDist(X.p,C(cable,:))-...
        omega.cables.pretension(cable)/omega.cables.stiffness(cable);
end
X_desired.RL = X_desired.RL'; % column vector

% initialize rod lengths
R = omega.R;
for rod = 1:size(R,1)
    X_desired.L(rod) = separationDist(X.p,R(rod,:));
end
X_desired.L = X_desired.L'; % column vector

% store neutral inputs of equally pretensioned robot
omega.X.RL0 = X_desired.RL;
omega.X.L0 = X_desired.L;

% perturb cable lengths around nominal pretensioned lengths
if(perturb_cables_percent~=0)
    disp(['~~"perturbCablesPercentage" parameter not equal to nominal value of "0":',...
        ' Perturbing nominal pretensioned cable restlengths~~'])
    X_desired.RL = (eye(numel(X_desired.RL))+...
        diag(perturb_cables_percent*rand(numel(X_desired.RL),1)-...
        perturb_cables_percent/2))*X_desired.RL;
    % handle constrained cables
    if(~isempty(omega.cables.passive))
        X_desired.RL(omega.cables.passive) = omega.X.RL0(...
            omega.cables.passive);
    end
    % handle paired cables
    if(~isempty(omega.cables.paired))
        X_desired.RL(omega.cables.paired(:,2)) = ...
            -X_desired.RL(omega.cables.paired(:,1))+...
            omega.X.RL0(omega.cables.paired(:,1))+...
            omega.X.RL0(omega.cables.paired(:,2));
    end
    % handle similarly actuated cables
    if(~isempty(omega.cables.similar))
        X_desired.RL(omega.cables.similar(:,2)) = ...
            X_desired.RL(omega.cables.similar(:,1));
    end
end
if(perturb_rods_percent~=0)
    X_desired.L = (eye(numel(X_desired.L))+...
        diag(perturb_cables_percent*rand(numel(X_desired.L),1)-...
        perturb_cables_percent/2))*X_desired.L;
    % constrained rods
    if(~isempty(omega.rods.constrained))
        X_desired.L(omega.rods.constrained) = omega.X.L0(...
            omega.rods.constrained);
    end
end

% rotate robot about vertical Z-axis
if(xy_random_rotate==true)
    disp(['~~"XYrandomRotate" parameter set to TRUE: ',...
        'Randomly rotating nodal positions about vertical Z-axis~~'])
    % center XY
    centeredX = X.p(1:3:end)-mean(X.p(1:3:end));
    centeredY = X.p(2:3:end)-mean(X.p(2:3:end));
    randRotAngle = rand*(2*pi);
    rotatedXY = [cos(randRotAngle),-sin(randRotAngle);...
        sin(randRotAngle),cos(randRotAngle)]*[centeredX';centeredY'];
    X.p(1:3:end) = rotatedXY(1,:); % rotated X values
    X.p(2:3:end) = rotatedXY(2,:); % rotated Y values
    structurePlot(X,omega,constraints,camera_view_initial,1,0,1,1,1);
end

% rotate robot about X-axis
if(yz_random_rotate==true)
    disp(['~~"YZrandomRotate" parameter set to TRUE: Randomly ',
        'rotating nodal positions about horizontal X-axis~~'])
    % center YZ
    minZ = min(X.p(3:3:end));
    centeredY = X.p(2:3:end)-mean(X.p(2:3:end));
    centeredZ = X.p(3:3:end)-mean(X.p(3:3:end));
    randRotAngle = rand*(2*pi);
    rotatedYZ = [cos(randRotAngle),-sin(randRotAngle);...
        sin(randRotAngle),cos(randRotAngle)]*[centeredY';centeredZ'];
    X.p(2:3:end) = rotatedYZ(1,:); % rotated X values
    X.p(3:3:end) = rotatedYZ(2,:); % rotated Y values
    X.p(3:3:end) = X.p(3:3:end)-min(X.p(3:3:end)) + minZ;
    structurePlot(X,omega,constraints,camera_view_initial,1,0,0,1,1);
end

pause(1e-3) % give time to update and redraw structurePlot


% Instantiate Cost function handle & Plant,Controller,Observer objects:
% Handle and Objects created here after model instantiation
% Note: 'omega' parameter here now reflects updated cable/rod lengths after 
% any modifications to nominal model file (e.g., cable perturbations, model
% rotations, etc.

%% cost function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function encapsulates all of the necessary parameters which can be 
% tuned to produce a user-define desired behavior from the controller.
% Designate the desired cost function by filename in the 'User-Defined
% Parameters' section above.

% user-defined arguments to cost function
cost_args.RL_diff_weight = 1e-2;
cost_args.RL_actuation_weight = 0;
cost_args.L_diff_weight = 0;
cost_args.velocity_reward = 50;
cost_args.omega = omega; 
cost_args.stepDiscount = 0.95;
cost_args.N = controller_horizon;
costFcnHandle = @(X,U,runtimeArgs)feval(cost_name,X,U,cost_args,runtimeArgs);
disp('Generating Cost Function Handle...')

%% plant object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This object is the simulated object which represents the actual robot.
% Specifically, this object keeps track of the robot state throughout the
% simulation. This is also where any actuator input noise and sensor output
% noise should be incorporated. Define the desired Plant in the
% 'User-Defined Parameters' section above.

plant = feval(plant_name,X,omega,sim_time_step);

%% Forward Simulate to Equilibrium (with Kinetic Energy Damping)
% note: this occurs before Controller instantiation to update 'omega' struct
forwardSimIter = 1;
kineticEnergy_prev = 0;
cablesStillMoving = 1;

if(any(X_desired.RL<omega.cables.minLength))
    error(['Desired restlength is less than minimum cable ',...
        'length limits. Increase stiffness or decrease desired pretension'])
end

disp('Starting forward simulation to equilibrium...')
if(kinetic_energy_damping)
    disp('Kinetic Damping, Current KE:')
end
while((plant.getKineticEnergy()>1e-4)||... % non-negligible kinetic energy
        plant.getKineticEnergy()==0 ||... % KE just got reset
        cablesStillMoving) % cables haven't yet reached goal rest lengths
    
    % calculate appropriate cable length / rod length
    U.RLdot = sign(X_desired.RL - plant.Current_RL).*...
        min(omega.cables.linear_velocity,...
        abs(X_desired.RL - plant.Current_RL)/sim_time_step);
    U.Ldot = zeros(size(plant.Current_L));
    
    % advance simulation one timestep
    plant.stepForward(U,nominal_fcn,hFcns,true);
    X.p = plant.Current_p;
    X.pDOT = plant.Current_pDOT;
    
    % plot update
    if(show_initialization_graphics==true)
        structurePlot(X,omega,constraints,camera_view_initial,1,0,1,1,1);
        title(['t = ',num2str(plant.simulationTime)])
        drawnow()
    end
    
    % kinetic energy damping
    if(kinetic_energy_damping)
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
    
    % cables reach desired initial (perturbed) state
    if(norm(plant.Current_RL-X_desired.RL)<1e-6)
        cablesStillMoving = 0;
    end
    
    forwardSimIter = forwardSimIter+1;
end
plant.resetClock(); % reset simulation time for recordkeeping
omega.X.p0 = plant.Current_p; % record equilibrium state for reference


%% controller object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This object handles optimal control of the robot. A variety of approaches
% and algorithms can be plug-and-play inserted here. Examples include linear
% time-varying MPC, iLQR, neural network, etc. Define the desired
% Controller by filename in the 'User-Defined Parameters' section above.

controllerTargets = {[30,0]'};
targetIdx = 1;
controller = feval(controller_name,X,omega,sim_time_step,controller_horizon);
controller.setActuationMode(actuation_mode);
controller.setTargetDestination(controllerTargets{targetIdx});

%% observer object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This object is the observer/estimator, which takes in the full information
% from the Plant object, simulates partially observable sensor data, and 
% outputs the best estimate of the current state.
% Additionally, any information about uncertainty/variance should be handled
% by this object. Define the desired Observer by filename in the 'User-
% Defined Parameters' section above.

observer = feval(observer_name,nominal_fcn,hFcns,jacobian_fcns,...
    X,omega,sim_time_step);
Ydim = 84; % dimension of sensor observation vector

%% Instantiate Record Arrays for Data Storage
% **include 'record' in filename for automated storage**
simulationParameters_record.omega = omega;
simulationParameters_record.model = model_name;
simulationParameters_record.plant = plant_name;
simulationParameters_record.controller = controller_name;
simulationParameters_record.observer = observer_name;
simulationParameters_record.costName = cost_name;
simulationParameters_record.costArgs = cost_args;
simulationParameters_record.timestep = sim_time_step;
simulationParameters_record.controllerHorizon = controller_horizon;

% dimension of concatenated state vector
Xdim = numel(X.p)+numel(X.pDOT)+numel(X.RL)+numel(X.L);

X_record.p = zeros(size(X.p,1),total_sim_steps);
X_record.pDOT = zeros(size(X.pDOT,1),total_sim_steps);
X_record.RL = zeros(size(X.RL,1),total_sim_steps);
X_record.L = zeros(size(X.L,1),total_sim_steps);
U_record.RLdot = zeros(size(X.RL,1),total_sim_steps);
U_record.Ldot = zeros(size(X.L,1),total_sim_steps);
Xhat_record.p = zeros(size(X.p,1),total_sim_steps);
Xhat_record.pDOT = zeros(size(X.pDOT,1),total_sim_steps);
Xhat_record.RL = zeros(size(X.RL,1),total_sim_steps);
Xhat_record.L = zeros(size(X.L,1),total_sim_steps);
Z_record = zeros(Ydim,total_sim_steps);
Pm_record = zeros(Xdim,total_sim_steps);
Pp_record = zeros(Xdim,total_sim_steps);
X_openLoop_record = cell(1,total_sim_steps);
U_openLoop_record = cell(1,total_sim_steps);
Gamma_record = zeros(size(X.p,1),total_sim_steps);
GenForce_record = zeros(size(X.p,1),total_sim_steps);
ConstrForce_record = zeros(size(X.p,1),total_sim_steps);
cost_record = zeros(1,total_sim_steps);
controllerOutputArgs_record = cell(1,total_sim_steps);


%% Simulation Loop

LoopStart = tic;
resetCables = 0;

% initial zero input for timestep t=0;
U.RLdot = zeros(size(X.RL,1),1);
U.Ldot = zeros(size(X.L,1),1);

for iteration = 1:total_sim_steps
    disp('')
    disp(['Simulation Loop Timestamp: ',time_stamp])
    disp(['Iteration: ',num2str(iteration)])
    disp(['Total Elapsed Time: ',num2str(toc(LoopStart))])
    
    %% PLANT (takes U, outputs Y)
    % handle cable actuation, motor dynamics, input saturation, input/output
    % noise
    
    plantElapsed = tic;
    plant.stepForward(U,nominal_fcn,hFcns,false);
    Y = plant.outputSensors();    
    % plot results
    if(show_simulation_graphics)
        X.p = plant.Current_p;
        X.pDOT = plant.Current_pDOT;
        structurePlot(X,omega,constraints,camera_view_loop,1,0,1,1,1);
        title(['t = ',num2str(plant.simulationTime)])
    end
    disp(['Simulate Plant Elapsed Time: ',num2str(toc(plantElapsed))])
    
    %% OBSERVER (takes Y, outputs Xhat)
    % examples:full-state information, Luenberg Observer, Kalman Filter,
    % Extended Kalman Filter, Unscented Kalman Filter
    
    observerElapsed = tic;
    Qvar = (1e-2)^2*eye(size(X.p,1)+size(X.pDOT,1)+size(X.RL,1)+size(X.L,1));
    Rvar = diag([(1e-2)^2*ones(18,1);
        (10e-2)^2*ones(24,1);
        (1e-2)^2*ones(36,1);
        (1e-6)^2*ones(6,1)]);% observer vector dimension = 84
    [Xhat,Z,Pm,Pp] = observer.estimateState(Y,U,Qvar,Rvar);
    disp(['Observer Elapsed Time: ',num2str(toc(observerElapsed))])
    pause(1e-4)
    
    %% CONTROLLER (takes Xhat outputs U)
    % examples: constrained QP MPC, iLQR, LQR
    % input arguments - Xhat,Uhat,pDDOT,jacobians,hFcns
    
    controllerElapsed = tic;     
    [U,OL_states,OL_inputs,hVars,cost,controllerOutputs] =...
        controller.getOptimalInput(...
        Xhat,U,nominal_fcn,jacobian_fcns,hFcns,costFcnHandle,...
        debug_fcns,open_loop_flag); 
    disp(['Controller Elapsed Time: ',num2str(toc(controllerElapsed))])
    
    %% inner-loop checks
    xCOM = mean(plant.Current_p(1:3:end));
    yCOM = mean(plant.Current_p(2:3:end));
    zCOM = mean(plant.Current_p(3:3:end));
    centroid = [xCOM;yCOM;zCOM];
    if(norm(centroid(1:2)-controllerTargets{targetIdx})<0.5)
        targetIdx = targetIdx+1;
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
    if(log_toggle)
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
        X_openLoop_record{iteration} = OL_states; % save memory
        U_openLoop_record{iteration} = OL_inputs; % save memory
        Xbar.p = plant.Current_p;
        Xbar.pDOT = plant.Current_pDOT;
        Xbar.RL = plant.Current_RL;
        Xbar.L = plant.Current_L;
        Ubar = U;
        Gamma_record(:,iteration) = Gamma(Xbar,Ubar,hVars);
        GenForce_record(:,iteration) = gen_forces(Xbar,Ubar,hVars);
        ConstrForce_record(:,iteration) = constr_forces(Xbar,Ubar,hVars);
        cost_record(iteration) = cost;
        controllerOutputArgs_record{iteration} = controllerOutputs;
        if(mod(iteration,save_period)==0)
            workspaceVars = who;
            try
                save(['./Results/LatestResults_MPC_',num2str(time_stamp),'_',...
                    model_name,'_',controller_name,'_G1e-3'],...
                    workspaceVars{not(~contains(workspaceVars, 'record'))})
            catch
                disp('File Temporarily Unavailable')
            end
        end
    end    
end

