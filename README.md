# Tensegrity-Dynamics-and-Optimal-Control
Tensegrity Robots are a biologically-inspired approach to building robots based on the tension networks of tensegrity structures, which have no rigid connections between elements. This repository contains some of my research and relevant software tools related to optimal control of tensegrity robots!

Note - This is still an active work in progress. There still remains a lot of work to be done to ensure bug-free code and optimized execution. If you'd like to contribute to this Tensegrity research, please fork this repo and feel free to send useful pull requests. I'll be sure to add you as a contributor down below :)

![rolling tensegrity](/Images/rolling_tensegrity.png)

This repository is a simulation framework that enables rapid prototyping of tensegrity designs, actuation and sensing topologies, control policies, and state observers. To this end, the code in this repository follows object-oriented design principles, where the System Plant, Controllers, and Observers are all independently encapsulated to enable quick plug-and-play experimentation of any system component. As an example, you may choose to implement Model Predictive Control with full-state information or experiment with iLQR and an Unscented Kalman Filter with noisy sensor data. For a published work example of this simulation framework, check out this [IROS 2018 paper](https://ieeexplore.ieee.org/document/8594401) which uses some of the tools in this repository to design a contextual neural network policy using optimal state-action trajectories generated using Model Predictive Control.

![top-level design](/Images/System_Design.png)

## Tensegrity Dynamics

This work also features automatic generation of dynamics equations for any Type-1 tensegrity system (connection nodes are connected to at most *one* rigid body), based on a simplified pointmass dynamics formulation. While this trades off simulation accuracy, this approach enables rapid prototyping of different mechanical designs and sensor/actuator configurations to understand broad patterns of system feasibility, depending on design requirements. 

![dynamics](/Images/dynamics.png)

## Getting Started

This repository depends on some recent MATLAB packages. From our testing, the codebase should work for MATLAB 2018a onwards. After cloning the repo, everything MATLAB-related should work out of the box, once the dependencies (e.g., YALMIP, Guribi, MATLAB Optimization Toolbox Add-On) are properly installed. 

Note: Remember to add all folders and sub-folders in this repository on MATLAB's search path directory ('Set Path' in the 'Home' tab).



**TODO: Include information regarding YALMIP and GUROBI**

## Running Examples

In order to run the framework, call the top-level function:

```bash
>> TensegritySimulation
```

To examine and visualize a specific tensegrity topology, run the associate model file then run structurePlot.m:
```bash
>> [omega,X] = six_bar_model
>> structurePlot(X,omega,constraints,3,1,1,1,1,1)
```

To simply run a simulated drop test, call the top-level function which features no actuation or state estimation:
```bash
>> SimpleDynamicSimulation
```

**TODO: Explain User-Defined Parameters (Maybe consider a separate config file?)**


## Authors
* **Brian Cera** - *Initial Work*

