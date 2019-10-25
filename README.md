# Tensegrity-Dynamics-and-Optimal-Control
Tensegrity Robots are a biologically-inspired approach to building robots based on the tension networks of tensegrity structures, which have no rigid connections between elements. This repository contains some of my research and relevant software tools related to optimal control of tensegrity robots!

![rolling tensegrity](/Images/rolling_tensegrity.png)

This repository is a simulation framework that enables rapid prototyping of tensegrity designs, actuation and sensing topologies, control policies, and state observers. To this end, the code in this repository follows object-oriented design principles, where the System Plant, Controllers, and Observers are all independently encapsulated to enable quick plug-and-play experimentation of any system component. As an example, you may choose to implement Model Predictive Control with full-state information or experiment with iLQR and an Unscented Kalman Filter with noisy sensor data. For a published work example of this simulation framework, check out this [IROS 2018 paper](https://ieeexplore.ieee.org/document/8594401) which uses some of the tools in this repository to design a contextual neural network policy using optimal state-action trajectories generated using Model Predictive Control!

![top-level design](/Images/System_Design.png)

## Tensegrity Dynamics

This work also features automatic generation of dynamics equations for any Type-1 tensegrity system (connection nodes are connected to at most *one* rigid body), based on a simplified pointmass dynamics formulation. While this trades off simulation accuracy, this approach enables rapid prototyping of different mechanical designs and sensor/actuator configurations to understand broad patterns of system feasibility, depending on design requirements. 

![top-level design](/Images/dynamics.png)

## Getting Started

## Running Examples

## Authors
* **Brian Cera** - *Initial Work*

