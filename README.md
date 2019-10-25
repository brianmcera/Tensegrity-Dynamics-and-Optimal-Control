# Tensegrity-Dynamics-and-Optimal-Control
Tensegrity Robots are a biologically-inspired approach to building robots based on the tension networks of tensegrity structures, which have no rigid connections between elements. This repository contains some of my research and relevant software tools related to optimal control of tensegrity robots!

![rolling tensegrity](/Images/rolling_tensegrity.png)

In this repository, I'm working to build a simulation framework that will enable rapid prototyping of tensegrity designs, actuation and sensing topologies, control policies, and state observers. To this end, the code in this repository follows an object-oriented design principle, where the system Plant, Controllers, and Observers are all independently encapsulated to enable quick plug-and-play experimentation of any system component. As an example, you may choose to implement Model Predictive Control with full-state information or experiment with iLQR and an Unscented Kalman Filter with noisy sensor data. 

![top-level design](/Images/System_Design.png)

## Tensegrity Dynamics

This work also features automatic generation of dynamics equations for any Type-1 tensegrity system (connection nodes are connected to at most *one* rigid body), based on a simplified pointmass dynamics formulation.

![top-level design](/Images/Dynamics.png)

## Getting Started

## Running Examples

## Authors
* **Brian Cera** - *Initial Work*

