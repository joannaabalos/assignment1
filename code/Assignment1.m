%% Monte-Carlo Modeling of Electron Transport
% Assignment 1 - Joanna Abalos 100962263

close all
clear
clc

Assignment1_1
Assignment1_2
Assignment1_3

% In this assignment, 10 000 particles are modelled to calculate
% temperatures, make models and observations using Monte-Carlo modeling. 7
% particles are plotted to observe their trajectories.

%% 1 Electron Modelling

% Figure 1 displays a subset of the particle trajectories travelling at the 
% thermal velocity. The particles that reach the top and bottom boundary 
% have their velocities reversed.

% Figure 2 shows that the semiconductor temperature remains constant since 
% all particles are travelling at the same velocity.

%% 2 Collisions with Mean Free Path (MFP)

% Figure 3 displays the hystersis plots of the particle velocities. The X
% and Y components of velocity are a normal distribution scaled by the
% thermal velocity (Vth). The resulting overall velocity of each particle
% results in a Maxwell-Boltzmann distribution.

% Figure 4 displays a subset of the particle trajectories travelling at the
% velocities defined by the Maxwell-Boltzmann distribution. After
% calculating the probability of scattering and applying them to each
% particle, the scattered particles are assigned new velocities as defined 
% by the Maxwell-Boltzmann distribution.

% Figure 5 shows that the average semiconductor temperature over time
% averages at around 300K. This is because the particles are travelling at 
% an average velocity of Vth.

%% 3 Enhancements

% Figure 6 displays a subset of the particle trajectories that scatter 
% similar to Figure 4 except with rectangle boundaries at which the 
% particles bounce.

% The Electron Density Map in Figure 7 shows that the particles are not 
% present within the boxes. It appears that a few particles did penetrate 
% the edges slightly but box region should have no particle penetration.

% The Temperature Map in Figure 7 shows that the average overall
% temperature of the region is around 300K. There is no temperature
% calculated in the areas where there are no moving particles (within the
% boxes).





