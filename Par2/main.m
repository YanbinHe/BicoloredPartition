% reproducible code for paper 'Dynamic Bi-Colored Graph Partitioning'

% run simulation method
% This code should be loaded first 
clc
clear
% this is the data for test and used in the paper. You can generate your
% own data/use your own graph to test this method

load test.mat% predefined graph

Wini = Adjini; 
Lini = diag(Wini*ones(size(Wini,1),1)) - Wini;
newStatus = newStatusini;
UserStatusChanging;
