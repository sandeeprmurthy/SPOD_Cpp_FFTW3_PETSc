clear all;
clc;
close all;

addpath(genpath('./lib'));
addpath(genpath('./lib2'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Read the matrix Q ...')
Q = PetscBinaryRead('../Q_matrix.dat', 'indices','int32','precision','float64','complex',true);
Q=full(Q);
Q(1:5,1:5)













