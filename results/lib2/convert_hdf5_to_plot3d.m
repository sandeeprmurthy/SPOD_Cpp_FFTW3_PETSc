% main.m
clear all;
clc;
close all;

% Read the grid and solution file
fname_grid = '../hdf5_file.h5'; 

% Read the grid
disp('Reading the HDF5 grid ...')
[N,X,ib] = readHDF5Grid(fname_grid,3);

% Read the grid
disp('Writing the plt3d grid ...')
writePlot3DGrid('plot3d_file.xyz',1,3,N,X,ib);
