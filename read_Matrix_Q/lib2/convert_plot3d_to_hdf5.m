% main.m
clear all;
clc;
close all;

% Read the grid and solution file
fname_grid = './cd_nozzle_plot3d_file.xyz'; 

% Read the grid
disp('Reading the Plot3D grid ...')
[N,X,ib] = readPlot3DGrid(fname_grid,3);

% Modify iblanking
inner_wall_index = 159;
outer_wall_index = 208;
top_wall_index = 1752;
gp = 9;
ib(1,inner_wall_index+gp+1:outer_wall_index-gp-1,1:top_wall_index-gp-1) = 0;

% Read the grid
disp('Writing the HDF5 grid ...')
writeHDF5Grid('./cd_nozzle_hdf5_file.h5',1,3,N,X,ib);
