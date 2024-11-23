function [Q] = readHDF5Soln(fname,ND,test_name)

disp('readHDF5Soln: Reading the conserved data from the hdf5 file ...')
rho  = h5read(fname,['/PlasCom2/Simulation/',test_name,'/grid1/rho']);
rhoU = h5read(fname,['/PlasCom2/Simulation/',test_name,'/grid1/rhoV-1']);
rhoV = h5read(fname,['/PlasCom2/Simulation/',test_name,'/grid1/rhoV-2']);
rhoW = 0*rhoV;
rhoE = h5read(fname,['/PlasCom2/Simulation/',test_name,'/grid1/rhoE']);

N = [size(rho,1),size(rho,2),1];

disp('readHDF5Soln: Constructing an object that can be supplied to the refining algorithm ...')
Q = zeros(1,N(1),N(2),N(3),ND+2);
Q(1,:,:,1,1) = rho;
Q(1,:,:,1,2) = rhoU;
Q(1,:,:,1,3) = rhoV;
Q(1,:,:,1,4) = rhoW;
Q(1,:,:,1,5) = rhoE;

return
