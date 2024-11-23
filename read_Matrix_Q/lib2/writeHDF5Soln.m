function writeHDF5Soln(fname_write,Q,test_name)

disp('writeHDF5Soln: Read the modified solution from input ...')
rho  = squeeze(Q(1,:,:,1,1));
rhoU = squeeze(Q(1,:,:,1,2));
rhoV = squeeze(Q(1,:,:,1,3));
rhoW = squeeze(Q(1,:,:,1,4));
rhoE = squeeze(Q(1,:,:,1,5));

disp('writeHDF5Soln: Write the modified solution back to the hdf5 file ...')
h5write(fname_write,['/PlasCom2/Simulation/',test_name,'/grid1/rho'],rho);
h5write(fname_write,['/PlasCom2/Simulation/',test_name,'/grid1/rhoV-1'],rhoU);
h5write(fname_write,['/PlasCom2/Simulation/',test_name,'/grid1/rhoV-2'],rhoV);
h5write(fname_write,['/PlasCom2/Simulation/',test_name,'/grid1/rhoE'],rhoE);
  
return
