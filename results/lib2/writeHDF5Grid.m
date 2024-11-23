function writeHDF5Grid(fname,ngrd,ND,N,X2,IB2,flg)

if(flg)
	delete(fname);
end

h5create(fname,'/Group001/X',[N(1),N(2)])
h5create(fname,'/Group001/Y',[N(1),N(2)])
h5create(fname,'/Group001/IBLANK',[N(1),N(2)])

h5write(fname,'/Group001/X',squeeze(X2(1,:,:,1,1)))
h5write(fname,'/Group001/Y',squeeze(X2(1,:,:,1,2)))
h5write(fname,'/Group001/IBLANK',squeeze(IB2(1,:,:)))

h5writeatt(fname,'/','HEADER',[0.0])
h5writeatt(fname,'/','numberOfGrids',[1])

h5writeatt(fname,'/Group001','HEADER',[0.0, 0.0, 0.0, 0.0])
h5writeatt(fname,'/Group001','gridSize',[N(1),N(2)])
h5writeatt(fname,'/Group001','numberOfAuxVars',[0])
h5writeatt(fname,'/Group001','numberOfDependentVars',[0])
h5writeatt(fname,'/Group001','numberOfGradPhiVars',[0])
h5writeatt(fname,'/Group001','numberOfPhiVars',[0])
h5writeatt(fname,'/Group001','useIB',[1])
  
return
