function [N,X2,IB] = readHDF5GridMultiMesh(fname,gridName,ND)

X  = hdf5read(fname,['/PlasCom2/Geometry/cmgeom/',gridName,'/X']);
Y  = hdf5read(fname,['/PlasCom2/Geometry/cmgeom/',gridName,'/Y']);
%ib = hdf5read(fname,['/PlasCom2/Geometry/cmgeom/',gridName,'/IBLANK']);

N = [size(X,1),size(X,2),1];

X2 = zeros(1,N(1),N(2),N(3),ND);
X2(1,:,:,1,1) = X;
X2(1,:,:,1,2) = Y;
X2(1,:,:,1,3) = 0*Y;

IB = zeros(1,N(1),N(2));
%IB(1,:,:) = ib;

return
