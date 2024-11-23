clear all
clc


A = PetscBinaryRead('A.dat','complex',true);
spy(A)
% B = PetscBinaryRead('B_mod_corners_ss.dat','complex',true);

% break
% 
% i_cplx = sqrt(-1);
% 
% sigma = 0.2 - 0.8*i_cplx;
% 
% [V,D] = eigs(A,B,20,sigma);

%plot(real(V(5:5:1135)))
