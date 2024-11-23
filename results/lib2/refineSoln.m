function [ qnew ] = refineSoln( Nold, Xold, Qold, Nnew, Xnew )
%qnew = refineSoln(Nold, xold, qold, Nnew, xnew)
%   Given a solution qold defined on xold, interpolate a new
%   solution qnew on the grid xnew

  % array sizes
  size_xold = size(Xold);
  size_qold = size(Qold);
  size_xnew = size(Xnew);
  
  % number of grids -- interpolate grid-by-grid
  ngrd = size_xold(1);
  
  % error check array sizes
  if (size_xold(1:4) ~= size_qold(1:4))
      error('xold and qold have incompatible sizes.');
  end
  
  % error check on number of grids
  if (size_xnew(1) ~= ngrd)
      error('xnew and xold have different numbers of grids');
  end
  
  % set qnew
  qnew = repmat(0,[size_xnew(1:4) size_qold(end)]);
  
  % now interpolate
  for n = 1:ngrd
      vold1 = 1:Nold(n,1);
      vold2 = 1:Nold(n,2);
      vold3 = 1:Nold(n,3);
      vnew1 = 1:Nnew(n,1);
      vnew2 = 1:Nnew(n,2);      
      vnew3 = 1:Nnew(n,3);      
      xold = squeeze(Xold(n,vold1,vold2,vold3,1));
      yold = squeeze(Xold(n,vold1,vold2,vold3,2));
      zold = squeeze(Xold(n,vold1,vold2,vold3,3));
      xnew = squeeze(Xnew(n,vnew1,vnew2,vnew3,1));
      ynew = squeeze(Xnew(n,vnew1,vnew2,vnew3,2));          
      znew = squeeze(Xnew(n,vnew1,vnew2,vnew3,3));
      for l = 1:size_qold(5)
          qold = squeeze(Qold(n,vold1,vold2,vold3,l));
          if (Nold(n,3) ~= 1)
    		  qnew(n,vnew1,vnew2,vnew3,l) = ...
                griddata(xold.',yold.',zold.',qold.',xnew.',ynew.',znew.').';
          else
              %qnew(n,vnew1,vnew2,vnew3,l) = ...
              %  griddata(xold.',yold.',qold.',xnew.',ynew.').';
              F = scatteredInterpolant(xold(:), yold(:), qold(:), 'natural', 'linear');
              qnew(n,vnew1,vnew2,vnew3,l) = F(xnew,ynew);
          end
      end

  end
  
end

