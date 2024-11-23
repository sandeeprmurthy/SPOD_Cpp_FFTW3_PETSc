function [status] = writePlot3DGrid(fname,ngrd,ND,N,X,IB)
% [status] = writePlot3DGrid(fname,ngrd,ND,N,X,IB)
% fname = filename (string)
% ngrd  = number of grids (int)
% ND    = dimension (int, should always equal 3)
% N     = array of grid dimensions N(ng,:) = [N1 N2 N3] for grid ng
% X     = array of grid coordinates X(ng,:,:,:,:) for grid ng
% IB    = array of iblank values IB(ng,:,:,:) for grid ng

  % assume 32 bit
  sizeof_int = 4; % 4 byte integers
  sizeof_dbl = 8; % 8 byte doubles

  %%%% now read PLOT3D grid %%%%
  fid = fopen(fname,'wb','ieee-be');

  % number of meshes
  record = sizeof_int;
  fwrite(fid,record,'int32');
  fwrite(fid,ngrd,'int32');
  fwrite(fid,record,'int32');
  
    % mesh size
  record = ngrd*ND*sizeof_int;
  fwrite(fid,record,'int32');
  for n = 1:ngrd
    fwrite(fid,N(n,1:ND),'int32');
  end
  fwrite(fid,record,'int32');

  if (ND ~= 3) 
      error('ND ~= 3 is not supported.');
  end

  if (ND == 3)
      
    if (size(N,2) ~= 3)
        error('Inconsistent input: size(N,2) ~= 3.');
    end
      
    % write the data
    for n = 1:ngrd
        
      % write the grid coordinates
      record = prod(N(n,1:3))*ND*sizeof_dbl + prod(N(n,1:3))*sizeof_int;
      fwrite(fid,record,'int32');
      fwrite(fid,X(n,1:N(n,1),1:N(n,2),1:N(n,3),1:3),'double');
      fwrite(fid,IB(n,1:N(n,1),1:N(n,2),1:N(n,3)),'int32');
      fwrite(fid,record,'int32');

    end
      
  end


  fclose(fid);

  fprintf(1,'Wrote file %s.\n',fname);
  
return
