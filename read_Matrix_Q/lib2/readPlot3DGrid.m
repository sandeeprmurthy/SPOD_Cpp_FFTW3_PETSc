function [N,X2,IB] = readPlot3DGrid(fname,ND)
% [N,X2,IBLANK] = readPlot3DGrid(fname,ND)

sizeof_dbl = 8;
sizeof_int = 4;

  %%%% now read PLOT3D grid %%%%
  fid = fopen(fname,'rb','ieee-be');

  % number of meshes
  record = fread(fid,1,'int32');
  ngrd   = fread(fid,1,'int32');
  record = fread(fid,1,'int32');
  N = zeros(ngrd,ND);

  % mesh size
  record = fread(fid,1,'int32');
  for n = 1:ngrd
    N(n,:) = fread(fid,ND,'int32');
  end
  record = fread(fid,1,'int32');
  
  if (ND == 2) 

    % PLOT3D grid
    X2 = repmat(0,[ngrd max(N(:,1)) max(N(:,2)) ND]);
    IB = repmat(0,[ngrd max(N(:,1)) max(N(:,2))]);

    % read the data
    for n = 1:ngrd
      record = fread(fid,1,'int32');
      for k = 1:ND
        for j = 1:N(n,2)
          for i = 1:N(n,1)
            X2(n,i,j,k) = fread(fid,1,'double');
          end
        end
      end
      for j = 1:N(n,2)
        for i = 1:N(n,1)
          IBLANK(n,i,j) = fread(fid,1,'int32');
        end
      end
      record = fread(fid,1,'int32');
    end
  
  elseif (ND == 3) 
    
    % PLOT3D grid
    X2 = repmat(0,[ngrd max(N(:,1)) max(N(:,2)) 1 ND]);
    IB = repmat(0,[ngrd max(N(:,1)) max(N(:,2)) 1]);
    
    if (sum(N(:,3)) == size(N,1)) % 2D embedded in 3D

      % PLOT3D grid
      X2 = repmat(0,[ngrd max(N(:,1)) max(N(:,2)) 1 ND]);
      IB = repmat(1,[ngrd max(N(:,1)) max(N(:,2)) 1]);

      % read the data
      for ng = 1:ngrd
        record = fread(fid,1,'int32');        
        ftmp   = fread(fid,N(ng,1)*N(ng,2)*3,'double');
        X2(ng,1:N(ng,1),1:N(ng,2),1,:) = reshape(ftmp,[N(ng,1) N(ng,2) 3]);
        % check record to see if IB is included
        if (record ~= N(ng,1)*N(ng,2)*3*sizeof_dbl) 
          itmp = fread(fid,N(ng,1)*N(ng,2),'int32');
          IB(ng,1:N(ng,1),1:N(ng,2),1) = reshape(itmp,[N(ng,1) N(ng,2)]);
        end
        record = fread(fid,1,'int32');
      end
      
    end
  end


  fclose(fid);
  
  % X2 = squeeze(X2);

  %fprintf(1,'Read file %s.\n',fname);
  
return
