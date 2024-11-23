function [data,Q2] = readPlot3DSoln(fname,ND)
% [data,Q2] = readPlot3DSoln(fname,ND)
% data(1) = iter;
% data(4) = time;

  %%%% now read PLOT3D soln %%%%
  fid = fopen(fname,'rb','ieee-be');

  % number of grids
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
    
    Q2 = repmat(0,[ngrd max(N(:,1)) max(N(:,2)) (ND+2)]);

    % flow data
    for n = 1:ngrd
  
      % flow header
      record = fread(fid,1,'int32');
      data   = fread(fid,4,'double'); iter = data(1); t = data(4);
      record = fread(fid,1,'int32');

      record = fread(fid,1,'int32');
      for k = 1:(ND+2)
        for j = 1:N(n,2)
          for i = 1:N(n,1)
            Q2(n,i,j,k) = fread(fid,1,'double');
          end
        end
      end
      record = fread(fid,1,'int32');
    end
  end

  if (ND == 3) 
    
    Q2 = repmat(0,[ngrd max(N(:,1)) max(N(:,2)) max(N(:,3)) (ND+2)]);

    % flow data
    for n = 1:ngrd
  
      % flow header
      record = fread(fid,1,'int32');
      data   = fread(fid,4,'double'); iter = data(1); t = data(4);
      record = fread(fid,1,'int32');

      record = fread(fid,1,'int32');
      rbuf   = fread(fid,N(n,1)*N(n,2)*N(n,3)*(ND+2),'double');
      Q2(n,1:N(n,1),1:N(n,2),1:N(n,3),:) = reshape(rbuf,[N(n,1) N(n,2) N(n,3) (ND+2)]);
      record = fread(fid,1,'int32');
    end
  end

  
  % done
  fclose(fid);

  % Q2 = squeeze(Q2);
  
  fprintf(1,'Read file %s.\n',fname);


return