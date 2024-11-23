function [ Nnew, Xnew, IBnew ] = refineGrid( Nold, Xold, factor )
%[ Nnew, Xnew, IBnew ] = refineGrid( Nold, Xold, factor )
%   Produce a refined mesh with refinement factor = [f1 f2 f3];

    % number of grids
    ngrd = size(Xold,1);

    % refined grid-by-grid
    for n = 1:ngrd

        % size of this grid
        v1 = 1:Nold(n,1);
        v2 = 1:Nold(n,2);
        v3 = 1:Nold(n,3);

        % assume the grids are rectangles
        % find min / max of {x,y,z}
        xtmp = squeeze(Xold(n,v1,v2,v3,1));
        ytmp = squeeze(Xold(n,v1,v2,v3,2));
        ztmp = squeeze(Xold(n,v1,v2,v3,3));
        xmin = min(xtmp(:));
        xmax = max(xtmp(:));
        ymin = min(ytmp(:));
        ymax = max(ytmp(:));
        zmin = min(ztmp(:));
        zmax = max(ztmp(:));

        % size of new grid
        Nnew(n,1:3) = factor .* (Nold(n,1:3)-1) + 1;
        v1n = 1:Nnew(n,1);
        v2n = 1:Nnew(n,2);
        v3n = 1:Nnew(n,3);
        Xnew(n,v1n,v2n,v3n,1:3) = 0;
        IBnew(n,v1n,v2n,v3n) = 1;

        % linearly interpolate new spacing from old spacing
        % assume an orthogonal mesh
        xi_old = linspace(0,1,Nold(n,1)); eta_old = linspace(0,1,Nold(n,2));
        xi_new = linspace(0,1,Nnew(n,1)); eta_new = linspace(0,1,Nnew(n,2));
        
        [ETAold, XIold] = meshgrid(eta_old, xi_old);
        [ETAnew, XInew] = meshgrid(eta_new, xi_new);
        
        xx = interp2(XIold.', ETAold.', xtmp.', XInew.', ETAnew.').';
        yy = interp2(XIold.', ETAold.', ytmp.', XInew.', ETAnew.').';

        Xnew(n,v1n,v2n,v3n,1) = xx;
        Xnew(n,v1n,v2n,v3n,2) = yy;
        Xnew(n,v1n,v2n,v3n,3) = 0;

    end

end

