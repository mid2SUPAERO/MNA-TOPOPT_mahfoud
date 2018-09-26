%% FUNCTION prepcoarse - PREPARE MG PROLONGATION OPERATOR
function [Pu] = Prepcoarse(ney,nex)
% Assemble state variable prolongation
maxnum = nex*ney*20;
iP = zeros(maxnum,1); jP = zeros(maxnum,1); sP = zeros(maxnum,1);
nexc = nex/2; neyc = ney/2;
% weights for fixed distances to neighbors on a structured grid 
vals = [1,0.5,0.25];
cc = 0;
for nx = 1:nexc+1
    for ny = 1:neyc+1
        col = (nx-1)*(neyc+1)+ny;
        % Coordinate on fine grid
        nx1 = nx*2 - 1; ny1 = ny*2 - 1;
        % Loop over fine nodes within the rectangular domain
        for k = max(nx1-1,1):min(nx1+1,nex+1)
            for l = max(ny1-1,1):min(ny1+1,ney+1)
                row = (k-1)*(ney+1)+l;
                % Based on squared dist assign weights: 1.0 0.5 0.25
                ind = 1+((nx1-k)^2+(ny1-l)^2);
                cc = cc+1; iP(cc) = 2*row-1; jP(cc) = 2*col-1; sP(cc) = vals(ind);
                cc = cc+1; iP(cc) = 2*row; jP(cc) = 2*col; sP(cc) = vals(ind);
            end
        end
    end
end
% Assemble matrices
Pu = sparse(iP(1:cc),jP(1:cc),sP(1:cc));
end
