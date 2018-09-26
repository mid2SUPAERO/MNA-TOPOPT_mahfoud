%% ************ density field and it's gradient calculation ************ %%
% this function helps implementing MNA method in the 3D engine pylon codes,
% the inputs are design variables (of MNA) and outputs are density field
% and it's gradient.
function [den,grad_den] = mna_density(FEM_structure,X)
trans = 500 ; load V ; %#ok<LOAD> 
n_var = length(X) ;
n_c = n_var/9 ;
den = zeros(FEM_structure.n,1) ;
grad_den = zeros(FEM_structure.n,n_var) ;
for s = 1:n_c % loop on components
    ss = 9*(s-1)+1 ;
    lc = 2*sqrt(X(ss+3)^2+X(ss+4)^2+X(ss+5)^2) ;
    neigh_id = [] ;
    neigh = [] ;
    neigh_g =[] ;
    % component's neighborhood calculation
    for i = 1:FEM_structure.n % loop on elements
        xc = 1/8*sum([FEM_structure.x1(i),FEM_structure.x2(i),...
            FEM_structure.x3(i),FEM_structure.x4(i),FEM_structure.x5(i),...
            FEM_structure.x6(i),FEM_structure.x7(i),FEM_structure.x8(i)]) ;
        yc = 1/8*sum([FEM_structure.y1(i),FEM_structure.y2(i),...
            FEM_structure.y3(i),FEM_structure.y4(i),FEM_structure.y5(i),...
            FEM_structure.y6(i),FEM_structure.y7(i),FEM_structure.y8(i)]) ;
        zc = 1/8*sum([FEM_structure.z1(i),FEM_structure.z2(i),...
            FEM_structure.z3(i),FEM_structure.z4(i),FEM_structure.z5(i),...
            FEM_structure.z6(i),FEM_structure.z7(i),FEM_structure.z8(i)]) ;
        if norm([xc-X(ss);yc-X(ss+1);zc-X(ss+2)])<lc
            [ksi,eta,mu] = basis_change(xc,yc,zc,X(ss:ss+8)) ;
            if ksi<1 && eta<1 && mu<1
                neigh_id = [neigh_id,i] ; %#ok<AGROW>
                neigh = [neigh,[ksi;eta;mu]] ; %#ok<AGROW>
                neigh_g = [neigh_g,[xc;yc;zc]] ; %#ok<AGROW>
            end
        end
    end
    if ~isempty(neigh)
        ln = length(neigh_id) ;
        wx = zeros(ln,1) ; wy = zeros(ln,1) ; wz = zeros(ln,1) ;
        dwx_dksi = zeros(ln,1) ;dwy_deta = zeros(ln,1) ;dwz_dmu = zeros(ln,1) ;
        p1 = neigh(1,:)<0.5 ; p2 = neigh(2,:)<0.5 ; p3 = neigh(3,:)<0.5 ;
        % shape functions and their derivatives
        if FEM_structure.outit < trans
            wx(p1) = 1-6*neigh(1,p1).^2+6*neigh(1,p1).^3 ;
            wx(~p1) = 2-6*neigh(1,~p1)+6*neigh(1,~p1).^2-2*neigh(1,~p1).^3;
            wy(p2) = 1-6*neigh(2,p2).^2+6*neigh(2,p2).^3 ;
            wy(~p2) = 2-6*neigh(2,~p2)+6*neigh(2,~p2).^2-2*neigh(2,~p2).^3;
            wz(p3) = 1-6*neigh(3,p3).^2+6*neigh(3,p3).^3 ;
            wz(~p3) = 2-6*neigh(3,~p3)+6*neigh(3,~p3).^2-2*neigh(3,~p3).^3;
        else
            wx(p1) = 1-6*neigh(1,p1).^8+6*neigh(1,p1).^12 ;
            wx(~p1) = V(1)+V(2)*neigh(1,~p1).^2+V(3)*neigh(1,~p1).^4+...
                V(4)*neigh(1,~p1).^6+V(5)*neigh(1,~p1).^8+V(6)*neigh(1,~p1).^12 ;
            wy(p2) = 1-6*neigh(2,p2).^8+6*neigh(2,p2).^12 ;
            wy(~p2) = V(1)+V(2)*neigh(1,~p2).^2+V(3)*neigh(1,~p2).^4+...
                V(4)*neigh(1,~p2).^6+V(5)*neigh(1,~p2).^8+V(6)*neigh(1,~p2).^12 ;
            wz(p3) = 1-6*neigh(2,p3).^8+6*neigh(2,p3).^12 ;
            wz(~p3) = V(1)+V(2)*neigh(1,~p3).^2+V(3)*neigh(1,~p3).^4+...
                V(4)*neigh(1,~p3).^6+V(5)*neigh(1,~p3).^8+V(6)*neigh(1,~p3).^12 ;
        end
        W = wx.*wy.*wz ;
        if FEM_structure.outit < trans
            dwx_dksi(p1) = (-12*neigh(1,p1)+18*neigh(1,p1).^2) ;
            dwx_dksi(~p1) = (-6+12*neigh(1,~p1)-6*neigh(1,~p1).^2) ;
            dwy_deta(p2) = (-12*neigh(2,p2)+18*neigh(2,p2).^2) ;
            dwy_deta(~p2) = (-6+12*neigh(2,~p2)-6*neigh(2,~p2).^2) ;
            dwz_dmu(p3) = (-12*neigh(3,p3)+18*neigh(3,p3).^2) ;
            dwz_dmu(~p3) = (-6+12*neigh(3,~p3)-6*neigh(3,~p3).^2) ;
        else
            dwx_dksi(p1) = (-48*neigh(1,p1).^7+72*neigh(1,p1).^11) ;
            dwx_dksi(~p1) = 2*V(2)*neigh(1,~p1)+4*V(3)*neigh(1,~p1).^3+6*V(4)*neigh(1,~p1).^5+8*V(5)*neigh(1,~p1).^7+12*V(6)*neigh(1,~p1).^11 ;
            dwy_deta(p2) = (-48*neigh(2,p2).^7+72*neigh(2,p2).^11) ;
            dwy_deta(~p2) = 2*V(2)*neigh(1,~p2)+4*V(3)*neigh(1,~p2).^3+6*V(4)*neigh(1,~p2).^5+8*V(5)*neigh(1,~p2).^7+12*V(6)*neigh(1,~p2).^11 ;
            dwz_dmu(p2) = (-48*neigh(2,p2).^7+72*neigh(2,p2).^11) ;
            dwz_dmu(~p2) = 2*V(2)*neigh(1,~p2)+4*V(3)*neigh(1,~p2).^3+6*V(4)*neigh(1,~p2).^5+8*V(5)*neigh(1,~p2).^7+12*V(6)*neigh(1,~p2).^11 ;
        end
        % density and it's gradient computation
        den(neigh_id) = den(neigh_id) + W ;
        [grad_ksi,grad_eta,grad_mu] = grad_ksi_eta_mu(neigh_g,X(ss:ss+8)) ;
        grad_den(neigh_id,ss:ss+8) = repmat(dwx_dksi.*wy.*wz,1,9).*grad_ksi ...
            + repmat(wx.*dwy_deta.*wz,1,9).*grad_eta + ...
            repmat(wx.*wy.*dwz_dmu,1,9).*grad_mu ;
    end
end
% asymptotic density
grad_den(den<1,:) = diag((1+12*den(den<1).^2-28*den(den<1).^3+...
    15*den(den<1).^4))*grad_den(den<1,:);
grad_den(den>=1,:) = 0;
den = den+4*den.^3-7*den.^4+3*den.^5;
den = sparse(min(den,1));
grad_den = sparse(grad_den);
end