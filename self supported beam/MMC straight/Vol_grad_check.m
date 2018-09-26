%% ************* volume gradient check by finite differences *********** %%
clear ; close all ;
load grad_data ;
gradv = zeros(nn,1) ;
X = xy00 ; % design vector to be pertrubated
delta = 0.05 ;
for j = 1:N
    for ii = 1:Var_num
        k = ii+Var_num*(j-1) ;
        pert = zeros(nn,1) ;
        pert(k) = delta ;
        X1 = X-pert ; X2 = X+pert ;
        tmpPhiD1=tPhi(X1(Var_num*j-Var_num+1:Var_num*j),LSgrid.x,LSgrid.y,p);
        tempPhi_max1=tmpPhiD1;
        for ik=1:j-1
            tempPhi_max1=max(tempPhi_max1,Phi{ik});
        end
        for ik=j+1:N
            tempPhi_max1=max(tempPhi_max1,Phi{ik});
        end
        Phi_KS1=zeros(size(Phi{1}));
        for i=1:j-1
            Phi_KS1=Phi_KS1+1/N*exp(P*(Phi{i}-tempPhi_max1));
        end
        for i=j+1:N
            Phi_KS1=Phi_KS1+1/N*exp(P*(Phi{i}-tempPhi_max1));
        end
        Phi_KS1=Phi_KS1+1/N*exp(P*(tmpPhiD1-tempPhi_max1));
        Phi_KS1=tempPhi_max1+1/P*log(Phi_KS1);
        tmpPhiD2=tPhi(X2(Var_num*j-Var_num+1:Var_num*j),LSgrid.x,LSgrid.y,p);
        tempPhi_max2=tmpPhiD2;
        for ik=1:j-1
            tempPhi_max2=max(tempPhi_max2,Phi{ik});
        end
        for ik=j+1:N
            tempPhi_max2=max(tempPhi_max2,Phi{ik});
        end
        Phi_KS2=zeros(size(Phi{1}));
        for i=1:j-1
            Phi_KS2=Phi_KS2+1/N*exp(P*(Phi{i}-tempPhi_max2));
        end
        for i=j+1:N
            Phi_KS2=Phi_KS2+1/N*exp(P*(Phi{i}-tempPhi_max2));
        end
        Phi_KS2=Phi_KS2+1/N*exp(P*(tmpPhiD2-tempPhi_max2));
        Phi_KS2=tempPhi_max2+1/P*log(Phi_KS2);
        H1=Heaviside(Phi_KS1,alpha,nelx,nely,epsilon);
        H2=Heaviside(Phi_KS2,alpha,nelx,nely,epsilon);
        den1 = sum(H1(EleNodesID),2)/4 ; den2 = sum(H2(EleNodesID),2)/4 ;
        V1 = (sum(den1)/(nelx*nely)-volfrac)/volfrac*100 ;
        V2 = (sum(den2)/(nelx*nely)-volfrac)/volfrac*100 ;
        gradv(k) = (V2-V1)/(2*delta) ;
    end
end
dfdx1 = sum(dden)/nelx/nely/volfrac*100 ;
error = abs(dfdx1'-gradv)./abs(dfdx1') ;
figure ;
plot(error) ;
title('volume gradient error') ;
xlabel('rank of design variable') ;
ylabel('relative error value') ;