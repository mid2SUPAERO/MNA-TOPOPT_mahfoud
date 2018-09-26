%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function toprip(nelx,nely,volfrac,penal,rmin);
nelx=160;
nely=40;
volfrac=0.45;
penal=3.0;
rmin=1.2;
% INITIALIZE
x(1:nely,1:nelx) = volfrac;


for ely1 = 1:nely
    for elx1 = 1:nelx
         passivewhite(ely1,elx1) = 0;
         passiveblack(ely1,elx1) = 0;
     end
end


d = (3*nelx)/8;
R = 6*nely;                                   % radii
Rint1 = nely/4;
Rint2 = nely/4;
Rint3 = nely/6;
a = floor((sqrt((R)^2-(d)^2))-R+nely);        % edges
b = floor((sqrt((R)^2-(nelx-d)^2))-R+nely);
f = ((nelx/4)+Rint1+(((nelx/4)-Rint1-Rint2)/2)-(nelx/80)+1);  %white square dimensions
g = ((nelx/4)+Rint1+(((nelx/4)-Rint1-Rint2)/2)+(nelx/80));

 
% for ely1 = 1:nely            % black round
%   for elx1 = 1:nelx
%       if sqrt((ely1-R)^2+(elx1-d)^2) <= R
%           passiveblack(ely1,elx1) = 1;
%           x(ely1,elx1) = 1;
%       end
%   end
% end


for ely2 = 1:nely            % white external round
   for elx2 = 1:nelx
       if sqrt((ely2-R)^2+(elx2-d)^2) > R
           passivewhite(ely2,elx2) = 1;
           x(ely2,elx2) = 0.001;
       end
   end
end

for ely3 = 1:nely            % white internal round 1
   for elx3 = 1:nelx
       if sqrt((ely3-((3*nely)/5))^2+(elx3-(nelx/4))^2) <= Rint1
           passivewhite(ely3,elx3) = 1;
           x(ely3,elx3) = 0.001;
       end
   end
end

for ely4 = 1:nely            % white internal round 2
   for elx4 = 1:nelx
       if sqrt((ely4-((3*nely)/5))^2+(elx4-(nelx/2))^2) <= Rint2
           passivewhite(ely4,elx4) = 1;
           x(ely4,elx4) = 0.001;
       end
   end
end

for ely5 = 1:nely            % white internal round 3
   for elx5 = 1:nelx
       if sqrt((ely5-((3*nely)/5))^2+(elx5-((3*nelx)/4))^2) <= Rint3
           passivewhite(ely5,elx5) = 1;
           x(ely5,elx5) = 0.001;
       end
   end
end

for ely6 = ((nely/10)+1) : (9*nely)/10             % white square
     for elx6 = f : g
         passivewhite(ely6,elx6) = 1;
         x(ely6,elx6) = 0.001;
     end
end

%for elx7 = 1                               % forward edge
%     for ely7 = (nely-a+1) : nely
%        passiveblack(ely7,elx7) = 1;
%         x(ely7,elx7) = 1;
%     end
%end

%for elx8 = nelx                            % rear edge
%     for ely8 = (nely-b+1) : nely
%         passiveblack(ely8,elx8) = 1;
%        x(ely8,elx8) = 1;
%     end
%end


loop = 0;
change = 1.;
% START ITERATION
while change > 0.01
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal,d,a,b,f,g);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely;
      n2 = (nely+1)* elx   +ely;
      dc(ely,elx) = 0.;
      for i = 1:2
         Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],i);
         c = c + x(ely,elx)^penal*Ue'*KE*Ue;
   %     Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
   %     c = c + x(ely,elx)^penal*Ue'*KE*Ue;
        dc(ely,elx) = dc(ely,elx) -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; % aggiungi  dc(ely,elx)
      end
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc,passiveblack,passivewhite);       % passive
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
  %x1=x;             % symmetry
  %x2=flipdim(x,1);
  %x=(x1+x2)/2;
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end
print -dbmp16m result-toprip.bmp
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc,passiveblack,passivewhite)       % passive
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));

  xnew(find(passiveblack)) = 1;               % only with submitted regions
  xnew(find(passivewhite)) = 0.001;

  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0;
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal,d,a,b,f,g)
[KE] = lk;
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),2); U = zeros(2*(nely+1)*(nelx+1),2);   % cambia in 2
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely;
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS

dinx1 = 2* (                    nely-a +1) -1;   % above
diny1 = 2* (                    nely-a +1)   ;
dinx2 = 2* ((nely+1) * (d)             +1) -1;
diny2 = 2* ((nely+1) * (d)             +1)   ;
dinx3 = 2* ((nely+1) * (2*d)  + nely-a +1) -1;
diny3 = 2* ((nely+1) * (2*d)  + nely-a +1)   ;
dinx4 = 2* ((nely+1) * (nelx) + nely-b +1) -1;
diny4 = 2* ((nely+1) * (nelx) + nely-b +1)   ;

dout1 = 2* ((nely+1) * f   + (nely/10) +1)   ;
dout3 = 2* ((nely+1) * g   + (nely/10) +1)   ;


dinx5 = 2* ((nely+1)            ) -1;            % below
diny5 = 2* ((nely+1)            )   ;
dinx6 = 2* ((nely+1) * (d+1)    ) -1;
diny6 = 2* ((nely+1) * (d+1)    )   ;
dinx7 = 2* ((nely+1) * (2*d+1)  ) -1;
diny7 = 2* ((nely+1) * (2*d+1)  )   ;
dinx8 = 2* ((nely+1) * (nelx+1) ) -1;
diny8 = 2* ((nely+1) * (nelx+1) )   ;

dout2 = 2* ((nely+1) * f + (9*nely)/10);
dout4 = 2* ((nely+1) * g + (9*nely)/10);


F(dinx1,1) =  1;   % above
F(diny1,1) =  1;
F(dinx2,1) =  1;
F(diny2,1) =  1;
F(dinx3,1) =  1;
F(diny3,1) =  1;
F(dinx4,1) =  1;
F(diny4,1) =  1;

F(dout1,2) =  1;
F(dout3,2) =  1;


F(dinx5,1) = -1;   % below
F(diny5,1) = -1;
F(dinx6,1) = -1;
F(diny6,1) = -1;
F(dinx7,1) = -1;
F(diny7,1) = -1;
F(dinx8,1) = -1;
F(diny8,1) = -1;

F(dout2,2) = -1;
F(dout4,2) = -1;


fixeddofs1 = (2*(nely-a+1+1)-1):(2*(nely+1-1)-1);
fixeddofs2 = 2*(nely-a+1+1):2*(nely+1-1);
%fixeddofs1 = 2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1);
%fixeddofs2 = (2*(nely+1)-1):2*(nely+1):(2*(nely+1)*(nelx+1)-1);
fixeddofs = union(fixeddofs1,fixeddofs2);

alldofs     = (1:2*(nely+1)*(nelx+1));
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 1.0;
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
