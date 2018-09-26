clear ; close all ;
p = 6 ; nelx = 160 ; nely = 40 ; DW = 4 ; DH = 1 ;
EW = DW/nelx ; EH = DH/nely ; 
[meshx,meshy] = meshgrid(EW*(0:nelx),EH*(0:nely)) ;
x = [2 0.5 0 1 0.2 0.5] ;
phi = TDF_cc(x,meshx(:),meshy(:),p) ;
figure ;
contourf(meshx,meshy,reshape(phi,nely+1,nelx+1),[0,0]) ;
axis equal; axis([0 DW 0 DH]);  drawnow; title('component plot') ;