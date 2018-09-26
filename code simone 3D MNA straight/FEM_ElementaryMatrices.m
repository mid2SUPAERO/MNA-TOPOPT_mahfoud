%                     [ FEM_ElementaryMatrices.m ]
%
%     WHAT IT DOES ?
%     Computes the FE elementary matrix
%
%  Main variables
%  --------------
%
%  NodesLocation:      Location (x, y, z) of each node
%  Material:           Material characteristics (rigidity constants, density)
%  InterpolationType:  Interpolation order along X, Y and Z
%  ElNodIndCard:       Map of the nodes location within the element
function [K,M] = FEM_ElementaryMatrices(NodesLocation,Material,InterpolationType,ElNodIndCard)    

% Solid characteristics
  C1 = Material(1);
  C2 = Material(2);
  C3 = Material(3);
  rho = Material(4);

% Number of dof per node
  DofNbPerNode = 3;

% [NodesNumber] Number of nodes in the element
  NodesNumber = length(NodesLocation(:,1));

% Gauss points and weights for each direction
  [GaussPointsX,GaussWeightsX] = FEM_GaussIntegration(InterpolationType(1),1);
  [GaussPointsY,GaussWeightsY] = FEM_GaussIntegration(InterpolationType(2),1);
  [GaussPointsZ,GaussWeightsZ] = FEM_GaussIntegration(InterpolationType(3),1);

% Quadrature number of points (Gauss optimal)
  QuadratureNbPoints = InterpolationType + 1;
    
% Matrices initialization
  K = zeros(NodesNumber*DofNbPerNode,NodesNumber*DofNbPerNode);
  M = zeros(NodesNumber*DofNbPerNode,NodesNumber*DofNbPerNode);

% Shape function and its derivative initialization
  N = zeros(NodesNumber,1); 
  dN = zeros(NodesNumber,3);
  D = zeros(NodesNumber,3);  
 
% Shape functions [N] and Shape functions derivative [Nxsi]       
  DirList = [1 2 3];
      
% Loop over each direction
  for d = 1:length(DirList)

    % Direction  
      Dir = DirList(d);

    % Integration location
      if (Dir ==1)
          GP = GaussPointsX;
      elseif (Dir == 2)
          GP = GaussPointsY;          
      elseif (Dir == 3)
          GP = GaussPointsZ;
      end

    % Number of nodes per direction  
      NNb = length(GP);

    % Reference domain points
      P = linspace(-1,1,NNb);    

      SF =  zeros(NNb,NNb);
      DF = zeros(NNb,NNb);
      
    % Loop over all integration points  
      for alpha = 1:NNb
          xsi = GP(alpha);

        % Loop over all nodes in this direction  
          for j1 = 1:NNb

            % A(1) = xsi(alpha) - xsi_1  
              A = (xsi*ones(1,NNb)-P); A(j1) = []; 

            % (xsi_2 - xsi_1) x (xsi_2 - xsi_3) x (xsi_2 - xsi_4) x ... 
              DD = P(j1)*ones(1,NNb)-P; DD = prod(DD(find(DD)));

            % SF = Shape fct at integration point alpha and node j1  
              SF(j1,alpha) = prod(A)/DD;

            % DF = Derivative (xsi - xsi_2)(xsi - xsi_3)... + (xsi - xsi_1)(xsi - xsi_3)... + ...  
              R = 0;
              for j2 = 1:(NNb-1)
                  B = A; B(j2) = 1; R = R + prod(B);                       
              end
              DF(j1,alpha) = R/DD; 
          end  
      end

    % Save the results  
      if (Dir == 1)
          SFX = SF(ElNodIndCard(:,1),:); DFX = DF(ElNodIndCard(:,1),:);         
      elseif (Dir == 2)
          SFY = SF(ElNodIndCard(:,2),:); DFY = DF(ElNodIndCard(:,2),:);                    

      elseif (Dir == 3)
          SFZ = SF(ElNodIndCard(:,3),:); DFZ = DF(ElNodIndCard(:,3),:);
      end
  end 
  
% Loop over all 3 directions in the parent domain (xsi, eta, dzeta)
  for Xintegral = 1:QuadratureNbPoints(1)
      for Yintegral = 1:QuadratureNbPoints(2)
          for Zintegral = 1:QuadratureNbPoints(3)

            % Gauss points weighting 
              W = GaussWeightsX(Xintegral)*GaussWeightsY(Yintegral)*GaussWeightsZ(Zintegral);

            % Shape functions and derivative in each direction 
              for node =1:NodesNumber
                  N(node) =    SFX(node,Xintegral)*SFY(node,Yintegral)*SFZ(node,Zintegral);                                            
                  dN(node,1) = DFX(node,Xintegral)*SFY(node,Yintegral)*SFZ(node,Zintegral);                      
                  dN(node,2) = SFX(node,Xintegral)*DFY(node,Yintegral)*SFZ(node,Zintegral);                      
                  dN(node,3) = SFX(node,Xintegral)*SFY(node,Yintegral)*DFZ(node,Zintegral);                                            
              end                 

            % Computation of the Xxsi, Yxsi, Zxsi, Xeta, Yeta, Zeta, Xdze, Ydze, Zdze quantities
              Xxsi = 0.0; Xeta = 0.0; Xdze = 0.0; 
              Yxsi = 0.0; Yeta = 0.0; Ydze = 0.0; 
              Zxsi = 0.0; Zeta = 0.0; Zdze = 0.0; 

              for node = 1:NodesNumber
                  Xxsi = Xxsi + dN(node,1)*NodesLocation(node,1);
                  Xeta = Xeta + dN(node,2)*NodesLocation(node,1); 
                  Xdze = Xdze + dN(node,3)*NodesLocation(node,1);
                  Yxsi = Yxsi + dN(node,1)*NodesLocation(node,2);
                  Yeta = Yeta + dN(node,2)*NodesLocation(node,2);
                  Ydze = Ydze + dN(node,3)*NodesLocation(node,2);
                  Zxsi = Zxsi + dN(node,1)*NodesLocation(node,3);
                  Zeta = Zeta + dN(node,2)*NodesLocation(node,3);
                  Zdze = Zdze + dN(node,3)*NodesLocation(node,3);
              end

            % Cofactors computation / Matrix transformation jacobian
              Cof11 = Yeta*Zdze - Zeta*Ydze;  
              Cof12 = Ydze*Zxsi - Zdze*Yxsi;
              Cof13 = Yxsi*Zeta - Zxsi*Yeta;
              Cof21 = Zeta*Xdze - Zdze*Xeta;
              Cof22 = Zdze*Xxsi - Zxsi*Xdze;
              Cof23 = Zxsi*Xeta - Zeta*Xxsi;
              Cof31 = Xeta*Ydze - Xdze*Yeta;
              Cof32 = Xdze*Yxsi - Xxsi*Ydze;
              Cof33 = Xxsi*Yeta - Xeta*Yxsi;

              DetJ = Xxsi*Cof11 + Xeta*Cof12 + Xdze*Cof13;

            % Shape functions complete derivation
              for node = 1:NodesNumber
                 D(node,1)= (dN(node,1)*Cof11 + dN(node,2)*Cof12 + dN(node,3)*Cof13) / DetJ;
                 D(node,2)= (dN(node,1)*Cof21 + dN(node,2)*Cof22 + dN(node,3)*Cof23) / DetJ;
                 D(node,3)= (dN(node,1)*Cof31 + dN(node,2)*Cof32 + dN(node,3)*Cof33) / DetJ;
              end

            % Stiffness elementary matrix construction. Symmetrical storage.
              for no = 1:NodesNumber
                  for nd = no:NodesNumber

                      m = (no-1)*3;
                      j = (nd-1)*3;

                      Na1b1= D(no,1)*D(nd,1);
                      Na1b2= D(no,1)*D(nd,2);
                      Na1b3= D(no,1)*D(nd,3);
                      Na2b1= D(no,2)*D(nd,1);
                      Na2b2= D(no,2)*D(nd,2);
                      Na2b3= D(no,2)*D(nd,3);
                      Na3b1= D(no,3)*D(nd,1);
                      Na3b2= D(no,3)*D(nd,2);
                      Na3b3= D(no,3)*D(nd,3);

                      K(m+1,j+1) = K(m+1,j+1) + W*DetJ*(C1*Na1b1+C3*(Na2b2+Na3b3));
                      K(m+2,j+2) = K(m+2,j+2) + W*DetJ*(C1*Na2b2+C3*(Na1b1+Na3b3));
                      K(m+3,j+3) = K(m+3,j+3) + W*DetJ*(C1*Na3b3+C3*(Na1b1+Na2b2));
                      K(m+1,j+2) = K(m+1,j+2) + W*DetJ*(C2*Na1b2+C3*Na2b1);
                      K(m+1,j+3) = K(m+1,j+3) + W*DetJ*(C2*Na1b3+C3*Na3b1);
                      K(m+2,j+3) = K(m+2,j+3) + W*DetJ*(C2*Na2b3+C3*Na3b2);

                      if (nd > no)
                          K(m+2,j+1) = K(m+2,j+1)+W*DetJ*(C2*Na2b1+C3*Na1b2);
                          K(m+3,j+1) = K(m+3,j+1)+W*DetJ*(C2*Na3b1+C3*Na1b3);
                          K(m+3,j+2) = K(m+3,j+2)+W*DetJ*(C2*Na3b2+C3*Na2b3);
                      end
                  end
              end

            % Mass elementary matrix construction. Symmetrical storage.
              for no = 1:NodesNumber
                  for nd = no:NodesNumber
                      m = (no-1)*3;
                      j = (nd-1)*3;
                      M(m+1,j+1) = M(m+1,j+1) + W*DetJ*rho*N(no)*N(nd);
                      M(m+2,j+2) = M(m+2,j+2) + W*DetJ*rho*N(no)*N(nd);
                      M(m+3,j+3) = M(m+3,j+3) + W*DetJ*rho*N(no)*N(nd);
                  end
              end
          end % End of the integration loop in the 3 directions...
      end
  end    
end % end of function
