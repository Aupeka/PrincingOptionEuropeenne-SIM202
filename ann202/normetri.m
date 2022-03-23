function A = normetri(V)
% Aire d'un triangle
  M = [ V(3) - V(1) , V(5) - V(1) ; V(4) - V(2) V(6)-V(2) ] ;
  A = 0.5*abs(det(M)) ; 

