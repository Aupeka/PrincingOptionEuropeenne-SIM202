function A = normef(V)
  S1 = V(1,:);
  S2 = V(2,:);
  U = vec(S1,S2);
  A = norm(U) ;
