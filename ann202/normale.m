function V = normale(S1, S2)
% Normale à la face reliant les sommets S1 et S2
  U = vec(S1, S2);
  N = norm(U) ;
  V = (1/N)*[-U(2); U(1)];