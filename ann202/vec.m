function V = vec(S1, S2)
% Vecteur d'un sommet à un autre
  x = S2(1) -S1(1);
  y = S2(2) -S1(2);
  V = [x;y];