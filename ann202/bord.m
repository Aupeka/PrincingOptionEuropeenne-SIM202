function A = bord(I,J, Refneu)
% Vérifie si une face est sur la frontière 
  if ((Refneu(I) ==1)&&(Refneu(J) ==1))
    A = 1;
  else
    A=0;
  end ;
