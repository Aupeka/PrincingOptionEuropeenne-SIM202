function A = bord(I,J, Refneu)
% V�rifie si une face est sur la fronti�re 
  if ((Refneu(I) ==1)&&(Refneu(J) ==1))
    A = 1;
  else
    A=0;
  end ;
