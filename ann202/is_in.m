function A = is_in(F, Numface)
% Vérifie si la face F est déjà répertoriée
  S1 = F(1);
  S2 = F(2);
  B = 0;
for k= 1:length(Numface(:,1))
  x1 = Numface(k,1);
  x2 = Numface(k,2);
  if (x1 == S1)
    if (x2== S2)
      B = k;
    end ;
  end ;
end ;
if (B~=0)
  A = B;
else
  A = 0;
end ;
  
    
