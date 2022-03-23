function L = in_tri(i,Numtri, Nbtri)
  L= [];
  for l = 1:Nbtri
    K = Numtri(l,:);
    if (intersect(i,K)==i)
      L = [L;l];
    end ;
  end ;

  
  
