function val = f(x,y)
% Fonction second membre dans l'équation de Poisson
cas =1 ;
if cas == 1
    val = 2*pi^2*sin(pi*x)*sin(pi*y) ; 
end ;
if cas == 2
    val = 0 ;
end ; 


