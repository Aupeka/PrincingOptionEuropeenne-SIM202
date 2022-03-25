function val = f(x,y,cas)
% Fonction second membre dans l'équation de Poisson
if cas == 1
    val = 2*pi^2*sin(pi*x)*sin(pi*y) ; 
end ;
if cas == 2
    val = 0 ;
end ; 


