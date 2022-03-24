function vec = gradu(x,y,cas)
% Calcul du gradient de u 
if cas == 1 
    vec = pi*[ cos(pi*x)*sin(pi*y) ; cos(pi*y)*sin(pi*x) ] ; 
end ;
if cas == 2
    [t,r] = cart2pol(x,y) ;
    pol = (2/3)*r^(-1/3)*[sin(2*t/3) ; cos(2*t/3)] ;
    vec = [ cos(t)*pol(1) - sin(t)*pol(2) ; sin(t)*pol(1) + cos(t)*pol(2) ] ;
end ;
