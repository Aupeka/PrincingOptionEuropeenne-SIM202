function A = frontiere(x,y)
% Evalue la valeur de g en un point donné de la frontière
if x==0 && y==0
    A=0;
else
r = sqrt(x^2 + y^2) ;
if x>=0
    t = atan(y/x) ;
    A = r^(2/3)*sin(2*t/3) ;
else
    t= pi+atan(y/x);
    A = r^(2/3)*sin(2*t/3) ;
end
end
end
