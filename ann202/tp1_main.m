
%%
% -------------------------------- %
%        Lecture du maillage       %
% -------------------------------- %

clear all
global cas ;
cas = 2 ;
if cas == 1 
    nom_maillage = 'Cas1.msh' ;
end ;
if cas == 2 
    nom_maillage = 'Cas2.msh' ;
end ;


[Nbpt, Nbtri, Coorneu, Refneu, Numtri, ...
 Reftri, Nbaretes, Numaretes, Refaretes] = lecture_msh(nom_maillage);

%%
% ---------------------------------------------------------- %
%     Définition des données utiles pour notre problème      %
% ---------------------------------------------------------- %

    % Matrice des numéros des triangles des faces intérieures. Les dernières lignes (nulles) seront ignorées ensuite
Reffacetri = zeros(3*Nbtri, 2);     

    % Matrice des sommets des faces intérieures. Contient un premier élément pour les besoins du code, qui sera ignoré ensuite
Numface = [0 0];

    % Matrice des sommets des faces du bord, Nbaretes etant le nb de faces sur le bord
Reffacetri_bord = zeros(Nbaretes, 1); 

    %Nb de faces intérieures (Nb_faces_int)
s = 1;

    %Remplissage --> Boucle sur les triangles
for l=1:Nbtri
  % Pour les faces IJ et JK
    for i=1:2
        I = Numtri(l,i);
        J = Numtri(l,i+1);
        inNumface = is_in([min(I,J),max(I,J)], Numface);
        if (bord(I,J, Refneu) == 0 )            %Si F n'est pas sur le bord
          if (inNumface==0)                     %Si F n'est pas déjà dans la liste
            s=s+1;
            Numface = [Numface ; [min(I,J),max(I,J)]];
            Reffacetri(s,1) = l;
          else
            Reffacetri(inNumface,2) = l;
          end ;         
        end ;
    end
    
    % Pour la face KI
    I = Numtri(l,3);
    J = Numtri(l,1);
        inNumface = is_in([min(I,J),max(I,J)], Numface);
        if bord(I,J, Refneu) == 0             %Si F n'est pas sur le bord
          if (inNumface==0)                   %Si F n'est pas déjà dans la liste
            s=s+1;
            Numface = [Numface ; [min(I,J),max(I,J)]];
            Reffacetri(s,1) = l;
          else
            Reffacetri(inNumface,2) = l;
          end ;
        end ;
    
        
    for k = 1: Nbaretes %Aretes sur le bord
      I = min(Numaretes(k,1),Numaretes(k,2));
      J = max(Numaretes(k,1),Numaretes(k,2));
      if (length(intersect([I,J], Numtri(l,:)))==2)
          Reffacetri_bord(k,1) = l;
      end ;
    end ;
end % for l


%%
% ---------------------------------------------------------- %
%   Définition des matrices de masse et du second membre     %
% ---------------------------------------------------------- %


Nb_faces_int = length(Numface);
XF = [] ; % Matrice des coordonnées des milieux des faces


    %Définition de KK
KK = sparse(Nb_faces_int-1,Nb_faces_int-1);

    %Definition de BB
BB = [];

% /!\ Le premier élément de Numface n'est pas un vrai élément
for i = 2:Nb_faces_int %Face F
    
        %Id face "i"
    S = Numface(i,:); 
    S11 = Coorneu(S(1),:);  % Coordonnées du premier sommet de la face 
    S12 = Coorneu(S(2),:);  % Coordonnées du second sommet de la face
    xF1 = milieu(S11,S12);  % Milieu de la face
    nF1 = normale(S11,S12); % Normale à la face
          
    for j = 2:Nb_faces_int %Face F'
      
        %Les deux faces appartiennent-elles au même triangle ?
      triangle_commun = intersect(Reffacetri(i,:), Reffacetri(j,:)); % Numéro dans Numtri de l'éventuel triangle commun 
      if (isempty(triangle_commun) ~=1) % Si les deux faces sont sur le même triangle
      
          %Id face "j"
        S = Numface(j,:);
        S21 = Coorneu(S(1),:);
        S22 = Coorneu(S(2),:);
        nF2 = normale(S21,S22);
        xF2 = milieu(S21,S22);
            
          if (i ~= j)%Les deux faces sont pas les mêmes
            
              % Vérification de l'orientation de la normale --> Sinon réorientation vers l'extérieur du triangle
           if (nF1'*vec(xF2,xF1) <= 0) 
            nF1 = -nF1;               
          end ;
          if (nF2'*vec(xF1,xF2) <= 0)
           nF2 = -nF2;
          end ;
          
          %Calcul de KK
          
            T = Numtri(triangle_commun,:);    
            T1 = Coorneu(T(1),:);
            T2 = Coorneu(T(2),:);
            T3 = Coorneu(T(3),:);

            % Calcul du coefficient associé dans K
            KK(i-1,j-1) = ((normef([S11;S12])*normef([S21;S22]))/normetri([T1,T2,T3]))*(nF1'*nF2);
          
          else %i=j
            % Sinon, on ne considère qu'une seule face
        x = xF1(1);
        y = xF1(2);
        XF = [ XF ; x y ] ;
        
            % Calcul de l'aire des deux triangles adjacents à F
        tri = Reffacetri(i,:) ; 
        tri_1 = Numtri(tri(1), :) ;
        tri_2 = Numtri(tri(2), :) ;
        P1 = Coorneu(tri_1(1),:);
        P2 = Coorneu(tri_1(2),:);
        P3 = Coorneu(tri_1(3),:);
        M1 = Coorneu(tri_2(1),:);
        M2 = Coorneu(tri_2(2),:);
        M3 = Coorneu(tri_2(3),:);
        K_1 = normetri([P1,P2,P3]);
        K_2 = normetri([M1,M2,M3]);
        
            % Calcul de la diagonale de K
        KK(i-1,j-1) = normef([S11;S12])^2*(1/K_1 + 1/K_2) ;
              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calcul du second membre %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %On le calcul ici car la formule trouvée demande d'avoir des
        %informations sur les deux triangles annexes à la face qu'on
        %considère
        
        BB = [BB; (K_1 + K_2)*f(x,y,cas)/3];
        
          end;
      end;
              
    end; %j
end; %i

if cas==2
Surface = sparse(Nb_faces_int-1, Nbaretes) ;  % Matrice de rigidité surfacique 
  bb = [] ; % Fonction du bord aux milieux des faces du bord
  for j = 1: Nbaretes
  S = Numaretes(j,:); 
  S11 = Coorneu(S(1),:);  % Coordonnées du premier sommet de la face 
  S12 = Coorneu(S(2),:);  % Coordonnées du second sommet de la face
  xF1 = milieu(S11,S12);  % Milieu de la face
  nF1 = normale(S11,S12); % Normale à la face
  for i = 2:Nb_faces_int
    % Même procédure pour la deuxième face
    S = Numface(i,:);
    S21 = Coorneu(S(1),:);
    S22 = Coorneu(S(2),:);
    nF2 = normale(S21,S22);
    xF2 = milieu(S21,S22);
      if (nF1'*vec(xF2,xF1) <= 0) % Vérification de l'orientation de la normale
        nF1 = -nF1;               % Sinon réorientation vers l'extérieur du triangle
      end ;
      if (nF2'*vec(xF1,xF2) <= 0)
       nF2 = -nF2;
      end ;
    
      tricom = intersect(Reffacetri(i,:), Reffacetri_bord(j,:)); % Numéro dans Numtri de l'éventuel triangle commun 
      if (isempty(tricom) ~=1) % Si les deux faces sont sur le même triangle
     
        M = Numtri(tricom,:);    
        P1 = Coorneu(M(1),:);
        P2 = Coorneu(M(2),:);
        P3 = Coorneu(M(3),:);
        
        % Calcul du coefficient associé dans G
        Surface(i-1,j) = ((normef([S11;S12])*normef([S21;S22]))/normetri([P1,P2,P3]))*(nF1'*nF2);
      end ;  
    end ; 
    % Calcul de g(xF1)
    bb = [ bb ; frontiere(xF1(1), xF1(2)) ] ;
  end ;
    
  % Calcul de G
  G = Surface*bb ;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Resolution  du Problème        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cas ==1
U_h = KK\BB;
end;
if cas==2
    U_h= KK\(-G);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Calcul de l'erreur           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
Err_2 = 0; % epsilon = grad(u-uh)
GRAD_UH=sparse(Nbtri,2);
GRAD_U_UH=sparse(Nbtri,2);
Uaff = sparse(Nbtri, 3);
% Calcul de (grad(u-uh))^2 
for t = 1: Nbtri
  % Intégrale L2 sur tous les triangles
  Faces = [ ] ;
  Faces_bord = [ ] ;
  
  % Recherche des faces du triangle
  % Recherche des faces internes au domaine
  i = 2;
  while (length(Faces)<3 && i<= length(Reffacetri))
    if (any(Reffacetri(i,:) == t))
        Faces = [ Faces ; i] ; % Relevé du numéro de la face interne trouvée
    end ;
    i = i+1 ;
  end ;
  
  % Si le triangle a 2 faces internes, alors il a 1 face au bord 
  % Recherche de cette face
  if (length(Faces) == 2) 
    i = 1 ;
    while (length(Faces_bord)<1 && i<= length(Reffacetri_bord))
        if (Reffacetri_bord(i) == t)
            Faces_bord = [ Faces_bord ; i] ;
        end ;
        i = i+1 ;
    end ;
  end ;

  % Si le triangle a 1 face interne, alors il a 2 faces au bord 
  % Recherche de ces faces
  if (length(Faces) == 1)
  i = 1 ;
  while (length(Faces_bord)<2 && i<= length(Reffacetri_bord))
    if (Reffacetri_bord(i) == t)
      Faces_bord = [ Faces_bord ; i] ;
    end ;
    i = i+1 ;
   end ;
  end ;
  
  % Si le triangle est interne au domaine : 3 sommets/3 valeurs connus
  % -1 car numéros décalés, voir déclaration de Numface
  if (length(Faces) == 3)
    u1 = U_h(Faces(1)-1) ;
    x1 = XF(Faces(1)-1,1) ;
    y1 = XF(Faces(1)-1,2) ;
  
    u2 = U_h(Faces(2)-1) ;
    x2 = XF(Faces(2)-1,1) ;
    y2 = XF(Faces(2)-1,2) ;    
 
    u3 = U_h(Faces(3)-1) ;
    x3 = XF(Faces(3)-1,1) ;
    y3 = XF(Faces(3)-1,2) ;
  end ; 
  
  % Si le triangle a 1 face au bord : 2 sommets internes/2 valeurs connus
  % -1 car numéros décalés, voir déclaration de Numface
  if (length(Faces) == 2)
    u1 = U_h(Faces(1)-1) ;
    x1 = XF(Faces(1)-1,1) ;
    y1 = XF(Faces(1)-1,2) ;
  
    u2 = U_h(Faces(2)-1) ;
    x2 = XF(Faces(2)-1,1) ;
    y2 = XF(Faces(2)-1,2) ;    
     
    % Valeur au bord et coordonnées du centre de la face externe     
    arete = Faces_bord(1) ;
    Coor = milieu(Coorneu(Numaretes(arete,1),:), Coorneu(Numaretes(arete,2),:)) ;
    x3 = Coor(1) ;
    y3 = Coor(2) ;
    
    if cas == 1 
        u3 = 0 ;
    end ; 
    if cas == 2 
        u3 = frontiere(x3,y3) ;
    end ;
    
  end ;   
  
  % Si le triangle a 2 faces au bord : 1 sommet interne/1 valeur connus
  % -1 car numéros décalés, voir déclaration de Numface
  if (length(Faces) == 1)
    u1 = U_h(Faces(1)-1) ;
    x1 = XF(Faces(1)-1,1) ;
    y1 = XF(Faces(1)-1,2) ;
  
    % Valeur au bord et coordonnées du centre de la 1ere face externe
    arete_1 = Faces_bord(1) ;
    Coor = milieu(Coorneu(Numaretes(arete_1,1),:), Coorneu(Numaretes(arete_1,2),:)) ;
    x2 = Coor(1) ;
    y2 = Coor(2) ;
    if cas == 1 
        u2 = 0 ;
    end ; 
    if cas == 2 
        u2 = frontiere(x2,y2) ;
    end ;
    
    % Valeur nulle au bord, calcul des coordonnées du centre de la 2eme face
    arete_2 = Faces_bord(2) ;
    Coor = milieu(Coorneu(Numaretes(arete_2,1),:), Coorneu(Numaretes(arete_2,2),:)) ;
    x3 = Coor(1) ;
    y3 = Coor(2) ;
    if cas == 1 
        u3 = 0 ;
    end ; 
    if cas == 2 
        u3 = frontiere(x3,y3) ;
    end ;
  end ; 
  
  % Calcul du gradient de uh affine, connaissant 3 coordonnées et 3 valeurs
  a = ((u1-u2)*(y1-y3) - (u1-u3)*(y1-y2)) / ((x1-x2)*(y1-y3) - (x1-x3)*(y1-y2)) ;
  b = ((u1-u2)*(x1-x3) - (u1-u3)*(x1-x2)) / ((y1-y2)*(x1-x3) - (y1-y3)*(x1-x2)) ;
  c = u1-a*x1-b*y1;

  grad_uh = [ a ; b ] ;
  GRAD_UH(t,:)=grad_uh;
  Uaff(t,:) = [a , b , c];
  % Approximation de l'intégrale L2 par les valeurs aux 3 points-milieux
  grad_1 = gradu(x1,y1,cas)-grad_uh ;
  norm_1 = grad_1'*grad_1 ;
  grad_2 = gradu(x2,y2,cas)-grad_uh ;
  norm_2 = grad_2'*grad_2 ;
  grad_3 = gradu(x3,y3,cas)-grad_uh ;
  norm_3 = grad_3'*grad_3 ;
  
  % Pour le calcul de l'aire du triangle
  S1 = Coorneu(Numtri(t,1),:) ;
  S2 = Coorneu(Numtri(t,2),:) ;
  S3 = Coorneu(Numtri(t,3),:) ;
 
  % Ajout à la somme sur les triangles
  Err_2 = Err_2 + (normetri([S1,S2,S3])/3)*(norm_1 + norm_2 + norm_3) ;
  GRAD_U_UH(t)=sqrt((normetri([S1,S2,S3])/3)*(norm_1 + norm_2 + norm_3));
end ;

% Calcul de grad(u-uh)
Err = sqrt(Err_2) ;

[Err ; (Nb_faces_int -1)]

%Problème d'affichage des chiffres significatifs
% --> format long g

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Affichage Erreur           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Valeurs calculées pour 5 maillages différents
Err_4 = [1.1790 ; 0.5835 ; 0.2910 ; 0.145422223962109];
Err_4cas2=[ 0.1877 ; 0.1387 ; 0.0898 ];
eta_4 = [1.4774 ; 0.7253 ; 0.3581 ; 0.1783 ];
Ieff_4 = [1.2531 ; 1.2431 ; 1.2304 ; 1.2263 ];
h_4 = [ 1/4; 1/8; 1/16];
V_4cas2=[ 76 ; 272 ; 1120 ];
V_4 = [20; 88; 368; 1504];

    %Affichage
figure(1);
loglog(V_4, Err_4,'red');
hold on
loglog(V_4, eta_4,'blue');

legend('Erreur','eta');
figure(2);
semilogx(V_4,Ieff_4);

figure(8);
loglog(h_4,Err_4cas2)
title('Erreur en fonction de h');

figure(9);
loglog(V_4cas2,Err_4cas2)
title('Erreur en fonction de V')
    %Calcul de la pente
pente_Vcas2=(log(0.0898)-log(0.1877))/log(1120/76)
pente_hcas2=(log(0.1877)-log(0.0898))/log(8)
pente_V = (log(0.145422223962109) - log(1.1790))/ log(1504/20);
pente_h = (log(1.1790) - log(0.145422223962109))/log(8);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Affichage Solution          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure ;
uh = plot3(XF(:,1), XF(:,2), U_h, 'o') ;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Reconstruction potentiel     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_h = zeros(Nbpt, 1) ;
 
for i = 1:Nbpt
    if Refneu(i)==1
        if cas ==2
            s_h(i)=frontiere(Coorneu(i,1),Coorneu(i,2));
        else
            s_h(i)=0;
        end
    else
        L = in_tri(i,Numtri, Nbtri);
        T = length(L);
        S = Coorneu(i,:);
        S = [S, 1];
        for l = 1:T
            s_h(i) = s_h(i) + Uaff(L(l),:)*S';                   
        end ;
        s_h(i) = s_h(i)/T;
    end ;
end;


figure(15);
affiche(s_h,Numtri,Coorneu,"Potentiel");


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Reconstruction flux équilibré     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dans cette partie on a choisi de ne calculer que le vecteur f moyen sur
%tous les triangles sachant que calculer sigma_h ne sert pas dans la suite
%vu que la partie gradient de u_h s'annule dans l'erreur à posteriori.


    %Calcul de f_h
    
   f_h=sparse(Nbtri,1); % Vecteur des valeurs moyennes de f /d sur toutes les maille
   sigma_h=sparse(Nbtri,1);
   XK=zeros(Nbtri,2);
   for i=1:Nbtri
      x1= Coorneu(Numtri(i,1),1);
      y1= Coorneu(Numtri(i,1),2);
      x2= Coorneu(Numtri(i,2),1);
      y2= Coorneu(Numtri(i,2),2);
      x3= Coorneu(Numtri(i,3),1);
      y3= Coorneu(Numtri(i,3),2);
      f_h(i)=f((x1+x2)/2,(y1+y2)/2,cas)+f((x3+x2)/2,(y3+y2)/2,cas)+f((x1+x3)/2,(y1+y3)/2,cas);
      f_h(i)=f_h(i)/6; %on a divisé par 6 car d=2
      XK(i,1)= (x1+x2+x3)/3;
      XK(i,2)= (y1+y2+y3)/3;
   end

   

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Erreur à posteriori      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Calcul de l'erreur grad uh-sh

ss_h=zeros(Nb_faces_int-1,1);
for i=2:Nb_faces_int
    ss_h(i-1)=(s_h(Numface(i,1))+s_h(Numface(i,2)))/2;
end

%grad de sh

eta_s_h_2 = zeros(Nbtri,1);
GRAD_SH=sparse(Nbtri,2);
% Calcul de (grad(u-uh))^2 
for t = 1: Nbtri
  % Intégrale L2 sur tous les triangles
  Faces = [ ] ;
  Faces_bord = [ ] ;
  
  % Recherche des faces du triangle
  % Recherche des faces internes au domaine
  i = 2;
  while (length(Faces)<3 && i<= length(Reffacetri))
    if (any(Reffacetri(i,:) == t))
        Faces = [ Faces ; i] ; % Relevé du numéro de la face interne trouvée
    end ;
    i = i+1 ;
  end ;
  
  % Si le triangle a 2 faces internes, alors il a 1 face au bord 
  % Recherche de cette face
  if (length(Faces) == 2) 
    i = 1 ;
    while (length(Faces_bord)<1 && i<= length(Reffacetri_bord))
        if (Reffacetri_bord(i) == t)
            Faces_bord = [ Faces_bord ; i] ;
        end ;
        i = i+1 ;
    end ;
  end ;

  % Si le triangle a 1 face interne, alors il a 2 faces au bord 
  % Recherche de ces faces
  if (length(Faces) == 1)
  i = 1 ;
  while (length(Faces_bord)<2 && i<= length(Reffacetri_bord))
    if (Reffacetri_bord(i) == t)
      Faces_bord = [ Faces_bord ; i] ;
    end ;
    i = i+1 ;
   end ;
  end ;
  
  % Si le triangle est interne au domaine : 3 sommets/3 valeurs connus
  % -1 car numéros décalés, voir déclaration de Numface
  if (length(Faces) == 3)
    u1 = ss_h(Faces(1)-1) ;
    x1 = XF(Faces(1)-1,1) ;
    y1 = XF(Faces(1)-1,2) ;
  
    u2 = ss_h(Faces(2)-1) ;
    x2 = XF(Faces(2)-1,1) ;
    y2 = XF(Faces(2)-1,2) ;    
 
    u3 = ss_h(Faces(3)-1) ;
    x3 = XF(Faces(3)-1,1) ;
    y3 = XF(Faces(3)-1,2) ;
  end ; 
  
  % Si le triangle a 1 face au bord : 2 sommets internes/2 valeurs connus
  % -1 car numéros décalés, voir déclaration de Numface
  if (length(Faces) == 2)
    u1 = ss_h(Faces(1)-1) ;
    x1 = XF(Faces(1)-1,1) ;
    y1 = XF(Faces(1)-1,2) ;
  
    u2 = ss_h(Faces(2)-1) ;
    x2 = XF(Faces(2)-1,1) ;
    y2 = XF(Faces(2)-1,2) ;    
     
    % Valeur au bord et coordonnées du centre de la face externe     
    arete = Faces_bord(1) ;
    Coor = milieu(Coorneu(Numaretes(arete,1),:), Coorneu(Numaretes(arete,2),:)) ;
    x3 = Coor(1) ;
    y3 = Coor(2) ;
    
    if cas == 1 
        u3 = 0 ;
    end ; 
    if cas == 2 
        u3 = (s_h(Numaretes(arete,1))+s_h(Numaretes(arete,2)))/2 ;
    end ;
    
  end ;   
  
  % Si le triangle a 2 faces au bord : 1 sommet interne/1 valeur connus
  % -1 car numéros décalés, voir déclaration de Numface
  if (length(Faces) == 1)
    u1 = ss_h(Faces(1)-1) ;
    x1 = XF(Faces(1)-1,1) ;
    y1 = XF(Faces(1)-1,2) ;
  
    % Valeur au bord et coordonnées du centre de la 1ere face externe
    arete_1 = Faces_bord(1) ;
    Coor = milieu(Coorneu(Numaretes(arete_1,1),:), Coorneu(Numaretes(arete_1,2),:)) ;
    x2 = Coor(1) ;
    y2 = Coor(2) ;
    if cas == 1 
        u2 = 0 ;
    end ; 
    if cas == 2 
        u2 = (s_h(Numaretes(arete_1,1))+s_h(Numaretes(arete_1,2)))/2 ;
    end ;
    
    % Valeur nulle au bord, calcul des coordonnées du centre de la 2eme face
    arete_2 = Faces_bord(2) ;
    Coor = milieu(Coorneu(Numaretes(arete_2,1),:), Coorneu(Numaretes(arete_2,2),:)) ;
    x3 = Coor(1) ;
    y3 = Coor(2) ;
    if cas == 1 
        u3 = 0 ;
    end ; 
    if cas == 2 
        u3 = (s_h(Numaretes(arete_2,1))+s_h(Numaretes(arete_2,2)))/2 ;
    end ;
  end ; 
  
  % Calcul du gradient de sh affine, connaissant 3 coordonnées et 3 valeurs
  a = ((u1-u2)*(y1-y3) - (u1-u3)*(y1-y2)) / ((x1-x2)*(y1-y3) - (x1-x3)*(y1-y2)) ;
  b = ((u1-u2)*(x1-x3) - (u1-u3)*(x1-x2)) / ((y1-y2)*(x1-x3) - (y1-y3)*(x1-x2)) ;
  c = u1-a*x1-b*y1;

  grad_sh = [ a ; b ] ;
  GRAD_SH(t,:)=grad_sh;
end ;


for t=1:Nbtri
    S1 = Coorneu(Numtri(t,1),:) ;
    S2 = Coorneu(Numtri(t,2),:) ;
    S3 = Coorneu(Numtri(t,3),:) ;
    eta_s_h_2(t)=normetri([S1,S2,S3])*((GRAD_UH(t,1)-GRAD_SH(t,1))^2 + (GRAD_UH(t,2)-GRAD_SH(t,2))^2);
end

% Calcul de l'erreur grad uh +sigma_h
eta_sigma_2=sparse(Nbtri,1);
 if cas==1
     for i=1:Nbtri
         S1 = Coorneu(Numtri(i,1),:) ;
         S2 = Coorneu(Numtri(i,2),:) ;
         S3 = Coorneu(Numtri(i,3),:) ;
         x_12=(S1+S2)/2;
         x_13=(S1+S3)/2;
         x_23=(S2+S3)/2;
         S_k=(S1+S2+S3)/3;
         eta_sigma_2(i)=f_h(i)*f_h(i)*(normetri([S1,S2,S3])/3)*((x_12-S_k)*(x_12-S_k)'+(x_13-S_k)*(x_13-S_k)'+(x_23-S_k)*(x_23-S_k)');
     end
 end  
 
eta_2=zeros(Nbtri,1);
for i=1:Nbtri
    eta_2(i)=eta_sigma_2(i)+eta_s_h_2(i);
end
eta1=sqrt(sum(eta_2));
Ieff=eta1/Err;

eta=sqrt(eta_2);

for i=1:Nbtri
    a=Coorneu(Numtri(i,1),:);
    b=Coorneu(Numtri(i,2),:);
    c=Coorneu(Numtri(i,3),:);
    X = [a(1),b(1),c(1)];
    Y = [a(2),b(2),c(2)];
    fill(X,Y,eta(i))
    hold on
end



    
    
    
