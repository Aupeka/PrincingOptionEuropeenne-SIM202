
figure;
hold on

% maillage
trimesh(Numtri,Coorneu(:,1),Coorneu(:,2)); 
view(2);
axis('equal');

% ajouter eventuellement un titre
title('Maillage');

hold off;