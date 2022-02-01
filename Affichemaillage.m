
figure;
hold on

% maillage
trimesh(Numtri,Coorneu(:,1),Coorneu(:,2),'Color','k'); 
view(2);
axis('equal');

% ajouter eventuellement un titre
title('Maillage');

hold off;