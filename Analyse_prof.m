%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=64; % Doit etre consistant avec l'input du code C++
fichier = 'output';
fid=fopen([fichier,'_mesh.out']); % ouverture du fichier
x=str2num(fgetl(fid)); % lit la 1e ligne du fichier, qui contient les x_i
y=str2num(fgetl(fid)); % lit la 2e ligne du fichier, qui contient les y_j
u = load([fichier,'_u.out']); % lit tout le fichier des u(x_i,t_j)
data = load([fichier,'_E.out']); % fichier de (t,E)
t = data(:,1);
E = data(:,2);
data = load([fichier,'_f.out']);
% Le fichier contient, ligne par ligne, {t_k, {f(i,j),j=1,Nx}, i=i,Ny},
% k=1,nsteps
[ii,jj]=size(data);
Ny=jj-1;
nsteps=ii/Nx;

% Exemple 1: extraire le tableau des f(i,j) au temps final
istart=(nsteps-1)*Nx +1; % indice de la premiere ligne a lire du tableau global (data)
iend=istart+Nx-1; % indice de la derniere ligne a lire du tableau global (data)\
f=(data(istart:iend,2:jj))'; % extrait le f(i,j); le ' transpose pour avoir f(x_i,y_j)
fs=16; % font size
% Contour plot de f(x,y)
figure
contourf(x,y,f) % lignes de niveaux avec coloriage (le ' transpose)
% voir aussi contour, surf, surfc, ...
set(gca,'fontsize',fs) % set font size for labels, title, etc
xlabel('x [m]')
ylabel('y [m]')
colorbar

% % Exemple 2: une coupe a y=y0=const de f(x,y,t) --&gt; fcut(x_i,t_k)
% % extraire le tableau des f(i,j,k) : f(x,y=y0,t) a y0 fixe.
% i0=2; % indice pour y_0 fixe (ici on prend le 2e point de maillage en y)
% fcut=(data(i0:Ny:end,2:jj))'; % toutes les Ny lignes, on est au meme y
% t=data(i0:Ny:end,1); % la 1e colonne contient le temps
% 
% % contour plot de f(x,y0,t)
% 
% figure
% contourf(x,t,fcut')
% set(gca,'fontsize',fs) % set font size for labels, title, etc
% xlabel('x [m]')
% ylabel('t [s]')
% colorbar
