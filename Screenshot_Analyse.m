%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2; fs = 16; Nx = 10;

i = 6000; %Timestamp à analyser

repertoire = ''; % Chemin d'acces au code compile
executable = 'Exercice7'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree

%% Simulations et analyse %%
%%%%%%%%%%%%%%%%%

% cmd = sprintf('%s %s %s %s%s %s','set', 'path=%path:C:\Program Files\MATLAB\R2020b\bin\win64;=%', '&', repertoire, executable, input);
% system(cmd);
% disp('Done.')

data = load("output_mesh.out");
x = data(1,:);
y = data(2,:);

data = load("output_E.out");
t = data(:,1);

data = load("output_f.out");
f = data(:,2:end);

size(t)

%Dessin de la fonction
f = data(1+(i-1)*Nx:i*Nx,2:end);
surf(x,y,f);
xlim([0 5000])
ylim([0 2000])
zlim([-3 3])
xlabel('X', 'fontsize', fs)
ylabel('Y', 'fontsize', fs)
zlabel('Z', 'fontsize', fs)