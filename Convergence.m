%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2; fs = 16;
u = 4; nx = 64; ny = 64; Lx = 10; Ly = 6; m = 2; n = 2; hx = Lx/(nx-1); hy = Ly/(ny-1);
tfin = 1.286239;

repertoire = ''; % Chemin d'accès au code compilé
executable = 'Exercice7'; % Nom de l'exécutable
input = 'configuration.in'; % Nom du fichier d'entrée

nsimul = 3;

nper = round(logspace(2,3,nsimul)); % Nombre d'itérations entier de 10^2 à 10^4  
CFL = (63./nper);
dt = tfin./nper;

paramstr = 'CFL'; % Nom du paramètre à scanner  
param = CFL; % Valeurs du paramètre à scanner 

%% Simulations %%
%%%%%%%%%%%%%%%%%
output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    para = sprintf('%s%s%s', paramstr, '=', num2str(param(i)));
    cmd = sprintf('%s %s %s %s%s %s CFL=%.15g output_file=%s', 'set', 'path=%path:C:\Program Files\MATLAB\R2020b\bin\win64;=%', '&', repertoire, executable, input, CFL(i), output{i});
    system(cmd);
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%%%%%
data = load("output_mesh.out");
x = data(1,:);
y = data(2,:);

fana = zeros(nx,ny);
for j = 1:nx
    for l = 1:ny
        fana(j,l) = cos(pi*((x(j)/5) + (y(l)/3)));
    end
end

error = zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i});
    fnum = data(end-nx+1:end,2:end);

    result = 0;
    for j = 1:nx-1
        for l = 1:ny-1
            result = result + (fana(j,l) - fnum(j,l))^2 + (fana(j+1,l) - fnum(j+1,l))^2 + (fana(j,l+1) - fnum(j,l+1))^2 + (fana(j+1,l+1) - fnum(j+1,l+1))^2;
        end
    end
    
    error(i) = 0.5*sqrt(hx*hy*result);
end

%% Figures
p = polyfit(CFL,error, 1)

figure
loglog(CFL,error,'k+')
set(gca,'fontsize',fs)
xlabel('$\beta_{CFL}$', 'interpreter', 'latex','fontsize', fs)
ylabel('Error')
s = sprintf('y = (%.2f)x  + (%.2f)', p(1) ,p(2));
text(9, 1.5, s, 'Color','red', 'FontSize', 16)
grid on

figure
loglog(dt,error,'k+')
set(gca,'fontsize',fs)
xlabel('$\Delta t$', 'interpreter', 'latex','fontsize', fs)
ylabel('Error')
s = sprintf('y = (%.2f)x  + (%.2f)', p(1) ,p(2));
text(9, 1.5, s, 'Color','red', 'FontSize', 16)
grid on