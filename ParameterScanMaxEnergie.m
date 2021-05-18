%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = 'cmake-build-debug\'; % Chemin d'acces au code compile
executable = 'Exo7.exe'; % Nom de l'executable
input = 'cmake-build-debug\configuration.in'; % Nom du fichier d'entree de base

nsimul = 100; % Nombre de simulations ¨


%32 utilisé pour le mesh

pert_velocity = linspace(4,6,nsimul);

paramstr = 'pert_velocity';    % Nom du parametre a scanner 
param = pert_velocity;         % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

% output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
% for i = 1:nsimul
%     output{i} = ['cmake-build-debug\SimulationOut\',paramstr, '=', num2str(param(i),16), '.out']
%     if ~isfile(strcat (repertoire,executable))
%         'erreur' 
%     end
%     cmd = sprintf('%s%s %s %s=%.15g output_energy=%s', repertoire, executable, input, paramstr, param(i), output{i})
%     disp(cmd)
%     system(cmd);
% end
% 
% 
% %% Analyse 7.2 e)
% 
% emax = zeros(1,nsimul);
% for j = 1:nsimul
%     data = load(output{j});
%     
%     t = data(:,1);
%     E = data(:,2);  
%     emax(j) = max(E);
%    
% end

figure
plot(pert_velocity,emax,'.', 'linewidth', 2.5, 'markersize', 9);

hold on
plot([4.88492831206 4.88492831206],[0 140],'--','Color', [1 0 0], 'linewidth', 2)

grid on

legend('simulation','$\omega_{2,2} \approx 4.88$','interpreter','latex')


xlabel('$\omega[s^{-1}]$','interpreter', 'latex','fontsize',2)
ylabel('$\hat E(\omega)[m^{3}]$','interpreter', 'latex','fontsize',2)

ax=gca;
ax.FontName = 'Times';
ax.FontSize = 15;







