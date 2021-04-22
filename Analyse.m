%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2; fs = 16; Nx = 64;

data = load("output_mesh.out");
x = data(1,:);
y = data(2,:);

data = load("output_E.out");
t = data(:,1);

data = load("output_f.out");
f = data(:,2:end);

size(x)
size(y)
size(f)
h = figure;
for i = 1:size(t)
    f = data(1+(i-1)*Nx:i*Nx,2:end);
%     plot(y,f);
%     xlim([0 6])
%     ylim([-3 3])
%     
%     plot(x,f);
%     xlim([0 10])
%     ylim([-3 3])
%
%     plot3(x,y,f);
%     xlim([0 10])
%     ylim([0 6])
%     zlim([-3 3])
%
    surf(y,x,f);
    ylim([0 10])
    xlim([0 6])
    zlim([-3 3])
    xlabel('Y', 'fontsize', fs)
    ylabel('X', 'fontsize', fs)
    zlabel('Z', 'fontsize', fs)

%     view(90,0)
%     view(0, 90);
     drawnow
%      pause(1.5)
%     pause(0.005)

  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 

  if i == 1 
      imwrite(imind,cm,'2D_fixe.gif','gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,'2D_fixe.gif','gif','WriteMode','append'); 
  end 

end

% for i = 1:size(t)
%     f = data(1+(i-1)*Nx:i*Nx,2:end);
%     surf(x,y,f);
% %     xlim([0 10])
% %     ylim([0 6])
% %     zlim([-3 3])
% %     legend('X','Y','Z')
% %     view(90,0)
% %     view(0, 90);
%     drawnow
% %     pause(0.005)
% end
% data = load("tolo.out");
% f = data(:,2:end);
% surf(f);