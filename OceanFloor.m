% figure fond ocean cas 1 et 2 piur le  7.3)

h0=1000;
h1=20;
a=2000;
b=5000;


x1 = linspace(0,2000,30);
x=linspace(2000,5000,40);
y=linspace(0,2000,70);

[X,Y] = meshgrid(x,y);
[X1,Y1] = meshgrid(x1,y);


Z1=-1000*X1./X1;
Z=-(h0-(h0-h1)*sin((3.14159265358979323846*(X-a))/(b-a)).*sin((3.14159265358979323846*Y)/2000));



figure
surf(X,Y,Z)

hold on

surf(X1,Y1,Z1)



xlabel('X[m]', 'fontsize', fs)
ylabel('Y[m]', 'fontsize', fs)
zlabel('Z[m]', 'fontsize', fs)
colorbar
close

%figure en 2D pour le cas 1

z= -(h0-(h0-h1)*sin((3.14159265358979323846*(x-a))/(b-a)));
figure
plot(x,z,'-', 'linewidth', 2.5);


grid on


xlabel('$x[m]$','interpreter', 'latex','fontsize',2)
ylabel('$z[m]$','interpreter', 'latex','fontsize',2)

ax=gca;
ax.FontName = 'Times';
ax.FontSize = 15;





