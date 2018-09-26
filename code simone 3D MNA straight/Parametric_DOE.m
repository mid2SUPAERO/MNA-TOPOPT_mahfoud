T=20:2.5:80;
x=2500:100:4500;
[xx,TT]=meshgrid(x,T);
xx=xx';
TT=TT';
perfo=zeros(size(xx));
fan_disp=zeros(size(T));
mu=zeros(size(xx));
for l=1:length(x)
    for k=1:length(T)
        [perfo(l,k),fan_disp(l,k),mu(l,k)]=perfo_func(x(l),T(k));
    end
end
surf(xx,TT,perfo)
xlabel('x [mm]')
ylabel( '\tetha [°]')
zlabel('\Delta TSFC %')
figure (2)
surf(xx,TT,fan_disp)
xlabel('x [mm]')
ylabel( '\theta [°]')
zlabel('FAN center axial displacement [mm]%')
