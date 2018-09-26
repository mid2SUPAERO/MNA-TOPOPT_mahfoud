%% ********************* polynomials test script *********************** %%
clear ; close all ;
t = 0.5 ;
A1 = [  1   t^2   t^4    t^6     t^8     t^12 ;...
        0   2*t  4*t^3  6*t^5   8*t^7    12*t^11; ...
        0    2  4*3*t^2 6*5*t^4 8*7*t^6  12*11*t^10] ;
A2 = [1 1 1 1 1 1 ; 0 2 4 6 8 12 ; 0 2 4*3  6*5 8*7 12*11] ;
A = [A1;A2] ;
B = [1-6*t^8+6*t^12;-6*8*t^7+6*12*t^11;-6*8*7*t^6+6*12*11*t^10;0;0;0] ;
V = A\B ;
%%
x = linspace(0,1,100) ;
f = zeros(1,100) ;
g = zeros(1,100) ;
for i = 1:100
    if x(i)<0.5
        f(i) = 1-6*x(i)^2+6*x(i)^3 ;
    else
        f(i) = 2-6*x(i)+6*x(i)^2-2*x(i)^3 ;
    end
    if x(i)<0.5
        g(i) = 1-6*x(i)^8+6*x(i)^12 ;
    else
        g(i) = V(1)+V(2)*x(i)^2+V(3)*x(i)^4+V(4)*x(i)^6+V(5)*x(i)^8+V(6)*x(i)^12 ;
    end
end
figure;
plot(x,g)
hold on 
plot(x,f)
hold off
legend('new function','old function')
xlabel('x') ; ylabel('new & old functions') ;