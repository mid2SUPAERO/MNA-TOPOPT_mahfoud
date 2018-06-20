close all;

%training data points
x = [333 0;666 0;1000 0; -333 0;-666 0;-1000 0;0 0;
    333 333;666 333;1000 333; -333 333;-666 333;-1000 333; 0 333;
    333 666;666 666;1000 666; -333 666;-666 666;-1000 666; 0 666;
    333 1000;666 1000;1000 1000; -333 1000;-666 1000;-1000 1000; 0 1000;
    333 -333;666 -333;1000 -333; -333 -333;-666 -333;-1000 -333; 0 -333;
    333 -666;666 -666;1000 -666; -333 -666;-666 -666;-1000 -666; 0 -666;
    333 -1000;666 -1000;1000 -1000; -333 -1000;-666 -1000;-1000 -1000; 0 -1000;];
%max principle stress for training data
y = [11.34;22.67;34.04; 11.23;22.46;33.73;0;
    13.63;22.50;31.44; 13.68;22.63;31.76; 11.45;
    22.40;27.27;36.10; 22.25;27.36;36.22; 22.90;
    31.32;35.83;40.94; 31.27;35.58;41.08; 34.38;
    13.55;22.39;31.47; 13.27;22.39;31.29; 11.35;
    22.49;27.09;35.84; 22.25;26.53;35.87; 22.69;
    31.58;35.94;40.68; 31.16;35.45;39.84; 34.07;];
 
 x_test=[
  -700  -400
   300   900
   400   500
   600  -400
  -300   200
  -500   600
  -300   700
  
   700  -500
   500   200
   600  -400
     0   600
   100  -500
   300  -1000
  -800  -400
  ];

y_realvalue=[
  24.13
  28.19
  18.82
  21.51
  10.87
  22.56
  22.65
 
  25.58
  16.31
  21.51
  20.63
  15.32
  31.14
  26.89
  ];
% begin a time mark to compute the training time
 tic
%prior
% meanfunc=  {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1];
% meanfunc=   @meanConst; hyp.mean = 1;
% meanfunc= {@meanPoly,2};hyp.mean =[1;1;1;1];
% meanfunc=  {@meanSum, {@meanConst, {@meanPoly,2}}}; hyp.mean = [10; 0.001;0.001;0.001;0.001];
% meanfunc=   @meanLinear; hyp.mean = [0.5; 1];
meanfunc=[]; hyp.mean=[];
%  covfunc = @covSEiso; hyp.cov = [9; 9];
%  covfunc = {@covMaterniso, 3}; ell = 50; sf = 50; hyp.cov = log([ell; sf]);
%  covfunc = @covRQiso; hyp.cov = [1; 1;0.01];%
%   covfunc=  {@covSum,{@covRQiso, {@covMaterniso, 3}}}; hyp.cov = [15;1;1;5;3];
%  covfunc=  {@covSum, {@covSEiso,{@covMaterniso, 3}}}; hyp.cov = [15;1;5;5];
%  covfunc=  {@covSum, {@covSEiso,{@covRQiso}}}; hyp.cov = [5;5;5;10;10];
covfunc=  {@covSum, {@covSEiso,{@covRQiso},{@covMaterniso, 3}}}; hyp.cov = [0;30;5;10;1;5;1;];

likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

%optimize the hyperparameters
hyp5 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x,y);
%prediction
[ymu5 ys55 fmu5 fs55] = gp(hyp5, @infGaussLik, meanfunc, covfunc,likfunc,x,y,x_test);
% end time mark
toc

% plot
figure('NumberTitle', 'off', 'Name', 'training,testing and evaluating');
plot3(x(:,1),x(:,2),y,'k.','MarkerSize',20)
title('training,testing and evaluating');
hold on;
plot3(x_test(:,1),x_test(:,2),fmu5,'r.','MarkerSize',20)
hold on;
plot3(x_test(:,1),x_test(:,2),y_realvalue,'b.','MarkerSize',20)
xlabel('Fx');
ylabel('Fy');
zlabel('Max Principle Stress');
legend('training data','prediction value','real value','Location','NorthEast');

% the RMSE between ptrfiction and real value
RMSE = sqrt(mean(( y_realvalue- ymu5).^2))
% the R square error 
R_Square = 1- (sum(( y_realvalue- ymu5).^2))/(sum(( y_realvalue- mean(y_realvalue)).^2))

% plot the the RMSE
p2=polyfit(ymu5,y_realvalue,1)
x1=1:1:40;  
y1=polyval(p2,x1);  

figure('NumberTitle', 'off', 'Name', 'Comparision between GP Prediction and Real Value');
plot(ymu5,y_realvalue,'*',x1,y1);
title('Comparision between GP Prediction and Real Value');
xlabel('Max Stress (GP Prediction)');
ylabel('Max Stress (Real Value)');
legend('','fitting line');