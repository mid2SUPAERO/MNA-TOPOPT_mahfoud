%all data points
x = [333 0;666 0;1000 0; -333 0;-666 0;-1000 0;
    333 333;666 333;1000 333; -333 333;-666 333;-1000 333; 0 333;
    333 666;666 666;1000 666; -333 666;-666 666;-1000 666; 0 666;
    333 1000;666 1000;1000 1000; -333 1000;-666 1000;-1000 1000; 0 1000;
    333 -333;666 -333;1000 -333; -333 -333;-666 -333;-1000 -333; 0 -333;
    333 -666;666 -666;1000 -666; -333 -666;-666 -666;-1000 -666; 0 -666;
    333 -1000;666 -1000;1000 -1000; -333 -1000;-666 -1000;-1000 -1000; 0 -1000;];
%all response
y = [11.34;22.67;34.04; 11.23;22.46;33.73;
    13.63;22.50;31.44; 13.68;22.63;31.76; 11.45;
    22.40;27.27;36.10; 22.25;27.36;36.22; 22.90;
    31.32;35.83;40.94; 31.27;35.58;41.08; 34.38;
    13.55;22.39;31.47; 13.27;22.39;31.29; 11.35;
    22.49;27.09;35.84; 22.25;26.53;35.87; 22.69;
    31.58;35.94;40.68; 31.16;35.45;39.84; 34.07;];
% position of max stress 
y_posi = [-12 2;-12 2;-12 2; -28 2;-28 2;-28 2;
    -25.305 7.98809;-23.2833 9.2952;-21.9145 9.76753; -14.8577 8.12836;-16.8314 9.34573;-17.7056 9.66392; -20 10;
    -12.695 -1.26147;-25.305 7.98809;-24.1387 8.84625; -27.3457 -1.16864;-14.8577 8.12836;-16 8.9282; -20 10;
    -12.2167 0.15059;-13.1537 -2.1387;-25.305 7.98809; -27.6533 -0.32972;-26.9282 -2;-14.8577 8.12836; -20 10;
    -13.8716 7.1423;-23.1686 -5.3457;-22.2944 -5.66392;-14.695 -3.98809;-16.7167 -5.2952;-18.0855 -5.76753; -20 -6;
    -12.6543 5.16864;-13.8716 7.1423;-24 -4.9282; -27.2952 5.2833;-14.695 -3.98809;-15.8613 -4.84625; -20 -6 ;
    -12.3361 4.29443;-13.0718 6;-13.8716 7.1423; -27.7675 3.91453;-26.8463 6.1387;-14.695 -3.98809; -20 -6;];
% testing data
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
% max stress of testing
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
% position of testing
 y_test_posi=[
  -16.2822  -5.08365
   -12.2325 0.08547
   -13.4161 -2.54452
   -24      -4.9282
   -16      8.9282
   -26.417  -2.77727
   -27.5175 -0.736161
  
   -24.3961 -4.6839
   -22.38   9.63777
   -24      -4.9282
   -20      10
   -12.0541 2.92874
   -12.2156 3.84493
   -16.7167 -5.2952];

% the x y coordinate of max stress(training)
y_posi_x = y_posi(:,1);
y_posi_y = y_posi(:,2);
% the x y coordinate of max stress(testing)
y_test_posi_x = y_test_posi(:,1);
y_test_posi_y = y_test_posi(:,2);

% the final output list include the input, max stress and its position
prediction_list=zeros(14,5);
prediction_list(:,1:2)=x_test;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. prediction to max stress
meanfunc=  {@meanSum, {@meanConst, {@meanPoly,2}}}; hyp.mean = [10; 0.001;0.001;0.001;0.001];
%  meanfunc=   @meanLinear; hyp.mean = [0.5; 1];
covfunc = @covSEiso; hyp.cov = [0; 0];
% covfunc=  {@covProd, {@covSEiso,{@covMaterniso, 3}}}; hyp.cov = [0; 0;log(0.25); log(1)];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

hyp5 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x,y);
%prediction
[ymu_max ys_max fmu_max fs_max] = gp(hyp5, @infGaussLik, meanfunc, covfunc,likfunc,x,y,x_test);
prediction_list(:,3)=ymu_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. prediction of x coordinate
% meanfunc=  {@meanSum, {@meanConst, {@meanPoly,2}}}; hyp.mean = [-20;0.004;0.004;0.004;0.004];
% meanfunc=  {@meanSum, {@meanConst, {@meanLinear}}}; hyp.mean = [-22;1;1;];
% meanfunc=  {@meanSum, {@meanLinear, {@meanPoly,2}}}; hyp.mean = [0.1;0.1;0.004;0.004;0.004;0.004];
%  meanfunc=   @meanLinear; hyp.mean = [0.5; 1];
meanfunc = []; hyp.mean=[];
covfunc=  {@covSum, {@covSEiso,{@covRQiso},{@covMaterniso, 3}}}; hyp.cov = [3;3;0.01;0.01;0.01;log(0.25);log(1);];
% covfunc = @covSEiso; hyp.cov = [0; 0];
% covfunc=  {@covProd, {@covSEiso,{@covMaterniso, 3}}}; hyp.cov = [0; 0;log(0.25); log(1)];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

hyp6 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x,y_posi_x);
%prediction
[ymu_x ys_x fmu_x fs_x] = gp(hyp6, @infGaussLik, meanfunc, covfunc,likfunc,x,y_posi_x,x_test);
prediction_list(:,4)=ymu_x;
% the RMSE between ptrfiction x coordibate and real value
RMSE_x = sqrt(mean(( y_test_posi_x- ymu_x).^2))

% plot the the RMSE
p_x=polyfit(ymu_x,y_test_posi_x,1)

x1=-21:1:-17;  
y1=polyval(p_x,x1);  

figure('NumberTitle', 'off', 'Name', 'Comparision of x coordinate between GP Prediction and Real Value');
plot(ymu_x,y_test_posi_x,'*',x1,y1);
title('Comparision between GP Prediction and Real Value');
xlabel('Max Stress (GP Prediction)');
ylabel('Max Stress (Real Value)');
legend('','fitting line');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. prediction of y coordinate
% meanfunc=  {@meanSum, {@meanConst, {@meanPoly,2}}}; hyp.mean = [1; 0.001;0.001;0.001;0.001];
meanfunc=[];hyp.mean=[];
%  meanfunc=   @meanLinear; hyp.mean = [0.5; 1];
covfunc=  {@covSum, {@covSEiso,{@covRQiso},{@covMaterniso, 3}}}; hyp.cov = [3;3;0.01;0.01;0.01;log(0.25);log(1);];
% covfunc=  {@covProd, {@covSEiso,{@covMaterniso, 3}}}; hyp.cov = [0; 0;log(0.25); log(1)];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

hyp7 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x,y_posi_y);
%prediction
[ymu_y ys_y fmu_y fs_y] = gp(hyp7, @infGaussLik, meanfunc, covfunc,likfunc,x,y_posi_y,x_test);
prediction_list(:,5)=ymu_y;
% the RMSE between ptrfiction x coordibate and real value
RMSE_y = sqrt(mean(( y_test_posi_y- ymu_y).^2))

% plot the the RMSE
p_y=polyfit(ymu_y,y_test_posi_y,1)

x1=-10:1:10;  
y1=polyval(p_y,x1);  

figure('NumberTitle', 'off', 'Name', 'Comparision of y coordinate between GP Prediction and Real Value');
plot(ymu_y,y_test_posi_y,'*',x1,y1);
title('Comparision between GP Prediction and Real Value');
xlabel('Max Stress (GP Prediction)');
ylabel('Max Stress (Real Value)');
legend('','fitting line');

prediction_list;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%