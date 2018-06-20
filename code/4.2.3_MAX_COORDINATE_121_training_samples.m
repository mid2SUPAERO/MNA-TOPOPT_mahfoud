%all data points
x = [
    200 200;200 400;200 600;200 800;200 1000;200 -200;200 -400;200 -600;200 -800;200 -1000;200 0;
    400 200;400 400;400 600;400 800;400 1000;400 -200;400 -400;400 -600;400 -800;400 -1000;400 0;
    600 200;600 400;600 600;600 800;600 1000;600 -200;600 -400;600 -600;600 -800;600 -1000;600 0;
    800 200;800 400;800 600;800 800;800 1000;800 -200;800 -400;800 -600;800 -800;800 -1000;800 0;
    1000 200;1000 400;1000 600;1000 800;1000 1000;1000 -200;1000 -400;1000 -600;1000 -800;1000 -1000;1000 0;
    -200 200;-200 400;-200 600;-200 800;-200 1000;-200 -200;-200 -400;-200 -600;-200 -800;-200 -1000;-200 0;
    -400 200;-400 400;-400 600;-400 800;-400 1000;-400 -200;-400 -400;-400 -600;-400 -800;-400 -1000;-400 0;
    -600 200;-600 400;-600 600;-600 800;-600 1000;-600 -200;-600 -400;-600 -600;-600 -800;-600 -1000;-600 0;
    -800 200;-800 400;-800 600;-800 800;-800 1000;-800 -200;-800 -400;-800 -600;-800 -800;-800 -1000;-800 0;
    -1000 200;-1000 400;-1000 600;-1000 800;-1000 1000;-1000 -200;-1000 -400;-1000 -600;-1000 -800;-1000 -1000;-1000 0;
    0 200;0 400;0 600;0 800;0 1000;0 -200;0 -400;0 -600;0 -800;0 -1000;];

%all response
y = [
    8.19;13.45;18.80;23.28;30.55;8.136;13.51;18.95;25.46;30.64;6.809;
    13.51;16.380;21.50;26.90;32.49;13.45;16.27;21.57;27.01;32.55;13.62;
    18.87;21.67;24.57;29.55;34.83;18.89;21.51;24.41;29.62;34.96;20.43;
    25.38;27.02;29.90;32.76;37.64;25.39;26.89;29.61;32.64;37.80;27.23;
    30.70;32.6285;35.02;38.10;40.49;30.56;32.40;34.81;37.71;40.68;34.04;
    8.217;13.36;18.77;25.21;30.32;7.968;13.37;18.70;25.13;30.34;6.746;
    13.59;16.43;21.35;26.72;32.18;13.45;15.94;21.28;26.73;32.32;13.49;
    19.06;21.73;24.65;29.34;34.59;18.78;21.52;23.90;29.15;34.56;20.24;
    25.62;27.18;29.91;32.87;37.34;25.24;26.89;29.61;31.87;36.99;26.98;
    30.86;32.75;35.19;38.16;41.08;30.51;32.46;34.83;37.71;39.84;33.73;
    6.877;13.75;20.63;27.51;34.38;6.815;13.63;20.44;27.26;34.07;];
% position of max stress 
y_posi = [
    -25.305 7.98809;-12.7048 -1.2833;-12.2325 0.0855;-12.131 0.55796;-12.0583 1.03571;-13.8805 7.15288;-12.6543 5.1686;-12.3361 4.29443;-12.1215 3.3892;-12.0541 2.92874;-12 2;
   -23.2833 9.2952;-25.305 7.98809;-13.1537 -2.1387;-12.7048 -1.2833;-12.3622 -0.38;-23.1686 -5.34573;-13.8716 7.1423;-13.0718 6;-12.6543 5.16864;-12.3361 4.29443;-12 2;
    -21.9145 9.7675;-24.1387 8.84625;-25.305 7.98809;-13.4161 -2.54452;-12.9164 -1.71779;-22.2944 -5.66392;-24 -4.9282;-13.8716 7.1423;-13.3161 6.39607;-12.8509 5.5904;-12 2;
    -21.442 9.86896;-23.2833 9.2952;-24.5445 8.5838;-25.305 7.98809;-13.4161 -2.54452;-21.3892 -5.87846;-23.1686 -5.34573;-24.3961 -4.6839;-13.8716 7.1423;-13.3161 6.39607;-12 2;
    -20.9643 9.94167;-22.38 9.63777;-23.7178 9.08365;-24.5445 8.5838;-25.305 7.98809;-20.9287 -5.94591;-22.2944 -5.66392;-23.5904 -5.149;-24.3961 -4.6839;-13.8716 7.1423;-12 2;
    -14.8577 8.12836;-27.3457 -1.16864;-27.6639 -0.294426;-27.8785 0.610815;-27.9459 1.07126;-14.695 -3.98809;-27.2952 5.2833;-27.7675 3.91453;-27.869 3.44204;-27.9417 2.9643;-28 2;
    -16.8314 9.3457;-14.8577 8.12836;-26.9282 -2;-27.3457 -1.16864;-27.6639 -0.294426;-16.7167 -5.2952;-14.695 -3.98809;-26.8463 6.1387;-27.2952 5.2833;-27.6378 4.38;-28 2;
    -17.7056 9.66392;-16 8.9282;-14.8577 8.12836;-26.6839 -2.396;-27.1491 -1.5904;-18.0855 -5.76753;-15.8613 -4.8462;-14.695 -3.98809;-26.8463 6.1387;-27.0836 5.71779;-28 2;
    -18.6108 9.87846;-16.8314 9.3457;-15.6039 8.6839;-14.8577 8.12836;-26.6839 -2.396;-18.558 -5.86896;-16.7167 -5.2952;-15.4555 -4.58387;-14.695 -3.98809;-26.5839 6.54452;-28 2;
    -19.0713 9.94591;-17.7056 9.6639;-16.4096 9.1491;-15.6039 8.6839;-14.8577 8.12836;-19.0357 -5.94167;-17.62 -5.63777;-16.2822 -5.08365;-15.4555 -4.58387;-14.695 -3.98809;-28 2;
    -20 10;-20 10;-20 10;-20 10;-20 10;-20 -6;-20 -6;-20 -6;-20 -6;-20 -6;];
% the center point is (-20,2) so make a translation to (0,0) to do training
% and moving back afterwards
y_posi_trans= zeros(length(y_posi),2);
for i=1:length(y_posi)
    y_posi_trans(i,1)=y_posi(i,1)-(-20);
    y_posi_trans(i,2)=y_posi(i,2)-(2);
end
% to make sure that all points in the cycle
% try_list=zeros(length(y_posi),1)
% for i=1:length(y_posi)
%     try_list(i,1)=sqrt(y_posi_trans(i,1)^2+y_posi_trans(i,2)^2);
% end
% try_list
% create the angle list to save the angle of max stress point coordinate
angle_training= zeros(length(y_posi),1);
for i=1:length(y_posi)
    if(y_posi_trans(i,1)>=0)
       angle_training(i)=atan(y_posi_trans(i,2)/y_posi_trans(i,1));
    elseif (y_posi_trans(i,1)<0 && y_posi_trans(i,2)>=0)
       angle_training(i)=atan(y_posi_trans(i,2)/y_posi_trans(i,1))+pi;
        else
            angle_training(i)=atan(y_posi_trans(i,2)/y_posi_trans(i,1))-pi;
    end
end

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

% the center point is (-20,2) so make a translation to (0,0) to do training
% and moving back afterwards
y_test_posi_trans= zeros(length(y_test_posi),2);
for i=1:length(y_test_posi)
    y_test_posi_trans(i,1)=y_test_posi(i,1)-(-20);
    y_test_posi_trans(i,2)=y_test_posi(i,2)-(2);
end


% create the test angle list to save the angle of max stress point coordinate
angle_test= zeros(length(y_test_posi),1);
for i=1:length(y_test_posi)
    if(y_test_posi_trans(i,1)>=0)
       angle_test(i)=atan(y_test_posi_trans(i,2)/y_test_posi_trans(i,1));
    elseif (y_test_posi_trans(i,1)<0 && y_test_posi_trans(i,2)>=0)
       angle_test(i)=atan(y_test_posi_trans(i,2)/y_test_posi_trans(i,1))+pi;
        else
            angle_test(i)=atan(y_test_posi_trans(i,2)/y_test_posi_trans(i,1))-pi;
    end
end


% the x y coordinate of max stress(training)
y_posi_x_trans = y_posi_trans(:,1);
y_posi_y_trans = y_posi_trans(:,2);
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
% 2. prediction of angle
% meanfunc=  {@meanSum, {@meanConst,{@meanPoly,3}}}; hyp.mean = [0;0.000001;0.000001;0.00001;0.0001;0.00001;0.0001;]
% meanfunc=  {@meanSum, {@meanConst,{@meanWSPC,1},{@meanPoly,2}}}; hyp.mean = [0;0;0;0;0;0.00001;0.0001;0.0001;0.0001];
% meanfunc=  {@meanSum, {@meanConst,{@meanPoly,3}}}; hyp.mean = [0;0.00001;0.00001;0.00001;0.00001;0.00001;0.00001;]
% meanfunc=  {@meanSum, {@meanConst,{@meanPoly,4}}}; hyp.mean = [0;0.0000001;0.0000001;0.0000001;0.0000001;0.0000001;0.0000001;0.0000001;0.000001;]
% meanfunc=  {@meanSum, {@meanConst, {@meanLinear}}}; hyp.mean = [0;5.36;5.36];
% meanfunc=  {@meanSum, {@meanLinear, {@meanPoly,2}}}; hyp.mean = [0.3;0;0.0001;0.0001;0.0001;0.0001];
% meanfunc=   @meanLinear; hyp.mean = [0; 0];
  meanfunc=[];hyp.mean = [];
%  meanfunc=   @meanConst; hyp.mean = [0;0];
%  covfunc = @covSEiso; hyp.cov = [5.2; 5.2];
% covfunc=  {@covSum, {@covSEiso,{@covMaterniso, 3}}}; hyp.cov = [5; 0.5;1; -0.5];
covfunc=  {@covSum, {@covSEiso,{@covRQiso},{@covMaterniso, 3}}}; hyp.cov = [10;10;0.1;0.1;0.01;log(0.25);log(1);];
% covfunc=  {@covSum,{@covRQiso, {@covMaterniso, 3}}}; hyp.cov = [0;0.1;0.1;0.1;0.1];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

hyp6 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x,angle_training);
%prediction
[ymu_angle ys_angle fmu_angle fs_angle] = gp(hyp6, @infGaussLik, meanfunc, covfunc,likfunc,x,angle_training,x_test);

% the RMSE between ptrfiction angle and real angle
RMSE_angle = sqrt(mean((angle_test- ymu_angle(:,1)).^2))
% plot the the RMSE
p_angle=polyfit(angle_test,ymu_angle,1)

% x1=-3:1:2;  
% y1=polyval(p_angle,x1);  
% 
% figure('NumberTitle', 'off', 'Name', 'Comparision of angle coordinate between GP Prediction and Real Value');
% plot(ymu_angle,angle_test,'*',x1,y1);
% title('Comparision between GP Prediction and Real Value');
% xlabel('Max Stress (GP Prediction)');
% ylabel('Max Stress (Real Value)');
% legend('','fitting line');


% recover the xy coordiante from angle
ymu_coordinate= zeros(length(ymu_angle),2);
for i=1:length(ymu_angle)
    ymu_coordinate(i,1)=-20+8*cos(ymu_angle(i)) ;
    ymu_coordinate(i,2)=2+8*sin(ymu_angle(i)) ;
  
end
prediction_list(:,3)=ymu_coordinate(:,1);
prediction_list(:,4)=ymu_coordinate(:,2);

% the RMSE between ptrfiction x coordibate and real value
RMSE_x = sqrt(mean(( y_test_posi_x- ymu_coordinate(:,1)).^2))
RMSE_y = sqrt(mean(( y_test_posi_y- ymu_coordinate(:,2)).^2))

% plot the the RMSE
p_x=polyfit(ymu_coordinate(:,1),y_test_posi_x,1)
x1=-30:1:-12;  
y1=polyval(p_x,x1);  
figure('NumberTitle', 'off', 'Name', 'Comparision of x coordinate between GP Prediction and Real Value');
plot(ymu_coordinate(:,1),y_test_posi_x,'*',x1,y1);
title('Comparision between GP Prediction and Real Value');
xlabel('Max Stress (GP Prediction)');
ylabel('Max Stress (Real Value)');
legend('','fitting line');


% plot the the RMSE
p_y=polyfit(ymu_coordinate(:,2),y_test_posi_y,1)
x1=-15:1:10;  
y1=polyval(p_y,x1);  
figure('NumberTitle', 'off', 'Name', 'Comparision of y coordinate between GP Prediction and Real Value');
plot(ymu_coordinate(:,2),y_test_posi_y,'*',x1,y1);
title('Comparision between GP Prediction and Real Value');
xlabel('Max Stress (GP Prediction)');
ylabel('Max Stress (Real Value)');
legend('','fitting line');


% new points to give a prediction
%  xs1 = randi([-1000 1000],7,1);
%  xs2 = randi([-1000 1000],7,1);
%  xs = [xs1 xs2]
 
 x_test2= [630    94;
   812   915;
  -746   930;
   827  -685;
   265   942;
  -805   915;
  -443   -29;];
 x_test2_coordinate= [-12    2;
   -13.7025   -2.93375;
  -26.6839   -2.39607;
   -24.7773  -4.41699;
   -12.131   0.55796;
   -26.417   -2.77727;
  -28   2;];

y_posi_trans1= zeros(length(x_test2_coordinate),2);
for i=1:length(x_test2_coordinate)
    y_posi_trans1(i,1)=x_test2_coordinate(i,1)-(-20);
    y_posi_trans1(i,2)=x_test2_coordinate(i,2)-(2);
end

% try_list=zeros(length(y_posi_trans1),1);
% for i=1:length(y_posi_trans1)
%     try_list(i,1)=sqrt(y_posi_trans1(i,1)^2+y_posi_trans1(i,2)^2);
% end
% try_list