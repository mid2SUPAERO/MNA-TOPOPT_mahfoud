close all;

global MAXI_STRESS_COORD;
%all data points
x = [333 0;666 0;1000 0; -333 0;-666 0;-1000 0;0 0;
    333 333;666 333;1000 333; -333 333;-666 333;-1000 333; 0 333;
    333 666;666 666;1000 666; -333 666;-666 666;-1000 666; 0 666;
    333 1000;666 1000;1000 1000; -333 1000;-666 1000;-1000 1000; 0 1000;
    333 -333;666 -333;1000 -333; -333 -333;-666 -333;-1000 -333; 0 -333;
    333 -666;666 -666;1000 -666; -333 -666;-666 -666;-1000 -666; 0 -666;
    333 -1000;666 -1000;1000 -1000; -333 -1000;-666 -1000;-1000 -1000; 0 -1000;];
%all response
y = [11.34;22.67;34.04; 11.23;22.46;33.73;0;
    13.63;22.50;31.44; 13.68;22.63;31.76; 11.45;
    22.40;27.27;36.10; 22.25;27.36;36.22; 22.90;
    31.32;35.83;40.94; 31.27;35.58;41.08; 34.38;
    13.55;22.39;31.47; 13.27;22.39;31.29; 11.35;
    22.49;27.09;35.84; 22.25;26.53;35.87; 22.69;
    31.58;35.94;40.68; 31.16;35.45;39.84; 34.07;];

%we want to pick out 7 in the x, which means we pick one in every 7.
test_list_subscript=zeros(7,1);
test_list_input=zeros(7,2);
test_list_real_value=zeros(7,1);

for n=1:7
    test_list_subscript(n)=randi([1 7],1)+7*(n-1);
    test_list_input(n,:)=x(test_list_subscript(n),:);
    test_list_real_value(n)=y(test_list_subscript(n));
end

% delete the test data and leave the train data
x(test_list_subscript,:)=[];
y(test_list_subscript,:)=[];

train_list_input=x;
train_list_output=y;

%prior
meanfunc=  {@meanSum, {@meanConst, {@meanPoly,2}}}; hyp.mean = [10; 0.001;0.001;0.001;0.001];
% meanfunc=  {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1];
% meanfunc= {@meanPoly,1};hyp.mean =[1;1];
covfunc = @covSEiso; hyp.cov = [0; 0];
% covfunc=  {@covProd, {@covSEiso,{@covMaterniso, 3}}}; hyp.cov = [0; 0;log(0.25); log(1)];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

%optimize the hyperparameters
hyp4 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, train_list_input,train_list_output);
%prediction
[ymu4 ys44 fmu4 fs44] = gp(hyp4, @infGaussLik, meanfunc, covfunc,likfunc,train_list_input,train_list_output,test_list_input);

% plot
figure('NumberTitle', 'off', 'Name', 'training,testing and evaluating');
plot3(x(:,1),x(:,2),y,'k.','MarkerSize',20)
title('training,testing and evaluating');
hold on;
plot3(test_list_input(:,1),test_list_input(:,2),fmu4,'r.','MarkerSize',20)
hold on;
plot3(test_list_input(:,1),test_list_input(:,2),test_list_real_value,'b.','MarkerSize',20)
xlabel('Fx');
ylabel('Fy');
zlabel('Max Principle Stress');
legend('training data','prediction value','real value','Location','NorthEast');
% the RMSE between ptrfiction and real value
RMSE = sqrt(mean(( test_list_real_value- fmu4).^2))

% plot the the RMSE
p=polyfit(ymu4,test_list_real_value,1)

x1=5:1:45;  
y1=polyval(p,x1);  

figure('NumberTitle', 'off', 'Name', 'Comparision between GP Prediction and Real Value');
plot(ymu4,test_list_real_value,'*',x1,y1);
title('Comparision between GP Prediction and Real Value');
xlabel('Max Stress (GP Prediction)');
ylabel('Max Stress (Real Value)');
legend('','fitting line');

