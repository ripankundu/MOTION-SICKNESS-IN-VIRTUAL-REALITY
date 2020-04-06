%Data initialization%
data_camera = csvread('CAMERA_1.csv');
camera_time = data_camera(:,1);
camera_data_w = data_camera(:,2);
hold on
data_sensor = csvread('SENSOR_1.csv');
sensor_time = data_sensor(:,1);
sensor_data_w = data_sensor(:,2);

%Data processing&
p=2   % data length quadratic %
l=3 %data length cubic%
q=1  %data length Linear%
predict=[]% Array initialization of Cubic Data%
predicterror=[]% Array initialization of  predicted Error for Cubic Data%
predictlength=[]% Array initialization of predicted length for Cubic Data%
actuallength=[]% Array initialization of the actual length%

predictq=[] % Array initialization of Quadratic  Data%
predicterrorq=[]% Array initialization of predicted Error for Quadratic data Data%
predictlengthq=[]% Array initialization of predicted length of Quadratic Data%


predictl=[] % Array initialization of Linear   Data%
predicterrorl=[]% Array initialization of predicted Error for Linear data Data%
predictlengthl=[]% Array initialization of predicted length of Linear Data%
%noprediction=[]%Array initialization of no prediction%


d=[]; % Array initialization of difference between cubic and quadratic data%
n=4; % Data position%
i=1;
t=[]; % Array resizing for plot operation for camera data value%
k=[]; % Array resizing for plot operation for camera time value%
u=[];% Array resizing for plot operation for camera data value%
v=[];% Array resizing for plot operation for camera time value%

while n<311
    %camera Data Quadratic%
 x = camera_time(n-p:n,:); 
 y = camera_data_w(n-p:n,:);
 
%  Camera Data cubic% 
 xc = camera_time(n-l:n,:); 
 yc = camera_data_w(n-l:n,:);
%  
 %camera Data Linear %
 xl = camera_time(n-q:n,:); 
 yl = camera_data_w(n-q:n,:);
 
 u=camera_time(l:n+1,:);
 v=camera_data_w(l:n+1,:);
 
 %Cubic curve fit
parms = [xc.^3 xc.^2 xc ones(size(xc))]\yc;
xv = linspace(min(xc), max(xc))';
y_fit = [xv.^3 xv.^2 xv ones(size(xv))]*parms;


% %quadratic curve fit%
 parms1 = [x.^2 x ones(size(x))]\y;
 y_fit1 = [xv.^2 xv ones(size(xv))]*parms1;
 
%Linear curve fit%

  parms2 = [xl ones(size(xl))]\yl;
 y_fit2 = [xv ones(size(xv))]*parms2;
 


%Cubic parameters%
 cof1=parms(1);
 cof2=parms(2);
 cof3=parms(3);
 cof4=parms(4);
 
 %Quadratic parameter%
 coefq1=parms1(1);
 coefq2=parms1(2);
 coefq3=parms1(3);
 
 %Linear Parameters%
 coefl1=parms2(1);
 coefl2=parms2(2);


% Predicted Cubic data%
 predict(i)=cof1*(camera_time(n+1))^3+ cof2*(camera_time(n+1))^2+ cof3*camera_time(n+1)+ cof4;
 
%Predicted quadratic data%
 predictq(i)=coefq1*(camera_time(n+1))^2+ coefq2*(camera_time(n+1))+ coefq3;
 
%  
 %Predicted Linear data%
 predictl(i)=coefl1*(camera_time(n+1))+coefl2;
 
% No prediction Value%
noprediction(i)=abs(camera_data_w(i+1)-camera_data_w(i));
 
 
 %Difference between cubic and quadratic data%
 d(i)=abs(predict(i)-predictq(i));
 
 %Predicted error for cubic data%
 predicterror(i)=abs(camera_data_w(n+1)-predict(i));
 
 
 %predicted error for Quadratic data%
 predicterrorq(i)=abs(camera_data_w(n+1)-predictq(i));
 
 %predicted error for Linear data%
 predicterrorl(i)=abs(camera_data_w(n+1)-predictl(i));
 
 
 %Predicted length for cubic data%
 predictlength(i)=abs(predict(i)-camera_data_w(n));
 
 %Predicted length for quadratic data%
 predictlengthq(i)=abs(predictq(i)-camera_data_w(n));
 
 
 %Predicted length for quadratic data%
 predictlengthl(i)=abs(predictl(i)-camera_data_w(n));
 
 
 
 
 %Actual length of data%
 actuallength(i)=abs(camera_data_w(n+1)-camera_data_w(n));
  t(i)=camera_time(n+1);
  k(i)=camera_data_w(n+1);
n=n+1;
i=i+1;

%Figure 1 plotting%

%figure('Name', 'Curve fitting with predicted next value for cubic and quadratic function')
plot(x, y, 'bp')
hold on
grid on
plot(xv, y_fit, '-r')
plot(xv, y_fit1,'-b')
plot(xv,y_fit2,'-m')
plot(u,v,'-g')
plot(u,v,'-g','linewidth',1)   %plot data with big width %
plot( t(i-1),predict(i-1),'-x')
plot( t(i-1),predictq(i-1),'-o')
plot( t(i-1),predictl(i-1),'-s')
plot(t(i-1),k(i-1),'-d')
xlabel('Camera Time')
ylabel('camera data')

end


% Data analysis cubic data error%
maxerror=max(predicterror);
minerror=min(predicterror);
averageerorr=mean(predicterror);
standardev=std(predicterror);

%Data analysis quadratic data error%
maxerrorq=max(predicterrorq);
minerrorq=min(predicterrorq);
averageerorrq=mean(predicterrorq);
standardevq=std(predicterrorq);


%Data analysis quadratic data error%
maxerrorl=max(predicterrorl);
minerrorl=min(predicterrorl);
averageerorrl=mean(predicterrorl);
standardevl=std(predicterrorl);

% % The error in actual length 

maxerror_actual=max(actuallength);
minerror_actual=min(actuallength);
averageerorr_actual=mean(actuallength);
standardev_actual=std(actuallength);

% %Figure 1 plotting%
% 
% figure('Name', 'Curve fitting with predicted next value for cubic and quadratic function')
% plot(x, y, 'bp')
% hold on
% grid on
% plot(xv, y_fit, '-r')
% plot(xv, y_fit1,'-b')
% plot(xv,y_fit2,'-m')
% plot(t,k,'-g')
% plot( t(i-1),predict(i-1),'-x')
% plot( t(i-1),predictq(i-1),'-o')
% plot( t(i-1),predictl(i-1),'-s')
% plot(t(i-1),k(i-1),'-d')
% xlabel('Camera Time')
% ylabel('camera data')
% hold off

%Figure 2 plotting%
figure('Name', 'Predicted length comaparision with Actual length')
plot(t,predictlength,'-r')
hold on
plot(t,predictlengthq,'-b')
plot(t,predictlengthl,'-m')
plot(t,actuallength,'-g')
xlabel('Camera Time')
ylabel('Length of samples from previous point')
hold off

%Figure 3 plotting%
 figure('Name', 'Predicted Error comparision')
 plot(t,predicterror,'-r')
 hold on
 plot(t,predicterrorq,'-b')
  plot(t,predicterrorl,'-g')
  plot(t,actuallength,'-k')
%   plot(t,noprediction,'-c')
 %plot(t,d,'-m')
 xlabel('Camera Time')
 ylabel('Predicted Error')


 
 
