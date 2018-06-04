%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%> @file view_data.m  
%>
%> @brief Script for plotting the data from the zero-velocity aided inertial
%> navigations system.
%>
%> @authors Isaac Skog, John-Olof Nilsson
%> @modified by G.V. Prateek
%> @copyright Copyright (c) 2013 OpenShoe, ISC License (open source)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global simdata;


% Close all windows
close all

%% Straight Path
% Comment this next four lines if it is not U Path trajectory
a=[.3 .3 ];
b=[-.4243 34.9757];
c=[-.3 -.3 ];
d=[-.4243 34.9757];

%% L Path Trajectory
% Comment this next four lines if it is not U Path trajectory
a=[.3 .3 -23.7 ];
b=[-.4243 34.9757 34.9757 ];
c=[-.3 -.3 -23.1 ];
d=[-.4243 34.3757 34.3757 ];


%% U Path trajectory
% Comment this next four lines if it is not U Path trajectory
a=[.3 .3 -23.7 -23.7];
b=[-.4243 34.9757 34.9757 -.4243];
c=[-.3 -.3 -23.1 -23.1];
d=[-.4243 34.3757 34.3757 -.4243];


% Generate a vector with the time scale
N=min(length(cov_l),length(cov_r));
t=0:simdata.Ts:(N-1)*simdata.Ts;
foldername='';


%% Plot the IMU data
figure(1)
clf
subplot(2,1,1)
plot(t,u1(1:3,1:N)')
xlabel('time [s]')
ylabel('Specific force [m/s^2]')
title('Specific force (accelerometer) measurements')
legend('x-axis','y-axis','z-axis')
box on
grid on

subplot(2,1,2)
plot(t,u1(4:6,1:N)'*180/pi)
xlabel('time [s]')
ylabel('Angular rate  [deg/s]')
title('Angular rate measurements')
legend('x-axis','y-axis','z-axis')
box on
grid on

figure(2)
clf
subplot(2,1,1)
plot(t,u2(1:3,1:N)')
xlabel('time [s]')
ylabel('Specific force [m/s^2]')
title('Specific force (accelerometer) measurements')
legend('x-axis','y-axis','z-axis')
box on
grid on

subplot(2,1,2)
plot(t,u2(4:6,1:N)'*180/pi)
xlabel('time [s]')
ylabel('Angular rate  [deg/s]')
title('Angular rate measurements')
legend('x-axis','y-axis','z-axis')
box on
grid on


%% Plot the trajectory in the horizontal plane
figure(3)
clf
plot(x_h_r(1,:),x_h_r(2,:),'b')
hold on
plot(x_h_l(1,:),x_h_l(2,:),'r')
hold on
plot(x_h_r(1,1),x_h_r(2,1),'ms')
hold on
%Girish Starts
plot(a,b,'k')
hold on
plot(c,d,'k')
%Girish Ends
title('Trajectory')
legend('Right Leg Trajectory','Left Leg Trajectory','Start point','Actual Path')
xlabel('x [m]')
ylabel('y [m]')
axis equal
grid on
box on
saveas(gcf,[foldername 'trajectory.fig'])


%% Plot the height profile, the speed and when ZUPTs were applied

figure(4)
clf
subplot(3,1,1)
plot(t,-x_h_r(3,:))
title('Right foot Heigth')
xlabel('time [s]')
ylabel('z [m]')
grid on
box on


subplot(3,1,2)
plot(t,sqrt(sum(x_h_r(4:6,:).^2)))
title('Right foot Speed')
xlabel('time [s]')
ylabel('|v| [m/s]')
grid on
box on

subplot(3,1,3)
stem(t,zupt1(1:N))
title('Right foot Zupt applied')
xlabel('time [s]')
ylabel('on/off')
grid on
box on
saveas(gcf,[foldername 'zupt_r.fig'])

%%

figure(5)
clf
subplot(3,1,1)
plot(t,-x_h_l(3,:))
title('Left foot Heigth')
xlabel('time [s]')
ylabel('z [m]')
grid on
box on


subplot(3,1,2)
plot(t,sqrt(sum(x_h_l(4:6,:).^2)))
title('Left foot Speed')
xlabel('time [s]')
ylabel('|v| [m/s]')
grid on
box on

subplot(3,1,3)
stem(t,zupt2(1:N))
title('Left foot Zupt applied')
xlabel('time [s]')
ylabel('on/off')
grid on
box on
saveas(gcf,[foldername 'zupt_l.fig'])


%% Plot the attitude
%% function of time

figure(6)
clf
plot(t,unwrap(x_h_r(7:9,:)')*180/pi)
title('Right foot Attitude')
xlabel('time [s]')
ylabel('Angle [deg]')
legend('Roll','Pitch','Yaw')
grid on
box on
saveas(gcf,[foldername 'rpy_r.fig'])

%%
figure(7)
clf
plot(t,unwrap(x_h_l(7:9,:)')*180/pi)
title('Left foot Attitude')
xlabel('time [s]')
ylabel('Angle [deg]')
legend('Roll','Pitch','Yaw')
grid on
box on
saveas(gcf,[foldername 'rpy_l.fig'])


%% Plot the diagonal elements of the filter covariance matrices as a
%% function of time

figure(8)
clf

subplot(3,1,1)
plot(t,sqrt(cov_r(1:3,:))')
title('Right foot Position covariance')
ylabel('sqrt(cov) [m]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on

subplot(3,1,2)
plot(t,sqrt(cov_r(4:6,:))')
title('Right foot Velocity covariance')
ylabel('sqrt(cov) [m/s]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on

subplot(3,1,3)
plot(t,sqrt(cov_r(7:9,:))'*180/pi)
title('Right foot Heading covariance')
ylabel('sqrt(cov) [deg]')
xlabel('time [s]')
legend('Roll', 'Pitch','Yaw')
grid on
box on
saveas(gcf,[foldername 'cov_r.fig'])

%%
figure(9)
clf

subplot(3,1,1)
plot(t,sqrt(cov_l(1:3,:))')
title('Left foot Position covariance')
ylabel('sqrt(cov) [m]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on

subplot(3,1,2)
plot(t,sqrt(cov_l(4:6,:))')
title('Left foot Velocity covariance')
ylabel('sqrt(cov) [m/s]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on

subplot(3,1,3)
plot(t,sqrt(cov_l(7:9,:))'*180/pi)
title('Left foot Heading covariance')
ylabel('sqrt(cov) [deg]')
xlabel('time [s]')
legend('Roll', 'Pitch','Yaw')
grid on
box on
saveas(gcf,[foldername 'cov_l.fig'])


%% If the filter also estimates the sensor biases and/or the scalefactor,
%% plot these now

if (strcmp(simdata.scalefactors,'on') && strcmp(simdata.biases,'on'))
    %    Both scale and bias errors included
    figure(6)
    clf
    subplot(2,1,1)
    plot(t,x_h(10:12,:)')
    legend('x-axis','y-axis','z-axis')
    title('Accelerometer bias errors')
    xlabel('time [s]')
    ylabel('Bias [m/s^2]')
    grid on
    box on
    
    subplot(2,1,2)
    plot(t,x_h(13:15,:)'*180/pi)
    legend('x-axis','y-axis','z-axis')
    title('Gyroscope bias errors')
    xlabel('time [s]')
    ylabel('Bias [deg/s]')
    box on
    grid on
    
    figure(7)
    clf
    subplot(2,1,1)
    plot(t,x_h(16:18,:)')
    legend('x-axis','y-axis','z-axis')
    title('Accelerometer scale factor errors')
    xlabel('time [s]')
    box on
    grid on
    
    subplot(2,1,2)
    plot(t,x_h(19:21,:)')
    legend('x-axis','y-axis','z-axis')
    title('Gyroscope scale factor errors')
    xlabel('time [s]')
    box on
    grid on
    
    
elseif strcmp(simdata.scalefactors,'on') && strcmp(simdata.biases,'off')
    %    Scale errors included
    figure(6)
    clf
    subplot(2,1,1)
    plot(t,x_h(10:12,:)')
    legend('x-axis','y-axis','z-axis')
    title('Accelerometer scale factor errors')
    xlabel('time [s]')
    box on
    grid on
    
    subplot(2,1,2)
    plot(t,x_h(13:15,:)'*180/pi)
    legend('x-axis','y-axis','z-axis')
    title('Gyroscope scale factor errors')
    xlabel('time [s]')
    box on
    grid on
    
elseif strcmp(simdata.scalefactors,'off') && strcmp(simdata.biases,'on')
    %    Bias errors included
    figure(6)
    clf
    subplot(2,1,1)
    plot(t,x_h(10:12,:)')
    legend('x-axis','y-axis','z-axis')
    title('Accelerometer bias errors')
    xlabel('time [s]')
    ylabel('Bias [m/s^2]')
    grid on
    box on
    
    subplot(2,1,2)
    plot(t,x_h(13:15,:)'*180/pi)
    legend('x-axis','y-axis','z-axis')
    title('Gyroscope bias errors')
    xlabel('time [s]')
    ylabel('Bias [deg/s]')
    box on
    grid on
end








