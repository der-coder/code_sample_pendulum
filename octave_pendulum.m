%% Simple pendulum simulation
%% Author: Isaac Ayala
%% MIT License (c) 2020 

function octave_pendulum
  
%% Clear all stored or cached data
clear
clc
clf
close all

%% Set parameters
l = 0.15;
m = 0.3;
k = 0.35;
g = 9.81;
theta = pi/2;

x0 = [theta 0]; 
max_time = 10;
tspan = 0:0.01:max_time;

%% Solve system 
[t, sol] = ode45(@(t,y)pendulum(t,y,l,m,k, g), tspan, x0);

%% Calculate x and y
x = l * sin(sol(:, 1));
y = - l * cos(sol(:, 1));

%% Plot results
%% x-y phase plot
subplot(2, 2, 1)
plot(x, y,'k','linewidth', 2)
legend({'linear'},'location', 'north', 'orientation', 'horizontal')
legend('boxoff')

%% Cartesian coordinate plot
subplot(2, 2, 2)

xlabel('t')
ylabel('state')

plot(t, x,'k','linewidth', 2)
hold on
plot(t, y,'k--','linewidth', 2)
 
legend({'x','  y'},'location', 'northeast', 'orientation', 'horizontal')
legend('boxoff')

%% theta-omega phase plot
subplot(2, 2, 3)
plot(sol(:,1), sol(:,2),'k','linewidth', 2) 
legend({'angular'},'location', 'south', 'orientation', 'horizontal')
legend('boxoff')

%% Angular coordinate plot
subplot(2, 2, 4)
xlabel('t')
ylabel('state')

plot(t, sol(:,1),'k','linewidth', 2)
hold on
plot(t, sol(:,2),'k--','linewidth', 2)
 
legend({'\theta','  \omega'},'location', 'northeast', 'orientation', 'horizontal')
legend('boxoff')



set(gcf,'Color',[1 1 1])
set(gcf, 'fontsize', 12)


end

%% System dynamics equations
function dx = pendulum(~,x,l,m,k, g)
dx(1) = x(2);
dx(2) = -(g/l)*sin(x(1)) - (k/m)*x(2);
dx = dx';
end