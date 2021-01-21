close all;
clear all;
clc;

%% Initializations and constants

h  = 0.05;    % sampling time [s]
N = 6000;     % no. of samples

flag = 1;     % 1: nonlinear AOA, 0: linear AOA

% UAV parameters
b = 2.1;                 % wing span (m)
l = 0.79;                % fuselage length (m)
m = 2.5;                 % mass (kg)
S = 0.8;                 % wing area (m2)
c = 0.357;               % mean aerodynamic chord (m)
g = 9.81;
rho = 1.225;             % air density (kg/m3)

% Inertia matrix
Ix = 1.229; Iy = 0.1702; Iz = 0.8808;
Ixy = 0; Ixz = 0.9343; Iyz = 0;        % ref Gryte thesis

% I = [Ix, -Ixy, -Ixz;
%      -Ixy, Iy, -Iyz;
%      -Ixz, -Iyz, Iz];

% Inertia relationships
L = Ix*Iz-Ixz^2;
L1 = (Ixz*(Ix-Iz+Iz))/L;
L2 = (Iz*(Iz-Iy)+Ixz^2)/L;
L3 = Iz/L;
L4 = Ixz/L;
L5 = (Iz-Ix)/Iy;
L6 = Ixz/Iy;
L7 = ((Ix-Iy)*Ix+Ixz^2)/L;
L8 = Ix/L;

% Initialize states
pos = [0 0 0]';            % position [x y z] in NED
vel = [1 1 1]';            % velocity [u v w] in BODY
ang = [0 0 0]';            % euler angles [phi theta psi] in VEHICLE
angrate = [0 0 0]';        % angular rates [p q r] in BODY

% F = [1 1 1]'; % forces [X Y Z]
% M = [0 0 0]'; % moments [L M N]

% Wind speed
Vw = [3, 1, 2]';

delta_e = 0;
delta_a = 0;

%% FOR LOOP

for i = 1:N+1,
   t = (i-1)*h;      % time
   
   u_r = vel(1) - Vw(1);
   v_r = vel(2) - Vw(2);
   w_r = vel(3) - Vw(3);
   
   Va = sqrt(u_r^2 + v_r^2 + w_r^2);
   AOA = atan2(w_r,u_r);
   beta = asin(v_r/Va);
   
   % Gravity force
   fg = m*g*[-sin(ang(2));
            cos(ang(2))*sin(ang(1));
            cos(ang(2))*cos(ang(2))];
   % Aerodynamics
   [fa,ma] = aero(Va,AOA,beta,angrate,delta_e,delta_a,flag);
   
   F = fg + fa;
   M = [0 0 0]';  % Fungerer ikke med aerodynamic moment...
        
   
   % Rotation matrices
   R_trans = Rzyx(ang(1),ang(2),ang(3));  % BODY to NED 

   R_ang = [1, sin(ang(1))*tan(ang(2)), cos(ang(1))*tan(ang(2));
   0, cos(ang(1)), -sin(ang(1));
   0, sin(ang(1))/cos(ang(2)), cos(ang(1))/cos(ang(2))]; % BODY to VEHICLE

   % Kinematics (geometrical aspect of motion)
   pos_dot = R_trans*vel;
   ang_dot = R_ang*angrate;
   
   % Kinetics (forces causing the motion)
   vel_dot = [angrate(3)*vel(2)-angrate(2)*vel(3);
            angrate(1)*vel(3)-angrate(3)*vel(1);
            angrate(2)*vel(1)-angrate(1)*vel(2)] + 1/m * F;
        
   angrate_dot = [L1*angrate(1)*angrate(2)-L2*angrate(2)*angrate(3);
                L5*angrate(1)*angrate(3)-L6*(angrate(1)^2-angrate(3)^2);
                L7*angrate(1)*angrate(2)-L1*angrate(2)*angrate(3)] +... 
                [L3*M(1) + L4*M(3);
                 1/Iy*M(2);
                 L4*M(1) + L8*M(3)];
             
   table(i,:) = [t pos' vel' ang' angrate'];  % store data in table
   
   % Euler integration
   pos = pos + h*pos_dot;
   vel = vel + h*vel_dot;
   ang = ang + h*ang_dot;
   angrate = angrate + h*angrate_dot;
   
 
end 

t = table(:,1);
x = table(:,2);
y = table(:,3);
z = table(:,4);
u = table(:,5);
v = table(:,6);
w = table(:,7);
phi = table(:,8);
theta = table(:,9);
psi = table(:,10);
p = table(:,11);
q = table(:,12);
r = table(:,13);

%% Plots

% Position in NED
figure(1)
figure(gcf)
subplot(311)
plot(t(:),x(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$x$ $[m]$','Interpreter','latex');
legend('North')

subplot(312)
plot(t(:),y(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$y$ $[m]$','Interpreter','latex');
legend('East')

subplot(313)
plot(t(:),-z(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$z$ $[m]$','Interpreter','latex');
legend('Down')

% Velocities in BODY
figure(2)
figure(gcf)
subplot(311)
plot(t(:),u(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$u$ $[m/s]$','Interpreter','latex');
legend('u-vel')

subplot(312)
plot(t(:),v(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$v$ $[m/s]$','Interpreter','latex');
legend('v-vel')

subplot(313)
plot(t(:),w(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$w$ $[m/s]$','Interpreter','latex');
legend('w-vel')

% Position in NE
figure(3)
figure(gcf)
plot(y(:),x(:),'b','Linewidth',1); 
xlabel('$East$ $[m]$','Interpreter','latex');
ylabel('$North$ $[m]$','Interpreter','latex');
% legend('')
grid on
axis equal


