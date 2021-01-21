clear all;
clc

N = 50;
h = 0.05;

u = 0.5;
delta = 0;

for i = 1:N+1,
   t = (i-1)*h;      % time
 
   delta_dot = 10*u - 10*delta;
             
   table(i,:) = [t delta' ];  % store data in table
   
   % Euler integration
   delta = delta + h*delta_dot;
   
   
 
end 

t = table(:,1);
x = table(:,2);

figure(1)
figure(gcf)
subplot(311)
plot(t(:),x(:),'b','Linewidth',1);
xlabel('$time$ $[s]$','Interpreter','latex');
ylabel('$delta$ $[rad]$','Interpreter','latex');