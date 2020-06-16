tspan = [0 20];
[t,y] = ode45(@vdp1, tspan, [2;0]);
plot(t,y(:,1),'-o',t,y(:,2),'-o');
title('Solution of van der Pol equation (\mu=1) with ODE 45')
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')