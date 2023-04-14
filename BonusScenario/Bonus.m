clear all;
%{
Bonus Scenario: Group 21
This code aims to replicate Figure 7 from the paper. It is, in a sense, the
culmination of all the previous scenarios to provide a reasonable
approximation to the empirical self-heating of a lithium ion battery
%}

%Parameter Definition
n = 0.5;
K1 = 150;
E_1 = 1.4;
A_1 = 1.25e17;
k_B = 8.617e-5;
T0 = 80 + 273.15;
x_f = 0.1;
K2 = 325;
x_i = 0.75;
E_2 = 0.8;
A_2 = 1e8;
z0 = 0.15;
a_L = 1.1;
a = a_L;
a0 = 1;
m = 1;

%Solver Parameters
abserr = 1.0e-8;
relerr = 1.0e-6;
numpoints = 250;
t = linspace(0, 5000, numpoints);
options = odeset('RelTol',relerr,'AbsTol',abserr, 'NonNegative', 1);

%Parameter Array
p = [A_1, E_1, k_B, n, K1, A_2, E_2, m, K2, a, a0, z0];
%Initial Value Array
x0 = [x_f; T0; x_i; z0];
%Solve
xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);

%plot
figure;
ax1 = subplot(3, 1, 1);
hold on;
ax2 = subplot(3, 1, 2);
hold on;
ax3 = subplot(3, 1, 3);
hold on;
plot(ax1, xsol.y(2, :) - 273.15, xsol.y(1, :), 'LineWidth', 3)
plot(ax2, xsol.y(2, :) - 273.15, xsol.y(3, :), 'LineWidth', 3)
plot(ax3, xsol.y(2, :) - 273.15, xsol.y(4, :), 'LineWidth', 3)

%Get value of dTdt at temperature T
for i = 1:length(xsol.y(1, :))
    x = [xsol.y(1, i) xsol.y(2, i) xsol.y(3, i) xsol.y(4, i)];
    dTdt(:, i) = ODE(0, x, p); %#ok<*SAGROW> 

end
figure;
ax4 = axes();
hold on;
semilogy(ax4, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth', 3);

%Now solving for T_0 = 100
clear xsol dTdt;

%Change T0
T0 = 100 + 273.15;
%Reassign Initial Value 
x0 = [x_f; T0; x_i; z0];
%Solve
xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);

%Plot
plot(ax1, xsol.y(2, :) - 273.15, xsol.y(1, :), 'LineWidth', 3)
plot(ax2, xsol.y(2, :) - 273.15, xsol.y(3, :), 'LineWidth', 3)
plot(ax3, xsol.y(2, :) - 273.15, xsol.y(4, :), 'LineWidth', 3)

%Get value of dT/dt at temperature T
for i = 1:length(xsol.y(1, :))
    x = [xsol.y(1, i) xsol.y(2, i) xsol.y(3, i) xsol.y(4, i)];
    dTdt(:, i) = ODE(0, x, p);

end

semilogy(ax4, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth', 3);

%Repeated as prior
clear xsol dTdt;
T0 = 115 + 273.15;
x0 = [x_f; T0; x_i; z0];
xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
plot(ax1, xsol.y(2, :) - 273.15, xsol.y(1, :), 'LineWidth', 3)
plot(ax2, xsol.y(2, :) - 273.15, xsol.y(3, :), 'LineWidth', 3)
plot(ax3, xsol.y(2, :) - 273.15, xsol.y(4, :), 'LineWidth', 3)
for i = 1:length(xsol.y(1, :))
    x = [xsol.y(1, i) xsol.y(2, i) xsol.y(3, i) xsol.y(4, i)];
    dTdt(:, i) = ODE(0, x, p);

end
semilogy(ax4, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth', 3);

%Label and scale plots
xlim([80 220]);
xlabel(ax1, "Temperature (C)");
xlabel(ax2, "Temperature (C)");
xlabel(ax3, "Temperature (C)");
ylabel(ax1, "x_f");
ylabel(ax2, "x_i");
ylabel(ax3, "z");
xlim(ax1, [80 220]);
xlim(ax2, [80 220]);
xlim(ax3, [80 220]);
xlabel(ax4, "Temperature (C)");
ylabel(ax4, "dT/dt");
legend(ax1, "T0 = " + string([80 100 115]));
legend(ax2, "T0 = " + string([80 100 115]));
legend(ax3, "T0 = " + string([80 100 115]));
legend(ax4, "T0 = " + string([80 100 115]));
title(ax1, "x_f versus Temperature (Fig. 7)");
title(ax2, "x_i versus Temperature (Fig. 7)");
title(ax3, "z versus Tempearture (Fig. 7)");
title(ax4, "dT/dt versus Temperauture (Fig. 7)");
set(ax4, 'YScale', 'log')
set(ax1, 'FontSize', 15)
set(ax2, 'FontSize', 15)
set(ax3, 'FontSize', 15)
set(ax4, 'FontSize', 15)


%Function of Equations
function f = ODE(~, x, p)

    x_f = x(1);
    T = x(2);
    x_i = x(3);
    z = x(4);
    A_1 = p(1);
    E_1 = p(2);
    k_B = p(3);
    n = p(4);
    K1 = p(5);
    A_2 = p(6);
    E_2 = p(7);
    m = p(8);
    K2 = p(9);
    a = p(10);
    a0 = p(11);
    z0 = p(12);

    if x_f < 0
        x_f = 0;
    end

    f = [-exp(-E_1/(k_B*T))*A_1*(x_f^n);
        K1*exp(-E_1/(k_B*T))*A_1*(x_f^n) + (K2*A_2*x_i*(a/a0)*exp(-z/z0)*exp(-E_2/(k_B*T)));
        -exp(-E_2/(k_B*T))*A_2*(a/a0)*(x_i^m)*exp(-z/z0);
        A_2*exp(-E_2/(k_B*T))*(x_i^m)*exp(-z/z0)];

end