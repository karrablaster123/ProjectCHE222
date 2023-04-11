clear all;

n = 0.5;
K1 = 150;
E_1 = 1.4; % eV
A_1 = 1.25 * 10^17;
k_B = 8.617 * 10^-5;
T0 = 80 + 273.15;
x_f_L = 0.1;
x_f_H = 0.16;
x_f = x_f_L;
K2 = 325;
x_i = 0.75; %0.45 0.25
E_2 = 0.8;
A_2 = 1e8;
z0 = 0.15;
a_L = 1.1;
a_H = 6.9;
a = a_L;
a0 = 1;
m = 1;
abserr = 1.0e-10;
relerr = 1.0e-8;
numpoints = 250;
t = linspace(0, 5000, numpoints);
options = odeset('RelTol', relerr,'AbsTol', abserr, 'NonNegative', 1);

p = [A_1, E_1, k_B, n, K1, A_2, E_2, m, K2, a, a0, z0];
x0 = [x_f; T0; x_i; z0];
%warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);

figure;
ax1 = subplot(3, 1, 1);
hold on;
ax2 = subplot(3, 1, 2);
hold on;
ax3 = subplot(3, 1, 3);
hold on;
figure;
ax4 = axes();
hold on;

T0 = 100 + 273.15;
x0 = [x_f; T0; x_i; z0];
xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
plot(ax1, xsol.y(2, :) - 273.15, xsol.y(1, :))
plot(ax2, xsol.y(2, :) - 273.15, xsol.y(3, :))
plot(ax3, xsol.y(2, :) - 273.15, xsol.y(4, :))
for i = 1:length(xsol.y(1, :))
    x = [xsol.y(1, i) xsol.y(2, i) xsol.y(3, i) xsol.y(4, i)];
    dTdt(:, i) = ODE(0, x, p);

end
semilogy(ax4, xsol.y(2, :) - 273.15, dTdt(2, :));
output = [xsol.y; dTdt];
writematrix(output, "Low SA.xlsx")

clear dTdt xsol;
x_f = x_f_H;
a = a_L;
T0 = 100 + 273.15;
p = [A_1, E_1, k_B, n, K1, A_2, E_2, m, K2, a, a0, z0];
x0 = [x_f; T0; x_i; z0];
xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);


plot(ax1, xsol.y(2, :) - 273.15, xsol.y(1, :))
plot(ax2, xsol.y(2, :) - 273.15, xsol.y(3, :))
plot(ax3, xsol.y(2, :) - 273.15, xsol.y(4, :))
for i = 1:length(xsol.y(1, :))
    x = [xsol.y(1, i) xsol.y(2, i) xsol.y(3, i) xsol.y(4, i)];
    dTdt(:, i) = ODE(0, x, p);

end
semilogy(ax4, xsol.y(2, :) - 273.15, dTdt(2, :));
output = [xsol.y; dTdt];
writematrix(output, "High SA.xlsx")
xlabel(ax1, "Temperature (C)");
xlabel(ax2, "Temperature (C)");
xlabel(ax3, "Temperature (C)");
ylabel(ax1, "x_f");
ylabel(ax2, "x_i");
ylabel(ax3, "z");
xlim(ax1, [100 220]);
xlim(ax2, [100 220]);
xlim(ax3, [100 220]);
xlim(ax4, [100 220]);
ylim(ax4, [10^-1 10^2]);
xlabel(ax4, "Temperature (C)");
ylim(ax1, [0 0.25]);
ylabel(ax4, "dT/dt");
legend(ax1, "Surace Area: " + ["Low", "High"]);
legend(ax2, "Surace Area: " + ["Low", "High"]);
legend(ax3, "Surace Area: " + ["Low", "High"]);
legend(ax4, "Surace Area: " + ["Low", "High"]);
title(ax1, "Variation of Heating with respect to Surface Area");
title(ax2, "Variation of Heating with respect to Surface Area");
title(ax3, "Variation of Heating with respect to Surface Area");
title(ax4, "Variation of Heating with respect to Surface Area");
set(ax4, 'YScale', 'log')


function f = ODE(t, x, p)

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

    f = [-exp(-E_1/(k_B*T))*A_1*(x_f^n);
        K1*exp(-E_1/(k_B*T))*A_1*(x_f^n) + (K2*A_2*x_i*(a/a0)*exp(-z/z0)*exp(-E_2/(k_B*T)));
        -exp(-E_2/(k_B*T))*A_2*(a/a0)*(x_i^m)*exp(-z/z0);
        A_2*exp(-E_2/(k_B*T))*(x_i^m)*exp(-z/z0)];

end