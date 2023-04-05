clear all;

n = 0.5;
K = 150; % degrees Celsius
E_1 = 1.4; % eV
A_1 = 1 * 10^17;
k_B = 8.617 * 10^-5;
T0 = 80 + 273.15;
x_f = 0.1;
generate_array = @(CtrPt) [0.1*CtrPt, 0.5*CtrPt, CtrPt, 2.5*CtrPt, 5*CtrPt];
abserr = 1.0e-8;
relerr = 1.0e-6;
numpoints = 250;
t = linspace(0, 5000, numpoints);
x0 = [x_f, T0];
options = odeset('RelTol',relerr,'AbsTol',abserr, 'NonNegative', 1);

%Varying A1
A_1_array = generate_array(A_1);

figure;

ax1 = subplot(2, 1, 1);
hold on;
ax2 = subplot(2, 1, 2);
hold on;
figure; ax3 = gca;
hold on;
for i = 1:length(A_1_array)
    p = [A_1_array(i), E_1, k_B, n, K];
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
    plot(ax1, xsol.x, xsol.y(1, :));
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p); %#ok<*SAGROW> 
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :));
    clear dTdt;

end
title(ax1, "Dynamics when varying A");
title(ax2, "Dynamics when varying A");
xlabel(ax1, "Time");
ylabel(ax1, "x_f");
xlabel(ax2, "Time");
ylabel(ax2, "Temperature (C)");
legend(ax1, "A = " + string(A_1_array));
legend(ax2, "A = " + string(A_1_array));
title(ax3, "Dynamics when varying A");
xlabel(ax3, "Temperature");
ylabel(ax3, "dT/dt");
legend(ax3, "A = " + string(A_1_array));
set(ax3, 'YScale', 'log')
ylim(ax3, [10^-4 10^3])


%Varying K
clear xsol ax1 ax2 dTdt ax3;
K_array = generate_array(K);
figure;
ax1 = subplot(2, 1, 1);
hold on;
ax2 = subplot(2, 1, 2);
hold on;
figure; ax3 = gca;
hold on;
for i = 1:length(K_array)
    p = [A_1, E_1, k_B, n, K_array(i)];
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
    plot(ax1, xsol.x, xsol.y(1, :));
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15);
    %disp(xsol.y(2, :))
    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p); %#ok<*SAGROW> 
    end
    %disp(xsol.y(2, :))
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :)); clear dTdt;
end
title(ax1, "Dynamics when varying K");
title(ax2, "Dynamics when varying K");
xlabel(ax1, "Time");
ylabel(ax1, "x_f");
xlabel(ax2, "Time");
ylabel(ax2, "Temperature (C)");
legend(ax1, "K = " + string(K_array));
legend(ax2, "K = " + string(K_array));
title(ax3, "Dynamics when varying K");
xlabel(ax3, "Temperature");
ylabel(ax3, "dT/dt");
legend(ax3, "K = " + string(K_array));
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')

%Varying E_1
clear xsol ax1 ax2 dTdt ax3;
figure;
ax1 = subplot(2, 1, 1);
hold on;
ax2 = subplot(2, 1, 2);
hold on;
figure; ax3 = gca;
hold on;
E_1_array = linspace(1.38, E_1*5, 5);
for i = 1:length(E_1_array)
    p = [A_1, E_1_array(i), k_B, n, K];
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
    plot(ax1, xsol.x, xsol.y(1, :));
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :)); clear dTdt;
end
title(ax1, "Dynamics when varying E_1");
title(ax2, "Dynamics when varying E_1");
xlabel(ax1, "Time");
ylabel(ax1, "x_f");
xlabel(ax2, "Time");
ylabel(ax2, "Temperature (C)");
legend(ax1, "E_1 = " + string(E_1_array));
legend(ax2, "E_1 = " + string(E_1_array));
title(ax3, "Dynamics when varying E_1");
xlabel(ax3, "Temperature");
ylabel(ax3, "dT/dt");
legend(ax3, "E_1 = " + string(E_1_array));
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')

%Varying n
clear xsol ax1 ax2 dTdt ax3;
n_array = generate_array(n);
figure;
ax1 = subplot(2, 1, 1);

hold on;
ax2 = subplot(2, 1, 2);

hold on;
figure; ax3 = gca;
hold on;
for i = 1:length(n_array)
    p = [A_1, E_1, k_B, n_array(i), K];
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
    plot(ax1, xsol.x, xsol.y(1, :));
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :)); clear dTdt;
end
title(ax1, "Dynamics when varying n");
title(ax2, "Dynamics when varying n");
xlabel(ax1, "Time");
ylabel(ax1, "x_f");
xlabel(ax2, "Time");
ylabel(ax2, "Temperature (C)");
legend(ax1, "n = " + string(n_array));
legend(ax2, "n = " + string(n_array));
title(ax3, "Dynamics when varying n");
xlabel(ax3, "Temperature");
ylabel(ax3, "dT/dt");
legend(ax3, "n = " + string(n_array));
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')


%Varying x_f
clear xsol ax1 ax2 dTdt ax3;
x_f_array = generate_array(x_f);
figure;
ax1 = subplot(2, 1, 1);
hold on;
ax2 = subplot(2, 1, 2);
hold on;
figure; ax3 = gca;
hold on;
for i = 1:length(x_f_array)
    p = [A_1, E_1, k_B, n, K];
    x0 = [x_f_array(i), T0];
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
    plot(ax1, xsol.x, xsol.y(1, :));
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15);
    
    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :)); clear dTdt;

end
title(ax1, "Dynamics when varying x_f0");
title(ax2, "Dynamics when varying x_f0");
xlabel(ax1, "Time");
ylabel(ax1, "x_f");
xlabel(ax2, "Time");
ylabel(ax2, "Temperature (C)");
legend(ax1, "x_f = " + string(x_f_array));
legend(ax2, "x_f = " + string(x_f_array));
xlabel(ax3, "Temperature");
title(ax3, "Dynamics when varying x_f0");
ylabel(ax3, "dT/dt");
legend(ax3, "x_f = " + string(x_f_array));
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')

%Varying T0
clear xsol ax1 ax2 dTdt ax3;
T0_array = linspace(20+273.15, 150+273.15, 7);
figure;
ax1 = subplot(2, 1, 1);
hold on;
ax2 = subplot(2, 1, 2);
hold on;
figure; ax3 = gca;
hold on;
for i = 1:length(T0_array)
    p = [A_1, E_1, k_B, n, K];
    x0 = [x_f, T0_array(i)];
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);
    plot(ax1, xsol.x, xsol.y(1, :));
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :)); clear dTdt;
end
title(ax1, "Dynamics when varying T0");
title(ax2, "Dynamics when varying T0");
xlabel(ax1, "Time");
ylabel(ax1, "x_f");
xlabel(ax2, "Time");
ylabel(ax2, "Temperature (C)");
legend(ax1, "T0 = " + string(T0_array - 273.15));
legend(ax2, "T0 = " + string(T0_array - 273.15));
title(ax3, "Dynamics when varying T0");
xlabel(ax3, "Temperature");
ylabel(ax3, "dT/dt");
legend(ax3, "T0 = " + string(T0_array - 273.15));
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')


function f = ODE(~, x, p)
        x1 = x(1);
        T1 = x(2);
        A_1 = p(1);
        E_1 = p(2);
        k_B = p(3);
        n = p(4);
        K = p(5);

        f = [-exp(-E_1/(k_B*T1))*A_1*(x1^n); K*exp(-E_1/(k_B*T1))*A_1*(x1^n)];

end
    