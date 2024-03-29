clear all;

%{
Scenario 5: Group 21
Kush Patel - 1007678826
Sachiel Malik - 1007692098
Gurjas Singh Chawla - 1007872189
Uddhava Swaminathan - 1007511122

The function from Scenario 1 has been taken and parameters have been varied
to provide insight on their effect on the self-heating of the battery.
%}
% Parameter Defintion
n = 0.5;
K = 150; 
E_1 = 1.4;
A_1 = 1.25 * 10^17;
k_B = 8.617 * 10^-5;
T0 = 80 + 273.15;
x_f = 0.1;

%Generate array of variables
generate_array = @(CtrPt) [0.1*CtrPt, 0.5*CtrPt, CtrPt, 2.5*CtrPt, 5*CtrPt];

%Absolute and relative error for solver
abserr = 1.0e-8;
relerr = 1.0e-6;
numpoints = 250;

%Time Range
t = linspace(0, 5000, numpoints);

%Setting initial value
x0 = [x_f, T0];

%Setting options
options = odeset('RelTol',relerr,'AbsTol',abserr, 'NonNegative', 1);

%Varying A1
A_1_array = generate_array(A_1);

%A1 Figure and axes setup
figure;
ax1 = subplot(2, 1, 1);
hold on;
ax2 = subplot(2, 1, 2);
hold on;
figure; ax3 = gca;
hold on;

%Solving for each value
for i = 1:length(A_1_array)
    
    %Parameters to be passed into the function
    p = [A_1_array(i), E_1, k_B, n, K];
    
    %Solve ODE
    xsol = ode15s(@(t,x)ODE(t, x, p), t, x0, options);

    %Plot xi and temperature
    plot(ax1, xsol.x, xsol.y(1, :), 'LineWidth', 2);
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15, 'LineWidth', 2);

    %For each temperature, we resolve the equation to get dT/dt at T
    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p); %#ok<*SAGROW> 
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth', 4);
    clear dTdt;

end
%Plot Setup
title(ax1, "Dynamics when varying A");
title(ax2, "Dynamics when varying A");
xlabel(ax1, "Time (min)", 'FontSize', 20);
ylabel(ax1, "x_f", 'FontSize', 20);
xlabel(ax2, "Time (min)", 'FontSize', 20);
ylabel(ax2, "Temperature (C)", 'FontSize', 20);
legend(ax1, "A = " + string(A_1_array), 'FontSize', 14, 'Location','best');
legend(ax2, "A = " + string(A_1_array), 'FontSize', 14, 'Location','best');
title(ax3, "Dynamics when varying A", 'FontSize', 20);
xlabel(ax3, "Temperature (C)", 'FontSize', 20);
ylabel(ax3, "dT/dt", 'FontSize', 20);
legend(ax3, "A = " + string(A_1_array), 'FontSize', 14, 'Location','best');
set(ax3, 'YScale', 'log')
ylim(ax3, [10^-4 10^3])
set(ax3,"FontSize", 20); set(ax1, "FontSize", 20); set(ax2, "FontSize", 20);


%Clear variables for reassigning
clear xsol ax1 ax2 dTdt ax3;
%The pattern above now repeats for the rest of code
%Varying K
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
    plot(ax1, xsol.x, xsol.y(1, :), 'LineWidth', 2);
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15, 'LineWidth', 2);
    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p); %#ok<*SAGROW> 
    end
    
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth',4); clear dTdt;
end

title(ax1, "Dynamics when varying K");
title(ax2, "Dynamics when varying K");
xlabel(ax1, "Time (min)", 'FontSize', 20);
ylabel(ax1, "x_f", 'FontSize', 20);
xlabel(ax2, "Time (min)", 'FontSize', 20);
ylabel(ax2, "Temperature (C)", 'FontSize', 20);
legend(ax1, "K = " + string(K_array), 'FontSize', 14, 'Location','best');
legend(ax2, "K = " + string(K_array), 'FontSize', 14, 'Location','best');
title(ax3, "Dynamics when varying K", 'FontSize', 20);
xlabel(ax3, "Temperature (C)", 'FontSize', 20);
ylabel(ax3, "dT/dt", 'FontSize', 20);
legend(ax3, "K = " + string(K_array), 'FontSize', 14, 'Location','best');
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')
set(ax3,"FontSize", 20); set(ax1, "FontSize", 20); set(ax2, "FontSize", 20);

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
    plot(ax1, xsol.x, xsol.y(1, :), 'LineWidth', 2);
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15, 'LineWidth', 2);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth',4); clear dTdt;
end

title(ax1, "Dynamics when varying E_1");
title(ax2, "Dynamics when varying E_1");
xlabel(ax1, "Time (min)", 'FontSize', 20);
ylabel(ax1, "x_f", 'FontSize', 20);
xlabel(ax2, "Time (min)", 'FontSize', 20);
ylabel(ax2, "Temperature (C)", 'FontSize', 20);
legend(ax1, "E_1 = " + string(E_1_array), 'FontSize', 14, 'Location','best');
legend(ax2, "E_1 = " + string(E_1_array), 'FontSize', 14, 'Location','best');
title(ax3, "Dynamics when varying E_1", 'FontSize', 20);
xlabel(ax3, "Temperature (C)", 'FontSize', 20);
ylabel(ax3, "dT/dt", 'FontSize', 20);
legend(ax3, "E_1 = " + string(E_1_array), 'FontSize', 14, 'Location','best');
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')
set(ax3,"FontSize", 20); set(ax1, "FontSize", 20); set(ax2, "FontSize", 20);

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
    plot(ax1, xsol.x, xsol.y(1, :), 'LineWidth', 2);
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15, 'LineWidth', 2);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth',4); clear dTdt;
end

title(ax1, "Dynamics when varying n");
title(ax2, "Dynamics when varying n");
xlabel(ax1, "Time (min)", 'FontSize', 20);
ylabel(ax1, "x_f", 'FontSize', 20);
xlabel(ax2, "Time (min)", 'FontSize', 20);
ylabel(ax2, "Temperature (C)", 'FontSize', 20);
legend(ax1, "n = " + string(n_array), 'FontSize', 14, 'Location','best');
legend(ax2, "n = " + string(n_array), 'FontSize', 14, 'Location','best');
title(ax3, "Dynamics when varying n", 'FontSize', 20);
xlabel(ax3, "Temperature (C)", 'FontSize', 20);
ylabel(ax3, "dT/dt", 'FontSize', 20);
legend(ax3, "n = " + string(n_array), 'FontSize', 14, 'Location','best');
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')
set(ax3,"FontSize", 20); set(ax1, "FontSize", 20); set(ax2, "FontSize", 20);


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
    plot(ax1, xsol.x, xsol.y(1, :), 'LineWidth', 2);
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15, 'LineWidth', 2);
    
    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth',4); clear dTdt;

end

title(ax1, "Dynamics when varying x_f0");
title(ax2, "Dynamics when varying x_f0");
xlabel(ax1, "Time (min)", 'FontSize', 20);
ylabel(ax1, "x_f", 'FontSize', 20);
xlabel(ax2, "Time (min)", 'FontSize', 20);
ylabel(ax2, "Temperature (C)", 'FontSize', 20);
legend(ax1, "x_f = " + string(x_f_array), 'FontSize', 14, 'Location','best');
legend(ax2, "x_f = " + string(x_f_array), 'FontSize', 14, 'Location','best');
xlabel(ax3, "Temperature (C)", 'FontSize', 20);
title(ax3, "Dynamics when varying x_f0", 'FontSize', 20);
ylabel(ax3, "dT/dt", 'FontSize', 20);
legend(ax3, "x_f = " + string(x_f_array), 'FontSize', 14, 'Location','best');
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')
set(ax3,"FontSize", 20); set(ax1, "FontSize", 20); set(ax2, "FontSize", 20);

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
    plot(ax1, xsol.x, xsol.y(1, :), 'LineWidth', 2);
    plot(ax2, xsol.x, xsol.y(2, :) - 273.15, 'LineWidth', 2);

    for j = 1:length(xsol.y(2, :))
        dTdt(:, j) = ODE(0, [xsol.y(1, j) xsol.y(2, j)], p);
    end
    semilogy(ax3, xsol.y(2, :) - 273.15, dTdt(2, :), 'LineWidth',4); clear dTdt;
end

title(ax1, "Dynamics when varying T0");
title(ax2, "Dynamics when varying T0");
xlabel(ax1, "Time (min)", 'FontSize', 20);
ylabel(ax1, "x_f", 'FontSize', 20);
xlabel(ax2, "Time (min)", 'FontSize', 20);
ylabel(ax2, "Temperature (C)", 'FontSize', 20);
legend(ax1, "T0 = " + string(T0_array - 273.15), 'FontSize', 14, 'Location','best');
legend(ax2, "T0 = " + string(T0_array - 273.15), 'FontSize', 14, 'Location','best');
title(ax3, "Dynamics when varying T0", 'FontSize', 20);
xlabel(ax3, "Temperature (C)", 'FontSize', 20);
ylabel(ax3, "dT/dt", 'FontSize', 20);
legend(ax3, "T0 = " + string(T0_array - 273.15), 'FontSize', 14, 'Location','best');
ylim(ax3, [10^-4 10^3])
set(ax3, 'YScale', 'log')
set(ax3,"FontSize", 20); set(ax1, "FontSize", 20); set(ax2, "FontSize", 20);

%Differential Equation Function
function f = ODE(~, x, p)
    x1 = x(1);
    T1 = x(2);
    A_1 = p(1);
    E_1 = p(2);
    k_B = p(3);
    n = p(4);
    K = p(5);

    if x1 < 0
        x1 = 0;
    end

    f = [-exp(-E_1/(k_B*T1))*A_1*(x1^n); K*exp(-E_1/(k_B*T1))*A_1*(x1^n)];

end
    