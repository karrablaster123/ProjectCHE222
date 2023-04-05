clear all;

n = 0.5;
K = 150; % degrees Celsius
E_1 = 1.4; % eV
A_1 = 1 * 10^17;
k_B = 8.617 * 10^-5;
T0 = 80 + 273;
x_f = 0.1;
generate_array = @(CtrPt) linspace(0.1*CtrPt, 5*CtrPt, 50);

%Varying A1
A_1_array = generate_array(A_1);
for i = 1:length(A_1_array)
    p = [A_1_array(i), E_1, k_B, n, K];
    y1(i) = dTdt(T0, x_f, p);
end
plot(A_1_array, y1);
title("Initial dT/dt versus A");
xlabel("A");
ylabel("dT/dt");

%Varying K
clear y1 sol1;
K_array = generate_array(K);
for i = 1:length(K_array)
    p = [A_1, E_1, k_B, n, K_array(i)];
    y1(i) = dTdt(T0, x_f, p);
end
figure
semilogy(K_array, y1);
title("Initial dT/dt versus K");
xlabel("K");
ylabel("dT/dt");

%Varying E_1
clear y1 sol1;
E_1_array = linspace(1.38, E_1*5, 50);
for i = 1:length(E_1_array)
    p = [A_1, E_1_array(i), k_B, n, K];
    y1(i) = dTdt(T0, x_f, p);
end
figure
semilogy(E_1_array, y1);
title("Initial dT/dt versus E_1");
xlabel("E_1");
ylabel("dT/dt");

%Varying x_f
clear y1 sol1;
x_f_array = generate_array(x_f);
for i = 1:length(x_f_array)
    p = [A_1, E_1, k_B, n, K];
    y1(i) = dTdt(T0, x_f_array(i), p);
end
figure
semilogy(x_f_array, y1);
title("Initial dT/dt versus x_f");
xlabel("x_f");
ylabel("dT/dt");

%Varying n
clear y1 sol1;
n_array = generate_array(n);
for i = 1:length(n_array)
    p = [A_1, E_1, k_B, n_array(i), K];
    y1(i) = dTdt(T0, x_f, p);
end
figure
semilogy(n_array, y1);
title("Initial dT/dt versus n");
xlabel("n");
ylabel("dT/dt");

%Varying T0
clear y1 sol1;
T0_array = linspace(10+273, 150+273, 50);
for i = 1:length(T0_array)
    p = [A_1, E_1, k_B, n, K];
    y1(i) = dTdt(T0_array(i), x_f, p);
end
figure
semilogy(T0_array, y1);
title("Initial dT/dt versus T0");
xlabel("T0");
ylabel("dT/dt");

function y = dTdt(T1, x1, p)
    y = p(5)*exp(-p(2)/(p(3)*T1))*p(1)*(x1^p(4));
end
