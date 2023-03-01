n = [0.2, 0.333, 0.5, 1];


K = 150; % degrees Celsius
E_1 = 1.4; % eV
A_1 = 1.25 * 10^17;
k_B = 8.617 * 10^-5;

sol = solveODE(n(1), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y1(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(1));
end
T1 = T1 - 273;
plot(T1, y1);
hold on;

sol = solveODE(n(2), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y2(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(2));
end
T1 = T1 - 273;
plot(T1, y2);

sol = solveODE(n(3), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y3(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(3));
end
T1 = T1 - 273;
plot(T1, y3);

sol = solveODE(n(4), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y4(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(4));
end
T1 = T1 - 273;
plot(T1, y4);

legendString = "n = " + string(n);
legend(legendString);
title("T = 80");
figure;
plot(T1, y4);
legend("n = 1");

sol = solveODE(n(1), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y5(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(1));
end
T1 = T1 - 273;
figure;
hold on;
plot(T1, y5);


sol = solveODE(n(2), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y6(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(2));
end
T1 = T1 - 273;
plot(T1, y6);

sol = solveODE(n(3), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y7(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(3));
end
T1 = T1 - 273;
plot(T1, y7);

sol = solveODE(n(4), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y8(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(4));
end
T1 = T1 - 273;
plot(T1, y8);

legendString = "n = " + string(n);
legend(legendString);
title("T = 100");
figure;
plot(T1, y8);
legend("n = 1");

sol = solveODE(n(1), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y9(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(1));
end
T1 = T1 - 273;
figure;
hold on;
plot(T1, y9);


sol = solveODE(n(2), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y10(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(2));
end
T1 = T1 - 273;
plot(T1, y10);

sol = solveODE(n(3), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y11(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(3));
end
T1 = T1 - 273;
plot(T1, y11);

sol = solveODE(n(4), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y12(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(4));
end
T1 = T1 - 273;
plot(T1, y12);

legendString = "n = " + string(n);
legend(legendString);
title("T = 130");
figure;
plot(T1, y12);
legend("n = 1");

function y =  dTdt(K, E_1, k_B, T1, A_1, x1, n)
    y = K*exp(-E_1/(k_B*T1))*A_1*(x1^n);
end