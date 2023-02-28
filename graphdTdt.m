n = 0.5;
sol = solveODE(n, 160, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
K = 150; % degrees Celsius
E_1 = 1.4; % eV
A_1 = 1.25 * 10^17;
k_B = 8.617 * 10^-5;

for i = 1:length(T1)
    y(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n);
end
T1 = T1 - 273;
figure;
plot(T1, y); %This looks fine, the log one does not.
figure;
plot(T1, log10(y));

function y =  dTdt(K, E_1, k_B, T1, A_1, x1, n)
    y = K*exp(-E_1/(k_B*T1))*A_1*(x1^n);
end