function wsol = solveODE()
import matlab.*
import matlab.math.*
import matlab.graphics.*
import matlab.io.*

K = 150; % degrees Celsius
x_f = 0.1;
E_1 = 1.4; % eV
A_1 = 1.25 * 10^17;
k_B = 8.617 * 10^-5;
T0 = 80+273;
n = 0.5;

abserr = 1.0e-8;
relerr = 1.0e-6;
stoptime = 170;
numpoints = 250;
t = linspace(0, stoptime, numpoints);

p = [A_1, E_1, k_B, n, K];
w0 = [x_f, T0];

options = odeset('RelTol',relerr,'AbsTol',abserr);

wsol = ode45(@(t,w)vectorfield(w,t,p), t, w0, options);
plot(wsol.x, wsol.y(1, :));
figure;
plot(wsol.x, wsol.y(2, :));
    
    function f = vectorfield(w, t, p)
        x1 = w(1);
        T1 = w(2);
        A_1 = p(1);
        E_1 = p(2);
        k_B = p(3);
        n = p(4);
        K = p(5);

        f = [-exp(-E_1/(k_B*T1))*A_1*(x1^n); K*exp(-E_1/(k_B*T1))*A_1*(x1^n)];
    end

end