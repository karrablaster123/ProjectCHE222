function wsol = solveODE(n, stoptime, T0)


K = 150; % degrees Celsius
x_f = 0.1;
E_1 = 1.4; % eV
A_1 = 1 * 10^17;
k_B = 8.617 * 10^-5;
T0 = T0 + 273.15;

abserr = 1.0e-10;
relerr = 1.0e-8;
numpoints = 250;
t = linspace(0, stoptime, numpoints);

p = [A_1, E_1, k_B, n, K];
w0 = [x_f, T0];

options = odeset('RelTol',relerr,'AbsTol',abserr, 'NonNegative', 1);

wsol = ode45(@(t,w)vectorfield(w,t,p), t, w0, options);
%plot(wsol.x, wsol.y(1, :));
%figure;
%plot(wsol.x, wsol.y(2, :));
    
    function f = vectorfield(w, t, p) %#ok<INUSL> 
        x1 = w(1);
        T1 = w(2);
        A_1 = p(1);
        E_1 = p(2);
        k_B = p(3);
        n = p(4);
        K = p(5);
        
        if x1<0
            x1 = 0;
        end

        f = [-exp(-E_1/(k_B*T1))*A_1*(x1^n); K*exp(-E_1/(k_B*T1))*A_1*(x1^n)];
    end

end