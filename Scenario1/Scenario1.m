clear all;
%{
Scenario 1 - Group 21
Kush Patel - 1007678826
Sachiel Malik - 1007692098
Gurjas Singh Chawla - 1007872189
Uddhava Swaminathan - 1007511122

This code aims to predict the reaction order by matching different reaction
orders to the empirical data. 
%}

%Parameter Definition
n = [0.2, 0.333, 0.5, 1];
K = 150;
E_1 = 1.4;
A_1 = 1.25 * 10^17;
k_B = 8.617 * 10^-5;

%Get Solution
sol = solveODE(n(1), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
%For each solution point, get the dT/dt at the point
for i = 1:length(T1)
    y1(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(1)); %#ok<*SAGROW> 
end
T1 = T1 - 273.15;
%plot
semilogy(T1, y1, '-r', 'DisplayName', 'n = 0.2');
hold on;

%Repeat the same steps for different n values and different temperatures
%input of the form solveODE(n_val, max_time, temp0)
sol = solveODE(n(2), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y2(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(2));
end
T1 = T1 - 273.15;
semilogy(T1, y2, '--r', 'DisplayName', 'n = 0.333');

sol = solveODE(n(3), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y3(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(3));
end
T1 = T1 - 273.15;
semilogy(T1, y3, ':r', 'DisplayName', 'n = 0.5');

sol = solveODE(n(4), 5000, 80);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y4(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(4));
end
T1 = T1 - 273.15;
semilogy(T1, y4, '-.r', 'DisplayName', 'n = 1');


sol = solveODE(n(1), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y5(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(1));
end
T1 = T1 - 273.15;
semilogy(T1, y5, '-m', 'DisplayName', 'n = 0.2');


sol = solveODE(n(2), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y6(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(2));
end
T1 = T1 - 273.15;
semilogy(T1, y6, '--m', 'DisplayName', 'n = 0.333');

sol = solveODE(n(3), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y7(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(3));
end
T1 = T1 - 273.15;
semilogy(T1, y7, ':m', 'DisplayName', 'n = 0.5');

sol = solveODE(n(4), 5000, 100);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y8(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(4));
end
T1 = T1 - 273.15;
semilogy(T1, y8, '-.m', 'DisplayName', 'n = 1');

sol = solveODE(n(1), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y9(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(1));
end
T1 = T1 - 273.15;
semilogy(T1, y9, '-k', 'DisplayName', 'n = 0.2');


sol = solveODE(n(2), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y10(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(2));
end
T1 = T1 - 273.15;
semilogy(T1, y10, '--k', 'DisplayName', 'n = 0.333');

sol = solveODE(n(3), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y11(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(3));
end
T1 = T1 - 273.15;
semilogy(T1, y11, ':k', 'DisplayName', 'n = 0.5');

sol = solveODE(n(4), 5000, 130);
x1 = sol.y(1, :);
T1 = sol.y(2, :);
for i = 1:length(T1)
    y12(i) = dTdt(K, E_1, k_B, T1(i), A_1, x1(i), n(4));
end
T1 = T1 - 273.15;
semilogy(T1, y12, '-.k', 'DisplayName', 'n = 1');
%Format Plots
ylim([0.01 100]);
xlim([80 220]);
legendUnq();
legend show;
legend('Location', 'best');

xlabel("Temperature (°C)");
ylabel("dT/dt(°C/min) - Semilog Scale");
set(0,'DefaultLineLineWidth',4)

%Differential Function
function y =  dTdt(K, E_1, k_B, T1, A_1, x1, n)
    y = K*exp(-E_1/(k_B*T1))*A_1*(x1^n);
end

%ODE solver
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

%Ignore this code but don't delete it
function unqLegHands = legendUnq(h, sortType)
% unqLegHands = legendUnq(h, sortType)
%   Run this function just before running 'legend()' to avoid representing duplicate or missing
% DisplayNames within the legend. This solves the problem of having a cluttered legend with 
% duplicate or generic values such as "data1" assigned by matalab's legend() function.  This 
% also makes is incredibly easy to assign one legend to a figure with multiple subplots. 
% Use the 'DisplayName' property in your plots and input the axis handle or the figure handle 
% so this code can search for all potential legend elements, find duplicate DisplayName strings, 
% and remove redundent components by setting their IconDisplayStyle to 'off'.  Then call 
% legend(unqLegHands) to display unique legend components. 
% INTPUT
%       h: (optional) either a handle to a figure, an axis, or a vector of axis handles. The code 
%           will search for plot elements in all axes belonging to h.  If h is missing, gca is used.
%       sort: (optional) can be one of the following strings that will sort the unqLeghands.
%           'alpha': alphabetical order. 
% OUTPUT
%       unqLegHands: a list of handles that have unique DisplayNames; class 'matlab.graphics.chart.primitive.Line'.
%           ie: unqLegHands = legendUnq(figHandle); legend(unqLegHands)
% EXAMPLE 1: 
%         figure; axis; hold on
%         for i=1:10
%             plot(i,rand(), 'ko', 'DisplayName', 'randVal1');        % included in legend
%             plot(i+.33, rand(), 'ro', 'DisplayName', 'randVal2');   % included in legend       
%         end
%         plot(rand(1,10), 'b-'); 	% no DisplayName so it is absent from legend
%         legend(legendUnq())
% EXAMPLE 2: 
%         fh = figure; subplot(2,2,1); hold on
%         plot(1:10, rand(1,10), 'b-o', 'DisplayName', 'plot1 val1')
%         plot(1:2:10, rand(1,5), 'r-*', 'DisplayName', 'plot1 val2')
%         subplot(2,2,2); hold on
%         plot(1:10, rand(1,10), 'm-o', 'DisplayName', 'plot2 val1')
%         plot(1:2:10, rand(1,5), 'g-*', 'DisplayName', 'plot2 val2')
%         subplot(2,2,3); hold on
%         plot(1:10, rand(1,10), 'c-o', 'DisplayName', 'plot3 val1')
%         plot(1:2:10, rand(1,5), 'k-*', 'DisplayName', 'plot3 val2')
%         lh = legend(legendUnq(fh)); 
%         lh.Position = [.6 .2 .17 .21];
%
% Danz 180515

% Change history
% 180912 fixed error when plot is empty
% 180913 adapted use of undocumented function for matlab 2018b

persistent useOldMethod

% If handle isn't specified, choose current axes
if nargin == 0
    h = gca; 
end

% If user entered a figure handle, replace with a list of children axes; preserve order of axes
if strcmp(get(h, 'type'), 'figure')
    h = flipud(findall(h, 'type', 'Axes')); 
end

% set flag to use old method of obtaining legend children
% In 2018b matlab changed an undocumented function that obtains legend handles. 
useOldMethod = verLessThan('matlab', '9.5.0'); 

% Set the correct undocumented function 
if useOldMethod
    getLegendChildren = @(x) graph2dhelper('get_legendable_children', x);
else
    getLegendChildren = @(x) matlab.graphics.illustration.internal.getLegendableChildren(x);
end

% Get all objects that will be assigned to legend.
% This uses an undocumented function that the legend() func uses to get legend componenets.
legChildren = matlab.graphics.chart.primitive.Line; %initializing class (unsure of a better way)
for i = 1:length(h)
    temp = getLegendChildren(h(i));
    if ~isempty(temp)
        legChildren(end+1:end+length(temp),1) = temp; 
    end
end
legChildren(1) = [];
% Get display names
dispNames = get(legChildren, 'DisplayName');
if isempty(dispNames)
    dispNames = {''}; 
end
if ~iscell(dispNames)
    dispNames = {dispNames}; 
end
% Find the first occurance of each name 
[~, firstIdx] = unique(dispNames, 'first'); 
% Create an index of legend items that will be hidden from legend
legRmIdx = true(size(legChildren)); 
legRmIdx(firstIdx) = false; 
% Add any elements that have no displayName to removal index (legend would assign them 'dataX')
legRmIdx = legRmIdx | cellfun(@isempty,dispNames);
% get all annotations
annot = get(legChildren, 'Annotation'); 
% Loop through all items to be hidden and turn off IconDisplayStyle
for i = 1:length(annot)
    if legRmIdx(i)
        set(get(annot{i}, 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
end
% Output remaining, handles to unique legend entries
unqLegHands = legChildren(~legRmIdx); 

% Sort, if user requested
if nargin > 1 && ~isempty(sortType) && length(unqLegHands)>1
    [~, sortIdx] = sort(get(unqLegHands, 'DisplayName'));
    unqLegHands = unqLegHands(sortIdx); 
end

end


