% The 'p' argument is necessary if you want 
% To feed ODE45 a parameter argument
% not necessary for fixed parameters

function dy = predator_prey_ode45(t, y, p)

x  = y(1);
y  = y(2);


% Set up parameters
a = p(1);
b = p(2);
c = p(3);
d = p(4);

% Evaluate equations
dy    = zeros(2,1);    % a column vector

dy(1) = a*x   - b*x*y;
dy(2) = c*x*y - d*y;
end
