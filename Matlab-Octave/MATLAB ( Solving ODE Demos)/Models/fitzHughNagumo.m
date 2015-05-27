function dy = fitzHughNagumo(t, y, p)
 
% Use V = -1, R = 1 for initial conditions
% Set up states
V  = y(1);
R  = y(2);


% Parameters: use a  = 0.2, b = 0.2, c = 3
% Use time interval of [0 200]
a = p(1);
b = p(2);
c = p(3);


dy    = zeros(2,1);    % a column vector

dy(1) = c*(V - (V^3) / 3 + R);
dy(2) = - (V - a + b * R) / c;

end
