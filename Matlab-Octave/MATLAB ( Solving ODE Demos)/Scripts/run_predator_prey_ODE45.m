
%%% ODE45 with parameter Arguments tutorial %%%%%
%%% Colbert Sesanker, 5/26/2015, MATLAB R2015a

%%%%%%%%%%%% Model 1: Predator Prey %%%%%%%%%%%%%


% Add all files from folders and subfolders on path
% comes in handy when need to access lots of files 
% if you have lots of directories
addpath(genpath('./'))

% Close all Windows
close all;

% Initial Conditions
x =  3;
y =  3; 

initialValues = [x y]; 
equations     = @predator_prey_ode45;


% Parameters all set to 1
a = 1;
b = 1; 
c = 1; 
d = 1;
params = [a b c d];

% Set tolerances
odeOptions    = odeset('RelTol',1e-6,'AbsTol',1e-6);

%%%%%%%%% Basic ODE45 use %%%%%%%%%%%%%%%%%%%%%%%%
% Note: This is usefull if you are trying to fit data
% Hence the returned 'timeData' is exactly the points you put in
% Note you need odeOptions to use params

startTime     = 0;
endTime       = 100;
timeInterval  = [startTime endTime];

[timeData, stateEstimates] = ode45( @predator_prey_ode45,...
                                    timeInterval,...
                                    initialValues,... 
                                    odeOptions,...
                                    params);

% Note speciesEstimates is a solution matrix with
% Dimensions TxN where T is the number of time points and
% N is the number of state variables

% Extract out x and y
x_solved = stateEstimates(:, 1);  
y_solved = stateEstimates(:, 2);

%%%%%% Plot Solutions %%%%%%%%
plot(timeData, x_solved, 'r');
% This command is extremely usefull 
% Used to plot multiple things on same plot
hold on;
plot(timeData, y_solved, 'b');

% Plot in State Space
% Figure(i) initializes a new figure i
figure(2);
plot(x_solved, y_solved);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Specify Time Points Manually %%%%%%%%
% Note: This is usefull if you are trying to fit data
numTimePts    = 3000;
startTime     = 0;
endTime       = 500;
step          = (endTime - startTime) / numTimePts;
timePoints    = startTime: step :endTime;

[~, stateEstimates] = ode45( equations,...
                             timePoints,...
                             initialValues,... 
                             odeOptions,...
                             params);
% Plot in State Space
% Figure(i) initializes a new figure i
x_solved = stateEstimates(:, 1);  
y_solved = stateEstimates(:, 2);
figure(3);
plot(x_solved, y_solved);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%% Trade offs with Integration Tolerances and computation Time %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Abs and Rel Tolerances    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% RelTol: %%%%%%%%
% "This tolerance is a measure of the error relative to the size of each solution component. 
%  Roughly, it controls the number of correct digits in all solution components, except those smaller than thresholds AbsTol(i). 
%  The default, 1e-3, corresponds to 0.1% accuracy. " Mathworks Documentation

%%%% rel-Tol = abs(X - X*) / min(abs(X), abs(X*)), X is the numercal
%%%% estimate and X* is the exact solution

%%%%%% AbsTol %%%%%%
% "AbsTol(i) is a threshold below which the value of the ith solution component is unimportant. 
% The absolute error tolerances determine the accuracy when the solution
% approaches zero." Mathworks Documentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Here we iterate over different rel-tols and see how long
%%% ODE-45 takes
numTols = 10;
times = zeros(1, numTols);
startTime     = 0;
endTime       = 500;
timeInterval  = [startTime endTime];

for i = 1:numTols
relTol        = .1^(i);
% The  tic and toc pair allows you to time anything in between
tic;
[timeData, stateEstimates] = ode45( equations,...
                                    timeInterval,...
                                    initialValues,...      
                                    odeOptions,...                             
                                    params);
                           
times(i) = toc;
disp(['rel-tol and corresponding time: ' ...
       num2str(relTol) '    ' ...
       num2str(times(i))]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Iterative Solving: Predator Prey Solutions for differents ICs %%%
curves = 5;
x      =  1;
y      =  1; 
initialValues = [x y]; 

% Pick up plotting where last figure left off
figure(4);
for i = 1:curves

initialValues = initialValues + i;
startTime     = 0;
endTime       = 500;
[timeData, stateEstimates] = ode45( equations,...
                                    timeInterval,...
                                    initialValues,... 
                                    odeOptions,...
                                    params);                    
x_solved = stateEstimates(:, 1);  
y_solved = stateEstimates(:, 2);
plot(x_solved, y_solved);
hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%%%% Iterative Solving: Vary a parameter %%%%
%%%%       Fitz-Hugh Nagumo Model        %%%%
grid   =  50;
V      =  -1;
R      =   1; 
initialValues = [V R]; 
params        = [.2 .2 3];
startTime     = 0;
endTime       = 200;
timeInterval  = [startTime endTime];
step   =  .2;
equations = @fitzHughNagumo;
% Pick up plotting where last figure left off
figure(5);
for i  = 1:grid
params(2)     = params(2) + step;
[timeData, stateEstimates] = ode45(equations,...
                                   timeInterval,...
                                   initialValues,...  
                                   odeOptions,...
                                   params);                    
V_solved = stateEstimates(:, 1);  
R_solved = stateEstimates(:, 2);
plot(V_solved, R_solved);
hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
