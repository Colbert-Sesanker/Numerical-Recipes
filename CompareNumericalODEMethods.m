%Method Comparison Function on [0,1]
% Colbert Sesanker 4/2/11 
% This script generates 10 tables. 5 tables for each of two ODEs using
% Five approximation methods, involving explicit, implicit and multistep methods. 
% We are solving following linear ODEs: y' = 1000(y-t^2) + 2t and y' = -1000(y-t^2) + 2t
% with  exact solution y = t^2 on [0,1]. The sign change on 1000 causes some methods to blow up. 
% Runs through all the methods and plots the error tables
% Method with error term e/h^n is or order n when e/h^n  is close to 1 (see tables)

format shortG; %Formatting for tables
format compact;

% Initial Conditions:
v4(1)=0; w4(1)=0; v03=0; v02=0; v0=0; w03=0; w02=0; w0=0;
k=1;
while k<=8
    N=10^k;
    h=1/N;
    for i=1:3
        w4(i+1)=(h*i)^2; % The exact solutions
        v4(i+1)=(h*i)^2;
    end
    w5=w4;
    v5=v4;
  
    for i=0:N-1
       % Euler methods for equation with +1000 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
        W1=w0+h*(1000*(w0-(i*h)^2))+2*i*(h^2); % Forward Euler
        w0=W1;
        
        W2=(w02-(h^3)*(1000*(i+1)^2)+2*(i+1)*(h^2))/(1-1000*h); % Backward Euler
        w02=W2;
        
        W3=(w03+(h/2)*(i*h*(2-1000*i*h)+(i+1)*h*(2-1000*(i+1)*h)))...
            /(1-1000*h); % Trapezoidal
        w03=W3;
        
        % Euler methods for equation with -1000 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        V1=v0-h*(1000*(v0-(i*h)^2))+2*(i*h^2); % Forward Euler
        v0=V1;
        
        V2=(v02+(h^3)*(1000*(i+1)^2)+2*(i+1)*(h^2))/(1+1000*h); % Backward Euler
        w02=W2;
        v02=V2;
        
        V3=(v03+(h/2)*(i*h*(2+1000*i*h)+(i+1)*h*(2+1000*(i+1)*h)))...
            /(1+1000*h); % Trapezoidal
        v03=V3;
    end
    
    % Four Step Bashford and Moulton Methods
    for i=3:N+3
        %Adams-Moulton Four-Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W4=(w4(4)+(1/720)*h*(-19000*w4(1)+1440*h*i-720000*(h^2)*(i^2)+...
        646000*w4(4)+ 720*h-240000*(h^2)-720000*(h^2)*i-264000*w4(3)+...
        106000*w4(2)))/(1-(251000*h/720));
        for l=1:3
        w4(l)=w4(l+1);
        w4(4)=W4;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Adams Bashford Four-Step
        W5= w5(4)+(1/24)*h*(55000*w5(4)-24000*(h^2)*i^2+48*h*i-59000*w5(3)...
            -24000*(h^2)*i-8000*(h^2)+24*h+37000*w5(2)-9000*w5(1));
        for l=1:3
        w5(l)=w5(l+1);
        w5(4)=W5;
        end        
       
        %Adams-Moulton Four-Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        V4=(v4(4)+(1/720)*h*(1440*h*i+720000*(h^2)*(i^2)+264000*v4(3)+720*h+...
        240000*(h^2)+ 720000*(h^2)*i+19000*v4(1)-646000*v4(4)-106000*v4(2)))...
        /(1+(251000*h/720));
    
        for l=1:3
        v4(l)=v4(l+1);
        v4(4)=V4;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Adams Bashford Four-Step
        V5= v5(4)+(1/24)*h*(-55000*v5(4)+24000*(h^2)*(i^2)+48*h*i+...
            59000*v5(3)+24000*(h^2)*i+8000*(h^2)+24*h-37000*v5(2)+...
            9000*v5(1));
        for l=1:3
        v5(l)=v5(l+1);
        v5(4)=V5;
        end
    end
    
    % Create Tables For Global Error
    for j=1:5
        
         format shortG; % Formating for display in interpreter
         format compact;
        %Error tables for 5 methods using the equation with coeficient +10000
        ew1(k,j)=abs(W1-1)/(h^(j-1));
        ew2(k,j)=abs(W2-1)/(h^(j-1));
        ew3(k,j)=abs(W3-1)/(h^(j-1));
        ew4(k,j)=abs(W4-1)/(h^(j-1));
        ew5(k,j)=abs(W5-1)/(h^(j-1));
        
        
        %Error tables for 5 methods using the equation with coeficient -10000
        ev1(k,j)=abs(V1-1)/(h^(j-1));
        ev2(k,j)=abs(V2-1)/(h^(j-1));
        ev3(k,j)=abs(V3-1)/(h^(j-1));
        ev4(k,j)=abs(V4-1)/(h^(j-1));
        ev5(k,j)=abs(V5-1)/(h^(j-1));     
             
    end
    k=k+1;
end
    
        ew1  
        ew2  
        ew3
        ew4
        ew5                               
          
        ev1  
        ev2  
        ev3
        ev4
        ev5   
        
%{ 
ew1 =   e         e/h          e/h^2        e/h^3        e/h^4      
        
 N   9.1485e+017  9.1485e+018  9.1485e+019  9.1485e+020  9.1485e+021
 N^2 1.2607e+122  1.2607e+124  1.2607e+126  1.2607e+128  1.2607e+130
 N^3      Inf          Inf          Inf          Inf          Inf
 N^4      Inf          Inf          Inf          Inf          Inf
 .        Inf          Inf          Inf          Inf          Inf
 .        Inf          Inf          Inf          Inf          Inf
 .        Inf          Inf          Inf          Inf          Inf
 N^8      Inf          Inf          Inf          Inf          Inf

ew2 =
       0.0001        0.001         0.01          0.1            1
       1e-005        0.001          0.1           10         1000
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
ew3 =
       0.0951        0.951         9.51         95.1          951
      0.00996        0.996         99.6         9960    9.96e+005
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
          Inf          Inf          Inf          Inf          Inf
ew4 =
      0.68766       6.8766       68.766       687.66       6876.6
  2.8216e+034  2.8216e+036  2.8216e+038  2.8216e+040  2.8216e+042
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
ew5 =
  1.9434e+013  1.9434e+014  1.9434e+015  1.9434e+016  1.9434e+017
  3.5517e+064  3.5517e+066  3.5517e+068  3.5517e+070  3.5517e+072
  7.4548e+246  7.4548e+249  7.4548e+252  7.4548e+255  7.4548e+258
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
          NaN          NaN          NaN          NaN          NaN
ev1 =
  1.0672e+016  1.0672e+017  1.0672e+018  1.0672e+019  1.0672e+020
  2.8346e+111  2.8346e+113  2.8346e+115  2.8346e+117  2.8346e+119
        1.996         1996   1.996e+006   1.996e+009   1.996e+012
       19.978  1.9978e+005  1.9978e+009  1.9978e+013  1.9978e+017
        199.8   1.998e+007   1.998e+012   1.998e+017   1.998e+022
         1998   1.998e+009   1.998e+015   1.998e+021   1.998e+027
        19980   1.998e+011   1.998e+018   1.998e+025   1.998e+032
   1.998e+005   1.998e+013   1.998e+021   1.998e+029   1.998e+037
ev2 =
       0.0001        0.001         0.01          0.1            1
       1e-005        0.001          0.1           10         1000
       1e-006        0.001            1         1000       1e+006
       1e-007        0.001           10       1e+005       1e+009
       1e-008        0.001          100       1e+007       1e+012
  1.0001e-009    0.0010001       1000.1  1.0001e+009  1.0001e+015
  1.0011e-010    0.0010011        10011  1.0011e+011  1.0011e+018
  3.4501e-012   0.00034501        34501  3.4501e+012  3.4501e+020
ev3 =
       0.0949        0.949         9.49         94.9          949
      0.00994        0.994         99.4         9940    9.94e+005
    0.0009985       0.9985        998.5   9.985e+005   9.985e+008
  9.9895e-005      0.99895       9989.5  9.9895e+007  9.9895e+011
    9.99e-006        0.999        99900    9.99e+009    9.99e+014
    9.99e-007        0.999    9.99e+005    9.99e+011    9.99e+017
    9.99e-008        0.999    9.99e+006    9.99e+013    9.99e+020
  9.9965e-009      0.99965  9.9965e+007  9.9965e+015  9.9965e+023
ev4 =
      0.96624       9.6624       96.624       966.24       9662.4
      0.08669        8.669        866.9        86690   8.669e+006
     0.008509        8.509         8509   8.509e+006   8.509e+009
    0.0008493        8.493        84930   8.493e+008   8.493e+012
  8.4914e-005       8.4914  8.4914e+005  8.4914e+010  8.4914e+015
  8.4912e-006       8.4912  8.4912e+006  8.4912e+012  8.4912e+018
  8.4912e-007       8.4912  8.4912e+007  8.4912e+014  8.4912e+021
  8.4911e-008       8.4911  8.4911e+008  8.4911e+016  8.4911e+024
ev5 =
  1.5799e+011  1.5799e+012  1.5799e+013  1.5799e+014  1.5799e+015
  1.2341e+059  1.2341e+061  1.2341e+063  1.2341e+065  1.2341e+067
  1.5092e+054  1.5092e+057  1.5092e+060  1.5092e+063  1.5092e+066
    0.0010584       10.584  1.0584e+005  1.0584e+009  1.0584e+013
   0.00010581       10.581  1.0581e+006  1.0581e+011  1.0581e+016
  1.0581e-005       10.581  1.0581e+007  1.0581e+013  1.0581e+019
  1.0581e-006       10.581  1.0581e+008  1.0581e+015  1.0581e+022
  1.0581e-007       10.581  1.0581e+009  1.0581e+017  1.0581e+025

  %}

    


