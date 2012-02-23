%% Colbert Sesanker 4/1/11
%% Multistep Numerical Differential Equation Compared with exact in Tabular Form 
%% Adams-Moulton methods using 2,3 and 4 steps are used to solve the follwoing differential equations: y'= 1+y/t             (1)
%%                                                                                                     y'=t*exp(3*t) -2*y    (2)
%% 											               y'= cos(2t) +sin(3t)  (3)
%% The difference equation is explicitly solved and and simplified in Maple 
h=.02;
N=floor(1/h);

x=0:h:1;   
x2=1:h:2;
y1=x2.*log(x2)+2.*x2; % exact solution EQN (1)
y2=(1/5).*x.*exp(3.*x)-(1/25)*exp(3.*x)+(1/25)*exp(-2.*x); %exact solution EQN (2)
y3=(1/2)*sin(2.*x)-(1/3)*cos(3.*x)+(4/3);  %exact solution EQN (3)  


for i=1:4             %4 initial Conditions           
    w1(i)=y1(i); %Solve EQN (c) on [1,2]
    v1(i)=y2(i);  %Solve EQN (a) on [0,1]
    r1(i)=y3(i);  %Solve EQN (d) on [0,1]
end                        

%Initilize conditions for the 3 Adams-Moulton Methods
w2=w1; w3=w1;  v2=v1; v3=v1;  r2=r1; r3=r1; 

    
   
%%%%%%%%%%%%%%%%%%%%2-STEP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    for i=1:N+1    
        W1(i)=(w1(2)+...
        (h/12)*(12+8*w1(2)/(1+h*i)-w1(1)/(1+h*(i-1))))/(1-(5*h/(12*(1+h*(i+1))))); %%Eqn 1
    
        w1(1)=w1(2);
        w1(3)=W1(i);
        
        V1(i)=(v1(2)+(1/12)*h*(5*h*(i+1)*exp(3*h*(i+1)) + ...
            8*h*i*exp(3*h*i)-16*v1(2)-h*(i-1)*exp(3*h*(i-1))+2*v1(1)))/ (1+(10*h/12)); %%Eqn 2
        
        v1(1)=v1(2);
        v1(3)=V1(i);
        
        R1(i)=r1(2)+(1/12)*h*(5*cos(2*h*(i+1))+5*sin(3*h*(i+1))+...
            8*cos(2*h*i)+8*sin(3*h*i)-cos(2*h*(i-1))-sin(3*h*(i-1))); %%Eqn 3
        
        r1(1)=r1(2);
        r1(3)=R1(i);
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3-STEP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    for i=2:N+2    
            W2(i)=(w2(3)+(h/24)*(24+((19*w2(3))/(1+h*i))-(5*w2(2)/(1+h*(i-1)))+... %%Eqn 1
            (w2(1)/(1+h*(i-2)))))/ (1-(9*h/(24*(1+h*(i+1)))));
        
        w2(1)=w2(2);
        w2(2)=w2(3);
        w2(3)=W2(i);
        
        V2(i)=(v2(3)+(1/24)*h*(9*h*(i+1)*exp(3*h*(i+1)) +...
            19*h*i*exp(3*h*i)-38*v2(3)-5*h*(i-1)*exp(3*h*(i-1))+...  %%Eqn 2
            10*v2(2)+h*(i-2)*exp(3*h*(i-2))-2*v2(2)))/ (1+(18*h/24));

        v2(1)=v2(2);
        v2(2)=v2(3);
        v2(3)=V2(i);
        
        R2(i)=r2(2)+(1/24)*h*(9*cos(2*h*(i+1))+9*sin(3*h*(i+1))+...  %%Eqn 3
            19*cos(2*h*i)+19*sin(3*h*i)-5*cos(2*h*(i-1))-...
            5*sin(3*h*(i-1))+cos(2*h*(i-2))+sin(3*h*(i-2)));
        
        r2(1)=r2(2);
        r2(2)=r2(3);
        r2(3)=R2(i);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4-STEP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
       for i=3:N+3
            W3(i)=(w3(4)+(h/720)*(720+(646*w3(4)/(1+h*i))-(264*w3(3)/(1+h*(i-1)))+...
            (106*w3(2)/(1+h*(i-2)))-(19*w3(1)/(1+h*(i-3)))))/(1-251*h/(720*(1+h*(i+1)))); %%Eqn 1
        
        for l=1:3
        w3(l)=w3(l+1);
        w3(4)=W3(i);
        end   
        
        V3(i)=(v3(4)+(1/720)*h*(251*h*(i+1)*exp(3*h*(i+1)) +...
            646*h*i*exp(3*h*i)-1292*v3(3)-264*h*(i-1)*exp(3*h*(i-1))+...  %%Eqn 2
            528*v3(3)+106*h*(i-2)*exp(3*h*(i-2))-212*v3(3)-...
            19*h*(i-3)*exp(3*h*(i-3))+38*v3(1)))/ (1+(502*h/720));
        
        for l=1:3
        v3(l)=v3(l+1);
        v3(4)=V3(i);
        end      
        
        R3(i)=r3(4)+(1/720)*h*(251*cos(2*h*(i+1))+251*sin(3*h*(i+1))+...
            646*cos(2*h*i)+646*sin(3*h*i)-264*cos(2*h*(i-1))-...
            264*sin(3*h*(i-1))+106*cos(2*h*(i-2))+106*sin(3*h*(i-2))-...  %%Eqn 3
            19*cos(2*h*(i-3))-19*sin(3*h*(i-3)));
        
        for l=1:3
        r3(l)=r3(l+1);
        r3(4)=R3(i);
        end   
        
       end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       y1N=y1(N);
       y2N=y2(N);
       y3N=y3(N);
    
%%%%%%%%%%%%%%%%%%%%Generate Error Tables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    for j=1:10
        %Errors for EQN (1) |e|/(h^j) index across jth coulmn
        ew(1,j)=abs(W1(N)-y1N)/(h^(j-1));
        ew(2,j)=abs(W2(N)-y1N)/(h^(j-1));
        ew(3,j)=abs(W3(N)-y1N)/(h^(j-1));
       
                       
       %Errors for EQN (2) |e|/(h^j) index across jth coulmn
        ev(1,j)=abs(V1(N)-y2N)/(h^(j-1));
        ev(2,j)=abs(V2(N)-y2N)/(h^(j-1));
        ev(3,j)=abs(V3(N)-y2N)/(h^(j-1));
             
      %Errors for EQN (3) |e|/(h^j) index across jth coulmn
        er(1,j)=abs(R1(N)-y3N)/(h^(j-1));
        er(2,j)=abs(R2(N)-y3N)/(h^(j-1));
        er(3,j)=abs(R3(N)-y3N)/(h^(j-1));
   
                       
        
    end

%%%Simply Print the matrix to view the table%%%%

subplot(1,3,1)
 plot(R1,'r');
 hold on
 plot(R2,'b');
 hold on
 plot(R3,'g');
 hold on  
 plot(y3,'k');
 hold on

 
subplot(1,3,2)
 plot(V1,'r');
 hold on
 plot(V2,'b');
 hold on
 plot(V3,'g');
 hold on
 plot(y2,'k') 
 hold on

 
subplot(1,3,3)
 plot(W1,'r');
 hold on
 plot(W2,'b');
 hold on
 plot(W3,'g');
 hold on
 plot(y1,'k'); 

 
 
 
