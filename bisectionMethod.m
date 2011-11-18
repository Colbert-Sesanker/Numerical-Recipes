%Bisection Method
function [k er pn c]=bi(i,it,tol) %enter i as a row vector
%Colbert Sesanker
a=i(1);
b=i(2);
flagg=0;
if ((f(a))||(f(b)))==0
    'the root is on the interval'
  flagg=1;
end
k=1;
while (k<it)&&(flagg==0)
    pn=(a+b)/2;
    if f(pn)==0
        'root is pn'
       flagg=1;
    elseif sign(f(a))*sign(f(pn))<0
       b=pn;
    else
        a=pn;
    end
    er=abs(a-b);%/abs(f(a));
    if er<tol
        flagg=1;
    end
    k=k+1;
end
c=f(pn);        
    
    
    
    function F=f(x)
     F=(x+2)*((x+1)^2)*x*((x-1)^3)*(x-2); %Problem 10
%     F=3*(x+1)*(x-(1/2))*(x-1); %Problem 2
%F=(x^3)+x-4; %Problem 14