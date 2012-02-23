%%Author: Colbert Sesanker 
%%2/1/2011
function [k er] =it(n,alpha,IT,TOL,nm,method,omega)

%%% IT is the total number of iterations before you stop, TOL is the
%%% the upper bound for error stopping criteria. nm is the choice of norm for measuring error
%%% method 1,2 and any other number represents Jacobi, Gauss Sieldel and
%%% SOR methods respectively. See procedure for automatic table generation

%%%Intitialize variables and Matricies 
x=zeros(2,n); %X0 is initialized as x(1,:) and 
b=zeros(1,n); % 
b(1)=b(1)+alpha; 
A=zeros(n,n);  %Tridiagonal n*n Matrix
k=1; %Initialize counter to 1

%%%Method Choosing Scheme: method 1 is Jacobi, method 2 is Gauss Siedel and
%%%method 3 is SOR.
if method==1
    g=1;
    omega=1;
elseif method==2
    g=2;
    omega=1;
else g=2;
end

%%%Create the tridiagonal Matrix A%%%
for i=1:n-1
A(i,i)=A(i,i)+1;
A(i+1,i)=A(i+1,i)-alpha;
A(i,i+1)=A(i,i+1)-1+alpha;
end
A(n,n)=A(n,n)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flagg=0;
%%%Algorythym for Jacobi, Gauss Siedel or SOR depending on request%%%
er=TOL+1; 
while (k<IT)&&(flagg==0)
    for i=1:n
        temp=0;
         if i~=1
        for j=1:i-1
            temp=temp+A(i,j)*x(g,j);
        end
        end
        for j=i+1:n
            temp=temp+A(i,j)*x(1,j);
        end
        x(2,i)=(1-omega)*x(1,i)+(omega/A(i,i))*(b(i)-temp);
    end
    if k==1
        x1=x(2,:);
    end
    er=norm((x(2,:)-x(1,:)),nm)/norm(x1,nm);
    if (er<TOL)
        flagg=1;
    elseif (isnan(er))||(isinf(er))
        k='div';
    end
    x(1,:)=x(2,:);
    k=k+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
