%%%%%%Dynamic Inverse Power Method%%%%%%%
%%Colbert Sesanker 3/2/11

function [lM]= IP(x,A,it,tol) %enter x as coulm vector, A is the matrix, it the iterations and tol the tolerance

q=(x'*(A*x))./(x'*x);

p=find(abs(x)==norm(x,inf),1); %Find p; the 1 gets the least integer
x=x./x(p);

k=1;
n=size(A);
flagg=0;

while (k<it)&&(flagg==0)
    if det((A-(q.*eye(n,n))))==0
    'q is an eigenvalue'
    flagg=1;
    end

y=(A-q.*eye(n,n))\x; %Matlab's backslash operator
mu=y(p);

p=find(abs(y)==norm(y,inf),1);
er=norm((x-(y./y(p))),inf);
x=y./y(p);

%output all data into a matrix where first three coulmns are correspond to the 
%eigenvector approximations located across the rows and the fourth coulmn
%coresponds to the eigenvalue approximations the most recent is on the last
%row
l(k,:)=x;
M(k,1)=mu; 
lM=[l M]; 
if er<tol
    flagg=1;
end
k=k+1;
end
