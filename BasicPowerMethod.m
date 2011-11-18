%%%%% Basic Power method%%%%%% 
%%Author: Colbert Sesanker
function [lM]= P(x,A,it,tol)
p=find(abs(x)==norm(x,inf),1);
x=x./x(p);
k=1;
flagg=0;
while (k<it)&&(flagg==0)
y=A*x;
mu=y(p);
p=find(abs(y)==norm(y,inf),1);
if y(p)==0
    'eigenvector x eigenvalue 0'
end
er=norm((x-(y./y(p))),inf);
x=y./y(p);
l(k,:)=x;
M(k,1)=mu; 
lM=[l M]; 
if er<tol
    flagg=1;
end
k=k+1;
end

% For Example: type [lM]=P([1 2 1]',[4 2 1; 0 3 2; 1 1 4],4,10^-3) into the Interpreter
% [lM]=P([1 2 1]',[4 2 1; 0 3 2; 1 1 4],4,10^-3)
% 
%                           2a
% lM =
% 
%     1.0000    0.8889    0.7778    4.0000         This row corresponds to the eigen vector with the largest eigenvalue
%     1.0000    0.6441    0.7627    6.5556
%     1.0000    0.5714    0.7759    6.0508
% 
% >> [lM]=IP([1 2 1]',[4 2 1; 0 3 2; 1 1 4],4,10^-3)
% 
%                           4a
% lM =
% 
%     1.0000    0.4211    0.8195    0.6412
%     1.0000    0.5752    0.8013    1.6170
%     1.0000    0.5526    0.8130    1.7097
