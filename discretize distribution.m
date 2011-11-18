#!/usr/bin/octave

#Author Colbert Sesanker 9/15/11
#

function thrown = basket(points,a,b,bins)
dx = (b-a)/bins;
% x = a:1/grid:b;
total = quad("f",a,b);              % Integrate distribution f, defined below, over [a,b]
k=0;

for i=1:bins
  v(i)=quad("f",a+(i-1)*dx,a+i*dx); % Vector of integrals of each grid interval of distribution
  w(i)=round(v(i)*points/total);
    for j= k+1 : k + w(i)
	r(j)= a + (2*i-1)*dx/2;
	k=length(r);
	endfor
endfor

basket = f(r);
pick = round(1-rand() + rand()*length(basket));

thrown = basket(:);
endfunction

function y = f(x)                  % Distribution to be disrcretized
y = exp(-x.^2);
endfunction


