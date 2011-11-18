%Automatic table generation procedure for GS, Jacobi, SOR function
%Colbert Sesanker
clear; clc; close

nV=[10 20 40 80 100 200];

wV=[-0.4 .4 .8 1.2 1.6 2.4];

alphaV=[.5 .3 -.3];


    for l=1:length(alphaV)
            for j=1:length(wV)
                            for i=1:length(nV) 
                                alpha=alphaV(l);
                                 w=wV(j); 
                                 n=nV(i);
                                if l==1
                                [k(j,i) er(j,i)]=it(n,alpha,100000,10^(-6),inf,7,w);
                                elseif l==2
                                [k2(j,i) er2(j,i)]=it(n,alpha,100000,10^(-6),inf,7,w);
                                else 
                                [k3(j,i) er3(j,i)]=it(n,alpha,100000,10^(-6),inf,7,w);                                                         
                                end           
                                
                            end
            end
    end
    
   L=[1 2];
   for j=1:2
       method=L(j);
   for l=1:length(alphaV)
       for i=1:length(nV) 
             alpha=alphaV(l);
                 n=nV(i);
                 if l==1
                     [P(j,i) erP(j,i)]=it(n,alpha,100000,10^(-6),inf, method);
                     
                 elseif l==2
                     [P2(j,i) erP2(j,i)]=it(n,alpha,100000,10^(-6),inf,method);
                 else
                     [P3(j,i) erP3(j,i)]=it(n,alpha,100000,10^(-6),inf,method); 
                 end
       end   
   end
   end
    P=[['J' 'GS'] [nV P']'];
   P2=[['J' 'GS'] [nV P2']'];
   P3=[['J' 'GS'] [nV P3']'];
        
    



%Results=[nV; k; flagS]
