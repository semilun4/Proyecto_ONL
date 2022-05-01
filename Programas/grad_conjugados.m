% ------------------------------------------------------------------------
% %%%%% PROGRAMA QUE IMPLEMENTA EL ALGORITMO DE GRADIENTE CONJUGADO %%%%%%
% ------------------------------------------------------------------------
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai
% ------------------------------------------------------------------------
% Se hara uso del metodo  de seccion dorada para calcular el 
% tamaño de paso.
% Se detiene cuando la norma del gradiente es menor que 'eps1'.
% La función objetivo, el tamaño de paso y el gradiente (numerico) se 
% calculan en funciones auxiliares
%--------------------------------------------------------------------------
%   PARAMETROS DE ENTRADA:
%
%   xm      = punto inicial
%   eps     = tolerancia para la condición de paro,|(g(x)| < eps1
%   
%   SALIDA:
%   x0       = punto optimo 
%   f0       = funcion evaluada en el punto optimo
%   iter     = numero de iterciones con la que se optuvo el punto optimo
%   funcalls = numero de veces en que se evalua la funcion
%     
%   Ejemplo de uso: desde la línea de comandos introduce la definición de la
%                   función objetivo f y llama al método como se muestra: 
%
%     > [x0,fx0,num_iter,llamadas_f] = Grad_conj_Arm_bibi([1,1],0.01)

function [x0,fx0,llamadas_f] = grad_conjugados(x0,eps)
tic()
global fcalls
fcalls = 0;
h      = 0.01;

c      = 0.01;
alpha  = 1;
k      = 0;
gx0    = grad(x0,h);
d0 = -gx0;
while norm(gx0) >= eps
    gx0  = grad(x0,h); %gradiente en el punto incial
    nu   = -gx0; %direccion del nuevo paso
    x1   = x0 + alpha*nu; 
    fx1  = funcion(x1);
    fx0  = funcion(x0);
    %Calculo del tamaño de paso con Armijo
    
    while  fx1-fx0 > c*alpha*gx0'*nu
        alpha  = alpha * 0.5;
        x1     = x0 + alpha*nu;
        fx1    = funcion(x1);
    end
        x1   = x0 + alpha*nu;
        fx1  = funcion(x1);
        x0   = x1;
        fx0  = fx1;
        gx   = grad(x0,h);
        k    = k +1;
    if norm(gx0) <= 10e-05
        break
    else
        %d0   = -gx0; %direccion inicial
        %A    = (gx0*gx0');
        %B    = (gx*(gx-gx0)');
        %beta = B/A;% calculo de beta usando la formula Pollak-Ribiere
        beta = (gx*gx')/(gx0*gx0'); 
        d1   = -gx+beta*d0'; %calculo de la nueva direccion
        d0   = d1;
        if k>2
            d0 =-gx0;
            nu = -gx0;
            k=0;
        end
    end
    k   = k+1;
    %nu   = -gx0;
    x0  = x1;
    gx0 = gx;
    fx0 = fx1;
end
num_iter   = k; %Numero de iteraciones
llamadas_f = fcalls; 
toc()
end
