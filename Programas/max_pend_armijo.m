% Lenguaje: Matlab
% Programa para implementar el algoritmo de máxima pendiente
% usando el metodo de eliminacion de regiones Seccion Dorada
% para calcular el tamaño de paso.
% Nombre del archivo: max_pend_armijo.m
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai  
%
% Definicion:
%  
%  Entradas
%   f = funcion objetivo
%   x0 = punto inicial
%   eps = tolerancia
%  Salidas
%   opt = punto optimo local
%   fopt = evaluacion del optimo local
%   iter = iteraciones del algoritmo de max. pendiente
%
% Ejemplo de uso:
%  f esta en el archivo de nombre: funcion.m  y esta definida
%  como F = 3*x(1)^2 + 3*x(1)*x(2) + 4*x(2)^2, por ejemplo.
%  [opt,fopt,iter_pend,it_aplha,llamadas] = max_pend_armijo([10,10],0.1)
%
function [opt,fopt,llamadas_f2] = max_pend_armijo(x0,eps)
tic()
global fcalls
global ngrad
ngrad = 0;
 fcalls = 0;
 it_alpha = 0;
 it_pend = 0;
 h = 0.001;
 gx0 = grad(x0,h); %2
 fx0 = funcion(x0); %1
 cl=10^-4;
 alpha=1;
 while norm(gx0) > eps
     %if it_pend<=n
     nu = -gx0;  % Dirección de descenso.
     % f_nu  = @(x) funcion(x0+x*nu); % Composición de la funcion
     % alpha=1;
     while funcion(x0+alpha.*nu)>fx0+cl*(alpha.*gx0'*nu) %4+16
        alpha=alpha/2;
        it_alpha = it_alpha +1;
     end
     x0 = x0 + alpha*nu;
     fx0 = funcion(x0); %16
     gx0 = grad(x0,h); %2*16
     %else
        %break;;
     %end
     it_pend = it_pend+1;
 end
 opt = x0;
 fopt = fx0;
 llamadas_f2 = fcalls;
 llamadas_g2 = ngrad;
 toc()
end

