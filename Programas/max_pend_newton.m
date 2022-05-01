% Lenguaje: Matlab
% Programa para implementar el algoritmo de máxima pendiente
% usando el metodo de Newton para calcular el tamaño de paso.
% Nombre del archivo: max_pend_newton.m
% % García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai 
%
% Definicion: [opt,fopt,it_pend,it_alpha,llamadas_f] = max_pend_newton(x0,eps)
%  
%  Entradas
%   x0  = punto inicial
%   eps = tolerancia
%  Salidas
%   opt      = punto optimo local
%   fopt     = evaluacion del optimo local
%   it_pend  = iteraciones del algoritmo de max. pendiente
%   it_alpha = iteraciones para calcular el tamaño de paso
%
% Ejemplo de uso:
%  f esta en el archivo de nombre: funcion.m  y esta definida
%  como F = 3*x(1)^2 + 3*x(1)*x(2) + 4*x(2)^2, por ejemplo.
%  [opt,fopt,iter_pend,it_aplha,llamadas] = max_pend_newton([-4,1],0.001)
%
function [opt,fopt,llamadas_f] = max_pend_newton(x0,eps)
tic()
global fcalls
global paso
 fcalls = 0; paso = 0;
 it_pend = 0;
 h = 0.001;
 gx0 = grad(x0,h);
 while norm(gx0) > eps
     nu = -gx0;  % Dirección de descenso.
     f_nu  = @(x) funcion(x0+x*nu); % Composición de la funcion
   % Se construye la derivada para sacarle su raíz
     der_fnu = @(x) PDD_center(f_nu,x,h);
   % Se calcula el valor de alpha con el metodo de Newton
   % univariable, es decir, se estima el tamaño de paso
     alpha = newton_univar(der_fnu,0,h);
     x1 = x0 + alpha*nu;
     fx1 = funcion(x1);
     x0 = x1;
     gx0 = grad(x0,h);
     it_pend = it_pend+1;
 end
 opt = x1;
 fopt = fx1;
 it_alpha = paso;
 llamadas_f = fcalls;
 toc()
end
% Primeras diferencias divididas centradas
 function [df] = PDD_center(f,x,h)
% Funcion para aproximar la derivada por DD
 df=(f(x+h/2)-f(x-h/2))/h;
 end
