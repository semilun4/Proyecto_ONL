% Lenguaje: Matlab
% Programa del Método de Newton para optimización en varias variables.
% Nombre del archivo: newton_mult.m
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai  
%
% Definicion: [opt,fopt,cont_met,llamadas_f] = newton_mult(x0)
%   
%  Entradas 
%   x0  = punto inicial (ejemplo: [0,-1])
%
%  Salidas
%   opt      = punto optimo local
%   fopt     = evaluacion del optimo local
%   cont_met = veces que se realiza el metodo

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Querida Semi del futuro, el programa funciona bien si le das 
% un punto ANTES del optimo. Saludos cordiales.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [opt,fopt,llamadas_f1] = newton_mult(x0,eps1)
tic() 
global fcalls
 fcalls = 0;
  h = 0.001;
  %f0 = funcion(x0);
  g0 = grad(x0,h);
  if norm(g0) > eps1
      H0 = hessiana(x0,h);
      cont_met = 0;
      while norm(g0) >= eps1
          d = H0\g0';
          x1 = x0' - d;
          fx1 = funcion(x1);
          x0 = x1';
          g0 = grad(x0,h);
          H0 = hessiana(x0,h);
          cont_met = cont_met + 1;
      end
      opt = x0;
      fopt = fx1;
  else
      opt = [NaN, NaN];
      fopt = NaN;
  end
      llamadas_f1 = fcalls; 
      toc()
end
     