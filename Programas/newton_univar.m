%Lenguaje: Matlab
%Programa para calcular la raiz de una funcion aplicando el metodo
%de Newton para una variable. Para realizar la derivada numericamente,
%se usara la primer diferecia dividida centrada.
%Nombre del archivo: aprox_der.m
%Semiramis G. de la Cruz
%De la Torre Ortiz Bibiana
%Bautista Lopez Sara
%Nicolas Palacios Daniel Isai
%
%Definicion: [values,df] = exact_der(a,b,f,v)
%   
% Entradas 
%   x = punto en el que se quiere evaluar
%   h = tamaño del paso
% Salinas
%  df = valor aproximado de la derivada
%
%Ejemplo de uso: aprox_der(3,1)
%
function [xn] = newton_univar(f,x,eps)
global paso
     eps1 = eps;
     cond = true;
     % ciclo do-until, hecho con un while()
     while cond
        x  = x-f(x)/PDD_center(f,x,eps1);
        xn = x;
       cond = (abs(f(x))>eps);
     end
     paso = paso +1; 
end

function df = PDD_center(f,x,h)
  df=(f(x+h/2)-f(x-h/2))/h;
end