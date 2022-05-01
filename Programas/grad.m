% Lenguaje: Matlab
% Programa para calcular numericamente el gradiente de una funcion
% Se implementa el metodo diferencias divididas centradas
% Nombre del archivo: grad.m
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai
%
% Definicion: [grad] = grad(f,x,h)
%   
% Entradas 
%   x = punto en el que se quiere evaluar el gradiente
%   h = tamaño del paso
% Salidas
%  grad = vector gradiente que resulta de la aproximacion
%
% Ejemplo de uso: grad([1,1],0.001)
% NOTA: Se evalua f 2n veces cada que se calcula el gradiente.

function G = grad(x, h)
 global ngrad
    n = length(x);% Dimensión del vector gradiente
    I = (h/2)*eye(n);% Vectores con incrementos de x1,...,xn
    for i=1 : n
       % Formula de Diferencias Finitas
       G(i) = (funcion(x+I(i,:))-funcion(x-I(i,:)))/h;
    end
  ngrad = ngrad +1;
end
  