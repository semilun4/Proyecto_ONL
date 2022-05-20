% Lenguaje: Matlab
% Programa para evaluar una funcion
% Nombre del archivo: funcion.m
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai 
%
% Definición: [F] = funcion(x)
% 
% Entrada: 
%   x = valor en el que se quiere evaluar
% Salida
%   F = evaluacion de la funcion
%
% Ejemplo de uso: funcion(1)
%
function [F] = funcion(x)
% Es una variable global que cuenta num de evaluaciones
 global fcalls 
% Funciones univariables
  % F = (100-x)^2;
 
% Funciones Newton
  % F = (x(1)+10*x(2))^2 + 5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
  % F = x(1)^2 + 2*x(2)^2 + 3*x(3)^2 + 4*x(4)^2 + (x(1) +x(2) +x(3) +x(4))^2;
  % F = (x(1)-1)^2 +x(2)^3 -x(1)*x(2);
% Cuasi-Newton funciones
  % F = 100*(x(1)^2 - x(2)^2) + (1-x(1))^2;
 
% Proyecto
   % F = sin(x(1)+ x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1; %McCormick
   % F =((x(1)^2)/4000)+((x(2)^2)/4000)-cos(x(1))*cos((x(2))/(sqrt(2)))+1;%Griewank
    F=(-1/10000)*(abs(sin(x(1))*sin(x(2))*exp(abs(100-sqrt(x(1)^2+x(2)^2)/pi)))+1)^0.1;%
   % F = (x(1)-1)^2 + (x(2)-1)^2 - x(1)*x(2); %Tried
   % F =((4-2.1*x(1)^2+((x(1)^4)/3))*x(1)^2)+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2; %Camel
   % F = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2); %Easom
 fcalls = fcalls+1;
end
