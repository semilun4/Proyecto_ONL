% Lenguaje: Matlab
% Programa para calcular numericamente la matriz hessiana.
% Componemos la hessiana a partir de gradiente.
% Semi del futuro: f -> g (1x2) = f1 f2 -> g(f1 f2) = H (4x4)
% Nombre del archivo: hessiana.m
%
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai  
%
% Definicion: H = hessiana(x,h)
%   
% Entradas 
%   x = punto en el que se quiere evaluar
%   h = tamaño del paso
% Salidas
%  H = matriz hessiana
%
% Ejemplo de uso: H = hessiana(@(x) x(1)^2+x(2)^3,[1,1],0.001)
%
%  Matriz Hessiana
%  uso
%  y = hessi(x,h)
%  f funcion de la que se deasea conocer la hessiana
%  x vector donde se desea encontrar la hessiana
%  h tam. de paso
%
%  NOTA: Se llama al gradiente de la funcion 2n veces
%        donde n es length(punto).
%        La funcion se evalua (2n)*(2n)*n

function H = hessiana(x,h)
global nhess
    n=length(x); 
    for i=1:n
        x1 = x;
        x1(i) = x(i) - h ;
        df1 = grad(x1,h);
        x2 = x;
        x2(i) = x(i) + h ;
        df2 = grad(x2,h);
        d2f = (df2-df1) / (2*h);
        H(i,:) = d2f';
    end
 nhess = nhess +1;
end