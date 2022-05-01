%Lenguaje: Matlab
%Programa para calcular optimos locales utilizando el metodo de 
% Maxima pendiente con armijo en un mallado de [-10,10)x[-10,10]
%Nombre del archivo: max_pend_armijomultistar
% García de la Cruz Semiramís
%De la Torre Ortiz Bibiana
%Bautista Lopez Sara
%Nicolas Palacios Daniel Isai
%
%Ejemplo de uso: Correr el programa 
%Para elegir la funcion al final descomentar la funcion que se desea
%saber sus puntos optimos dentro del mallado.
%Se generara un archivo .txt donde se encuentran los 100 puntos hallados.

n=1;
for i=1:100 
  x1 = linspace(-10,10,100);
  x2 = linspace(-10,10,100);
  
  for j=1:100
        p(n,:)=([x1(i),x2(j)]);
        n=n+1; 
  end 
  [xp(i,:)]=max_pend_armijo(p(i,:),10^-4);
    xop(i,:)=round(xp(i,:),3);
    fop(i,:)=funcion(xop(i,:));
end
[c]=Filtro(fop);
fprintf("Se encontraron %d puntos óptimos convergentes al valor del global\n", c);
fprintf("los puntos están en el archivo generado Puntos.txt en la ruta de la carpeta\n");
%Se escribirán los puntos en un archivo .txt 
fileID = fopen('Puntos.txt','w');
[m,n]=size(xop);
for i=1:1:m
    fprintf(fileID,'%f\t%f\t%f\n',xop(i,1),xop(i,2),fop(i,1));
end
fclose(fileID);
function [opt] = max_pend_armijo(x0,eps)
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

end
function[contador]=Filtro(x)    
fop=0;                   %minimo global
tol=0.003;                     %tolerancia 
pr=fop-tol;                       %valor minimo
ps=fop+tol;                       %valor maximo
contador=0;                     %contador de puntos minimos 
[m,n]=size(x);                  
for i=1:1:m
    if x(i,:)>=pr && x(i,:)<=ps   %si el punto dado está dentro del rango
         contador=contador+1;
    end
end
end
function G = grad(x, h )
 global ngrad
    n = length(x);% Dimensión del vector gradiente
    I = (h/2)*eye(n);% Vectores con incrementos de x1,...,xn
    for i=1 : n
       % Formula de Diferencias Finitas
       G(i) = (funcion(x+I(i,:))-funcion(x-I(i,:)))/h;
    end
  ngrad = ngrad +1;
end
function [F] = funcion(x)
% Es una variable global que cuenta num de evaluaciones
 global fcalls 
% Funciones univariables
  % F = (100-x)^2;
 
% Ejemplo profesora.
  % F = 3*x(1)^2 + 3*x(1)*x(2) + 4*x(2)^2; 
  % F = 4*x(1)+4*x(2)^3-10*x(2);

% Funciones Tarea 4
  % F = 8*x(1)^2 + 4*x(1)*x(2) +5*x(2)^2;
  % F = 2*x(1)^3 +4*x(1)*x(2)^3 - 10*x(1)*x(2) + x(2)^2;
  % F = 2*x(1)^2 - 1.05*x(1)^4 + (x(1)^6)/6 + x(1)*x(2) + x(2)^2;
  % F = (x(1)+2*x(2) -7)^2 + (2*x(1) + x(2) -5)^2;
  
% Funciones Tarea 6
  % F = 3*x(1)^2+2*x(2)^2+x(3)^2-2*x(1)*x(2)-2*x(1)*x(3)+2*x(2)*x(3)-6*x(1)-4*x(2)-2*x(3);
  % F =((x(2)-((5.1/(4*pi^2))*x(1)^2)+((5*x(1))/pi)-6)^2+((10*(1-(1/(8*pi))))*cos(x(1)))+10);
  % F = (x(2) - (5.1)/(4*pi^2)*x(1)^2+(5/pi)*x(1)-6)^2+10*(1-(1/(8*pi)))*cos(x(1))+10;

% Funciones Newton
  % F = (x(1)+10*x(2))^2 + 5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
  % F = x(1)^2 + 2*x(2)^2 + 3*x(3)^2 + 4*x(4)^2 + (x(1) +x(2) +x(3) +x(4))^2;
  % F = (x(1)-1)^2 +x(2)^3 -x(1)*x(2);
% Cuasi-Newton funciones
  % F = 100*(x(1)^2 - x(2)^2) + (1-x(1))^2;
 
% Proyecto
    %F = sin(x(1)+ x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1; %McCormick
    F = ((x(1)^2)/4000)+((x(2)^2)/4000)-cos(x(1))*cos((x(2))/(sqrt(2)))+1;%Griewak
   % F = x(1) + 2*x(2); %Rotated Hyper-Ellipsoid
   %F = (x(1)-1)^2 + (x(2)-1)^2 - x(1)*x(2); %Tried
   %F = ((4-2.1*x(1)^2+((x(1)^4)/3))*x(1)^2)+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;%camellos
   %F = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2); %Easom
   %F=(-1/10000)*(abs(sin(x(1))*sin(x(2))*exp(abs(100-sqrt(x(1)^2+x(2)^2)/pi)))+1)^0.1;%Funcion 3
 fcalls = fcalls+1;
end