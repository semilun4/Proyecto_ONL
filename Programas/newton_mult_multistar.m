%Lenguaje: Matlab
%Programa para calcular optimos locales utilizando el metodo de 
%Newton en un mallado de [-10,10)x[-10,10]
%Nombre del archivo: newton_mult_multistar
%García de la Cruz Semiramís
%De la Torre Ortiz Bibiana
%Bautista Lopez Sara
%Nicolas Palacios Daniel Isai
%
%Ejemplo de uso: Correr el programa 
%Para elegir la funcion al final descomentar la funcion que se desea
%saber sus puntos optimos dentro del mallado.
%Se generara un archivo .txt donde se encuentran los 100 puntos hallados.
global fcalls
fcalls = 0;
n=1;
for i=1:100 
  x1 = linspace(-10,10,100);
  x2 = linspace(-10,10,100);
  
  for j=1:100
        p(n,:)=([x1(i),x2(j)]);
        n=n+1; 
  end 
  [xp(i,:)]=newton_mult(p(i,:));
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
function [opt] = newton_mult(x0) 
  eps1 = 10e-05;
  h = 0.001;
  %f0 = funcion(x0);
  g0 = grad(x0,h);
  H0 = hessiana(x0,h);
  cont_met = 0;
  while norm(g0) > eps1
      d = H0\g0';
      x1 = x0' - d;
      fx1 = funcion(x1);
      x0 = x1';
      g0 = grad(x0,h);
      H0 = hessiana(x0,h);
      cont_met = cont_met + 1;
  end
  opt = x1;
end
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
   % F = sin(x(1)+ x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1; %McCormick
    %F = ((x(1)^2)/4000)+((x(2)^2)/4000)-cos(x(1))*cos((x(2))/(sqrt(2)))+1;
   % F = x(1) + 2*x(2); %Rotated Hyper-Ellipsoid
   % F = (x(1)-1)^2 + (x(2)-1)^2 - x(1)*x(2); %Tried
   %%%%F = ((4-2.1*x(1)^2+((x(1)^4)/3))*x(1)^2)+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;
   %F = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2); %Easom
   F=(-1/10000)*(abs(sin(x(1))*sin(x(2))*exp(abs(100-sqrt(x(1)^2+x(2)^2)/pi)))+1)^0.1;%Funcion 3
 fcalls = fcalls+1;
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