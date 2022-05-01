%Lenguaje: Matlab
%Programa para calcular optimos locales utilizando el metodo de 
% Cuasi Newton en un mallado de [-10,10)x[-10,10]
%Nombre del archivo: cuasi_newton_multistar
%Semiramis G. de la Cruz
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
x=linspace(-10,10,100); %Generando los 100 valores en
y=linspace(-10,10,100); %Generando los 100 valores en
for i=1:1:100
    p(i,:)=[x(i),y(i)]; %Composicion del punto (x,y)
    [xp(i,:)]=cuasi_newton_armijo(p(i,:),10^-4);
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

%[xm,fx0,iter,llamadas_f]=cuasi_newton_armijo('funcion1',[1,1],0.01)
function [xm]=cuasi_newton_armijo(xm,epsilon)
global funcalls;
 tic()
funcalls =0;
K = length(xm); %mide le nuero de variables de a funcion
 
alpha=10
c=0.00001
h = epsilon; %tolerancia de error 
iter = 0;   
h0 = eye(K); %Matriz H inicial
grad0=grad(xm,h);
while norm (grad0)> epsilon
    grad0=grad(xm,h);
    d0 = (-h0*grad0')';
    xn   = xm + alpha*d0; 
    fxn  = funcion(xn);
    fxm  = funcion(xm);
    %Calculo del
    while  fxn-fxm > c*alpha*grad0'*d0
        alpha  = alpha * 0.5;
        xn     = xm + alpha*d0;
        fxn    = funcion(xn);
    end
    xn   = xm + alpha*d0;
    xm   = xn;
    fxm  = funcion(xm);
    grad1 = grad(xm,h);
    if norm(grad1)> epsilon
        deltax = alpha*d0';
        deltag = grad1 - grad0;
        %A = h0*deltag'
        %h1=h0+((deltax*deltax')/deltax'*deltag)-((A*A')/(deltag*A))
        A = deltax-h0*deltag';
        B = deltag*A;
        h1 = h0 +((A*A')/B);
        h0=h1;
    end
    iter=iter+1;
end
llamadas_f=funcalls;
toc()
end


%% FUNCION QUE CALCULA EL GRADIENTE
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

%%FUNCIONES DE PRUEBA%%
function y = funcion(x) %FUNCIÓN OBJETIVO
  global funcalls;
  %Funcion 3
  %y=(-1/10000)*(abs(sin(x(1))*sin(x(2))*exp(abs(100-sqrt(x(1)^2+x(2)^2)/pi)))+1)^0.1;
  % funcion camello de seis jorobas
  %y = (4-2.1*x(1)^2+(x(1)^4)/3)*(x(1)^2)+x(1)*x(2)+(-4+4*x(2)^2)*(x(2)^2);
 
  %funcion EASOM
  %y = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
  
  %Funcion Maccormick #26
  %y = sin(x(1)+x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1;
  
  % funcion Trid
  %y = (x(1)-1)^2 + (x(2)-1)^2 - x(1)*x(2); 
  
  %FUNCIÓN GRIEWANK #3
  %y = ((x(1)^2)/4000)+((x(2)^2)/4000)-cos(x(1))*cos((x(2))/(sqrt(2)))+1;
  %Funcion 3
  y=(-1/10000)*(abs(sin(x(1))*sin(x(2))*exp(abs(100-sqrt(x(1)^2+x(2)^2)/pi)))+1)^0.1;
   funcalls = funcalls + 1;
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