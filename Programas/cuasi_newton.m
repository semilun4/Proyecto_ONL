%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% García de la Cruz Semiramís
% De la Torre Ortiz Bibiana
% Bautista Lopez Sara
% Nicolas Palacion Daniel Isai 
% ------------------------------------------------------------------------
% %%%%% PROGRAMA QUE IMPLEMENTA EL ALGORITMO DE CUASI NEWTON RANGO 1%%%%%%%%%
% ------------------------------------------------------------------------
% Se hara uso del metodo del ranga 1 para hallar la inversa de la hessiana
% Se hara uso del metodo  de seccion dorada para calcular el 
% tamaño de paso.
% Se detiene cuando la norma del gradiente es menor que 'eps1'.
% La función objetivo, el tamaño de paso y el gradiente (numerico) se 
% calculan en funciones auxiliares
%--------------------------------------------------------------------------
%   PARAMETROS DE ENTRADA:
%
%   f       = función objetivo
%   xm      = punto inicial
%   epsilon = tolerancia para la condición de paro,|(g(x)| < eps1
%   
%   SALIDA:
%   xm       = punto optimo 
%   fm       = funcion evaluada en el punto optimo
%   iter     = numero de iterciones con la que se optuvo el punto optimo
%   funcalls = numero de veces en que se evalua la funcion
%     
%   Ejemplo de uso: desde la línea de comandos introduce la definición de la
%                   función objetivo f y llama al método como se muestra: 
%
%   >  [xm,fx0,llamadas_f]=cuasi_newton([1,1],0.01)
%

function [xm,fx0,llamadas_f]=cuasi_newton(xm,epsilon)
global fcalls;
fcalls = 0;
 tic()
K = length(xm); %mide le nuero de variables de a funcion 
h = epsilon; %tolerancia de error 
iter = 0;   
h0 = eye(K); %Matriz H inicial
grad0=grad(xm,h);
while norm (grad0)> epsilon
    grad0=grad(xm,h);
    d0 = (-h0*grad0')';
    f_nu = @(x) funcion(xm+x*d0);
    [alpha,itertp] = seccion_dorada(0,1,epsilon,f_nu); %se calcula el tamaño de paso con seccion dorada
    xn   = xm + alpha*d0;
    xm   = xn;
    fx0  = funcion(xm);
    grad1 = grad(xm,h);
    if norm(grad1)> epsilon
        deltax = alpha*d0';
        deltag = grad1 - grad0;
        A = deltax-h0*deltag';
        B = deltag*A;
        h1 = h0 +((A*A')/B);
        h0=h1;
    end
    iter = iter+1;
end
llamadas_f=fcalls;
toc()
end

%% FUNCION QUE CALCULA EL PUNTO MINIMO DE UNA FUNCION DE UNA VARIABLE
function [xn,itertp ]=seccion_dorada(x0,delta,epsilon,f)
itertp =0;
[a,b]=acotamiento(x0,delta,f);
L=b-a;
tau=(-1+sqrt(5))/2;
x1=b-tau*L;
fx1=f(x1);
x2=a+tau*L;
fx2=f(x2);
while (L>epsilon)
    if fx1<fx2
        b=x2;
        x2=x1;
        fx2=fx1;
        L=b-a;
        x1=b-tau*L;
        fx1=f(x1);
    else
        a=x1;
        x1=x2;
        fx1=fx2;
        L=b-a;
        x2=a+tau*L;
        fx2=f(x2);
    end
    itertp = itertp+1;
    L=b-a;
    xn=(b+a)/2;
end
end
%%
function [a,b]= acotamiento(x0,delta,f)
x1=x0-delta;
x2=x0+delta;
fx1=feval(f,x1);
fx2=feval(f,x2);
fx0=feval(f,x0);
k=2;
%determinamos el signo de delta 
if fx1>=fx0 & fx0>=fx2    % en este caso delta es positiva
    xk1 =x2;              % se calcula el nuevo el extremo izquierdo del intervalo
    fxk1=fx2;
    xk2 =xk1+2*delta;     % se calcula el extremo derecho del nuevo intervalo
    fxk2=feval(f,xk2);
    while fxk1>fxk2       % mientras fxk2 disminuya respecto a fxk1 se completara un ciclo nuevo
        xk1 =xk2;
        fxk1=fxk2;
        xk2 =(xk2+(2^k)*delta);
        fxk2=feval(f,xk2);
        k=k+1;
    end
    a=xk1-(2^(k-2))*delta;
    b=xk2;
    iter=k-2;
elseif fx1>=fx0 & fx0<=fx2 % en este caso el intervalo esta ente x1 y x2
    a=x1;
    b=x2;
    iter=k-1;
elseif fx1<=fx0 & fx0<=fx2 %en este caso delta es negativa
    xk1 =x1;              % se calcula el nuevo el extremo izquierdo del intervalo
    fxk1=fx1;
    xk2 =xk1-2*delta;      % se calcula el nuevo el extremo derecho del intervalo
    fxk2=feval(f,xk2);
    k=2;
    while fxk1>fxk2        % mientras fxk2 disminuya respecto a fxk1 se completara un ciclo nuevo    
        xk1 =xk2;           
        fxk1=fxk2;
        xk2 =(xk2-(2^k)*delta);
        fxk2=feval(f,xk2);
        k=k+1;
    end
    b=xk1+(2^(k-2))*delta;
    a=xk2;
    iter=k-2;
end
end



            

    
    
    
   
    
    
    
