# Acerca del proyecto

En optimización no lineal, estudiamos métodos que encuentran el conjunto de valores para obtener la mejor solución o la más aproximada, la cual puede estar sujeta a una serie de restricciones, donde una o más variables incluidas es no lineal. 

En el presente proyecto, nos centraremos en aquellos problemas que no tienen restricciones y compararemos la eficacia de los métodos vistos, que son los siguientes
\begin{enumerate}
    \item Método de Máxima Pendiente \textit{(Cauchy)}
    \item Condición de Armijo
    \item Método de Newton
    \item Método de Cuasi-Newton
    \item Método de Gradientes Conjugados
\end{enumerate}
Para el análisis, tomaremos cuatro funciones de prueba, de las cuales tenemos resultados esperados para luego compararlos con los resultados obtenidos en cada programa. Además, tendremos en cuenta el tiempo y número de evaluaciones de la función objetivo, puesto que no sólo se requiere que los métodos sean exactos sino también eficientes, es decir, preferiremos un método que disminuya el costo computacional.

Por ser un problema de optimización, es evidente que la finalidad es que encontremos el óptimo de la función, dicho de otro modo, el mínimo global si es que existe. Cabe resaltar, que cada algoritmo arrojará la mejor aproximación a dicho punto.
