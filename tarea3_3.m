% Matriz de adyacencia del grafo
M = [
    0 30 0 0 30 0 0 0 0 0 0 0;
    0 0 20 0 0 20 0 0 0 0 0 0;
    0 0 0 30 0 30 40 0 0 0 0 0;
    0 0 0 0 0 0 0 40 0 0 0 0;
    0 0 0 0 0 40 0 0 20 0 0 0;
    0 0 0 0 0 0 40 0 0 0 0 0;
    0 0 0 0 0 0 0 30 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 50;
    0 0 0 0 0 0 0 0 0 50 0 0;
    0 0 0 0 0 0 30 0 0 0 30 0;
    0 0 0 0 0 0 0 0 0 0 0 10;
    0 0 0 0 0 0 0 0 0 0 0 0
];

% Arcos y costos
[i,j,c] = find(M);

% Número de nodos
n = size(M,1);
% Número de arcos
m = length(c);

% Crear el vector de flujos
b = zeros(n,1);
b(1)  = 1; % Fuente
b(12) = -1; % Finalización

% Matriz de coeficientes
A = zeros(n,m);
for k = 1:m
    % Nodos de salida
    A(i(k),k) = 1;
    % Nodos de entrada
    A(j(k),k) = -1;
end

% Valor mínimo: 0
lb = zeros(m,1);

% linprog(f, A, b, Aeq, beq, lb, ub)
% Aquí solo hay restricciones de igualdad entonces de deja [] en A,b
[x,fval] = linprog(c,[],[],A,b,lb,[]);

% Mostrar valor óptimo
fprintf('Costo: %g\n', fval);

% Arcos seleccionados
arcos = find(x > 1e-5);
arcos = [i(arcos) j(arcos)];
fprintf('Arcos seleccionados: \n');
disp(arcos);

% Corroborar con Dikjstra
G = digraph(M);
[path, cost] = shortestpath(G, 1, 12);

fprintf('Costo (Dijkstra): %g\n', cost);
fprintf('Recorrido (Dijkstra): \n');
disp(path);