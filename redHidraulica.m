function [Caudal, perdidaCarga, LongitudDisco, ...
          AreaSeccionDisco] = redHidraulica(filename)


% filename = 'dumpUNIFORME';
[P, VERT, CON, J] = tet_voroVert(filename);

disp('Asociados vertices y tetraedros')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matriz de distancias entre poros [um]
disp('Calculando distancias entre poros')
Nvert = length(VERT);
Ncon  = length(CON);
L = sparse(Nvert,Nvert);
for nr = 1:Ncon
    i = CON(nr,1);
    j = CON(nr,2);
    L(i,j) = sqrt((VERT(i,1)-VERT(j,1)).^2 + ...
                  (VERT(i,2)-VERT(j,2)).^2 + ...
                  (VERT(i,3)-VERT(j,3)).^2);
              
    L(j,i) = L(i,j);
end
disp('Distancias entre poros calculadas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matriz con las conductividades entre poros
viscosidad = 1002; % [uPa*s] a 20 celsius
disp('Calculando conductividades entre poros')
[G, RadiosHidraulicos] = conductividad(P, VERT, J, viscosidad);
disp('Conductividades entre poros calculadas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Red hidraulica: analogia con una red electrica con
%%% resistencias entre nodos.
%%% [K][P] = [Q]
%%% donde K = G./L
%%% con las condiciones de contorno de tipo Dirichlet para los nodos de la
%%% cara superior e inferior, donde se conoce el valor de la presion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Construye la matriz de conductividad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = find(G);
K = sparse(Nvert, Nvert);
% Divide solo cuando hay algun elemento en G, para evitar
% una matriz resultado llena de NaN's.
K(w) = G(w) ./ L(w); % [um3 * uPa-1 * s-1]

% Los elementos calculados que no estan en la diagonal toman valor negativo
K = -K;
% Los elementos de la diagonal son la suma de todas la conductividades que
% llegan al nodo i, con valor positivo
for i = 1:Nvert
    K(i,i) = abs(sum(K(i,:)));
end
disp('Matriz de conductividad construida')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Presion = sparse(Nvert,1);
% Pverts es un vector de presiones que no se modificara al aplicar
% las condiciones de contorno, al contrario que P
Pverts  = zeros(Nvert,1); 
Q       = zeros(Nvert,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONDICIONES DE CONTORNO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Aplicando condiciones de contorno')
% Cuidado con no asignar una presion de cero en ninguna de las caras para
% poder diferenciar los poros que tienen asignadas condiciones de contorno
PcaraSup     = 3*10^(11); % 3e11 uPa -> 300000 Pa -> 3 bar
perdidaCarga = 1*10^(11); % 1e11 uPa -> 100000 Pa -> 1 bar
% PcaraSup     = 2e11; % 2 bar = 2*10^5 [Pa] => [Pa]*10^6 = [uPa]
% perdidaCarga = 1e11; % 1 bar = 1*10^5 [Pa] => [Pa]*10^6 = [uPa]
PcaraInf     = PcaraSup - perdidaCarga;

% Condiciones de contorno para los poros que se situen por debajo
% de la altura minima de una de las particulas,
% tambien para aquellos poros que superen la altura maxima de las
% particulas, los cuales tienen una presion fija correspondiente a la cara
% superior.
inputVoroFile = [filename,'InputVoro.txt'];
% lee inputVoroFile
DUMP  = dlmread(inputVoroFile);
% borra inputVoroFile
system(sprintf('rm %s', inputVoroFile));

x      = DUMP(:,2);         % posicion x del centro de la particula
y      = DUMP(:,3);         % posicion y del centro de la particula
z      = DUMP(:,4);         % posicion z del centro de la particula
r      = DUMP(:,5);         % radio de la particula
z_max = z+r;
z_min = z-r;
z_top = max(z_max);
z_inf = min(z_min);
LongitudDisco = z_top - z_inf;
distanciaCentroDisco = sqrt(x.^2 + y.^2);
diametroDisco = max(distanciaCentroDisco);
AreaSeccionDisco  = pi*(max(diametroDisco/2))^2;

Presion(VERT(:,3)>z_top) = PcaraSup;
Presion(VERT(:,3)<z_inf) = PcaraInf;
Pverts(VERT(:,3)>z_top)  = PcaraSup;
Pverts(VERT(:,3)<z_inf)  = PcaraInf;
NporosSuperiores = length(Presion(VERT(:,3)>z_top));
NporosInferiores = length(Presion(VERT(:,3)<z_inf));
NporosContorno   = NporosSuperiores + NporosInferiores;

indiceContornoSup = find(Presion == PcaraSup);
indiceContornoInf = find(Presion == PcaraInf);
indiceContorno    = [indiceContornoSup; indiceContornoInf];

% vector guardando los indices de las filas para poder reconocer
% la posicion de los poros despues de haber aplicado las condiciones de
% contorno
indices = (1:Nvert)';
for i = 1:Nvert
    for j = 1:NporosContorno
        % Q [um3 * s-1]
        Q(i) = Q(i) - K(i,indiceContorno(j)).*Presion(indiceContorno(j));
    end
end

% Elimina filas y columnas asociadas a los poros con condiciones de
% contorno
K(indiceContorno, :)    = [];
K(:, indiceContorno)    = [];
Presion(indiceContorno) = [];
Q(indiceContorno)       = [];
indices(indiceContorno) = []; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Mean normalization
% K_mean = sum(K(:))/(length(K(:)));
% Q_mean = sum(Q(:))/(length(Q(:)));
% 
% K = (K-K_mean)./(max(K(:))-min(K(:)));
% Q = (Q-Q_mean)./(max(Q(:))-min(Q(:)));


disp('Resolviendo el sistema de ecuaciones...')

% Resuelve el sistema de ecuaciones lineales
Presion = K \ Q;

disp('Sistema resuelto')

residual = K*Presion - Q; % calculo del residuo

% Recupera los valores de las presiones
% Presion = Presion*(1*10^(ordenQ-ordenK));
Pverts(indices) = Presion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Caudal = 0;
for i = 1:NporosInferiores
    Caudal = Caudal + K(indiceContornoInf(i),:)*Presion;
end
Caudal

indicePresionMedia = find(Presion >= (PcaraSup-PcaraInf)/2.5 & Presion <= (PcaraSup-PcaraInf)/1.5);
Caudal1 = 0;
for i = 1:length(indicePresionMedia)
    for j = 1:length(Presion)
        if(K(indicePresionMedia(i),j)*Presion(j) > 0)
            Caudal1 = Caudal1 + K(indicePresionMedia(i),j)*Presion(j);
        end
    end
end
Caudal1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dibuja todos los pares de uniones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Dibujando conexiones')
for i = 1:length(CON)
    plot3([VERT(CON(i,1),1),VERT(CON(i,2),1)],...
          [VERT(CON(i,1),2),VERT(CON(i,2),2)],...
          [VERT(CON(i,1),3),VERT(CON(i,2),3)], 'linewidth', 2.5);
    hold on
end
axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Asocia a cada presion un color dentro del mapa de colores VIRIDIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valoresPresion = sort(unique(Pverts));
Mcolor         = viridis();
PresionColor   = zeros(Nvert, 3);
for poro = 1:Nvert
    if(Pverts(poro)~=0)
        indicePresionOrdenada = find(valoresPresion == Pverts(poro));
        % redondea a uno de los 256 colores
        indiceColor = ceil((indicePresionOrdenada/length(valoresPresion))*255); 
        color = Mcolor(indiceColor, :);
        PresionColor(poro, :) = color;
    end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Asocia a cada union el radio hidraulico que le corresponde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radioCilindros    = zeros(length(CON), 1);
for union = 1:length(CON)
    radioCilindros(union) = RadiosHidraulicos(CON(union,1), CON(union,2)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


POROS = [VERT, Pverts, PresionColor];
CON   = [CON, radioCilindros];


dlmwrite('POROS.txt', POROS, 'delimiter', ' ');
dlmwrite('CON_poros.txt', CON, 'delimiter', ' ');

disp('Tareas completadas')

end
