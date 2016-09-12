function filename = dump2voroInput(filename)
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lee archivo con la informacion de las particulas de la simulacion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DUMP   = dlmread(filename,' ',9,0); % skip 9 rows, 0 columns

id     = DUMP(:,1);      
x      = DUMP(:,3);         % posicion x del centro de la particula
y      = DUMP(:,4);         % posicion y del centro de la particula
z      = DUMP(:,5);         % posicion z del centro de la particula
r      = DUMP(:,18);        % radio de la particula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Escribe el archivo de salida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M      = [id, x, y, z, r];

% Ordena la posicion de las particulas segun su indice, de mayor a menor
P = sortrows(M, 1);

% Write
filename = [filename,'InputVoro.txt'];
dlmwrite(filename, P, 'delimiter', ' ')

end