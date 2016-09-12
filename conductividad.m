function [G, RadiosHidraulicos] = conductividad(P, VERT, J, viscosidad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calcula las conductividades entre poros(nodos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nvert = length(VERT);

G                 = sparse(Nvert, Nvert);
RadiosHidraulicos = sparse(Nvert, Nvert);

[Nrows, ~] = size(J); 
for i = 1:Nrows
    
if(J(i,1)~=0) % si la fila contiene informacion de una conexion
    % Calcula el volumen de cada uno de los tetraedros formados por los
    % poros y los vertices que forman la cara perpendicular a la conexion
    % entre poros.

    M1 = [VERT(J(i,1), 1), VERT(J(i,1), 2), VERT(J(i,1), 3), 1; ...
          P(J(i,3), 1), P(J(i,3), 2), P(J(i,3), 3), 1; ...
          P(J(i,4), 1), P(J(i,4), 2), P(J(i,4), 3), 1; ...
          P(J(i,5), 1), P(J(i,5), 2), P(J(i,5), 3), 1];
    M2 = [VERT(J(i,2), 1), VERT(J(i,2), 2), VERT(J(i,2), 3), 1; ...
          P(J(i,3), 1), P(J(i,3), 2), P(J(i,3), 3), 1; ...
          P(J(i,4), 1), P(J(i,4), 2), P(J(i,4), 3), 1; ...
          P(J(i,5), 1), P(J(i,5), 2), P(J(i,5), 3), 1];
      
    Vtet_1 = (1/6)*abs(det(M1));
    Vtet_2 = (1/6)*abs(det(M2));
    
    % Calcula el angulo solido para cada tetraedro.
    [omega3_1] = anguloSolido(P(J(i,3), :), VERT(J(i,1), :), ...
                              P(J(i,4), :), P(J(i,5), :));
    [omega4_1] = anguloSolido(P(J(i,4), :), VERT(J(i,1), :), ...
                              P(J(i,3), :), P(J(i,5), :));
    [omega5_1] = anguloSolido(P(J(i,5), :), VERT(J(i,1), :), ...
                              P(J(i,4), :), P(J(i,3), :));
    [omega3_2] = anguloSolido(P(J(i,3), :), VERT(J(i,2), :), ...
                              P(J(i,4), :), P(J(i,5), :));
    [omega4_2] = anguloSolido(P(J(i,4), :), VERT(J(i,2), :), ...
                              P(J(i,3), :), P(J(i,5), :));
    [omega5_2] = anguloSolido(P(J(i,5), :), VERT(J(i,2), :), ...
                              P(J(i,4), :), P(J(i,3), :));
                          
    % Angulos solidos reunidos segun su tetraedro asociado.
    omegas1 = [omega3_1;omega4_1;omega5_1];
    omegas2 = [omega3_2;omega4_2;omega5_2];
    
    % Volumen de la garganta que conecta los dos poros.
    volumenGarganta1 = Vtet_1 - (1/3)*omegas1'*P(J(i,3:5), 4).^3;
    volumenGarganta2 = Vtet_2 - (1/3)*omegas2'*P(J(i,3:5), 4).^3;
    volumenGarganta  = volumenGarganta1 + volumenGarganta2;
    
    % Area de contacto entre esferas y el volumen de la garganta
    % = suma de los angulos solidos por los radios al cuadrado.
    area_gamma1 = omegas1'*P(J(i,3:5), 4).^2; 
    area_gamma2 = omegas2'*P(J(i,3:5), 4).^2;
    area_gamma  = area_gamma1 + area_gamma2;
    
    RadioHidraulico = volumenGarganta/area_gamma;
    RadiosHidraulicos(J(i,1), J(i,2)) = RadioHidraulico;
    RadiosHidraulicos(J(i,2), J(i,1)) = RadioHidraulico;
    
    % Area de seccion minima en la garganta.
    
    % Distancias entre centros de las particulas que forman el triangulo.
    dist3_4 = sqrt((P(J(i,3), 1) - P(J(i,4), 1))^2 ...
                 + (P(J(i,3), 2) - P(J(i,4), 2))^2 ... 
                 + (P(J(i,3), 3) - P(J(i,4), 3))^2);
    dist4_5 = sqrt((P(J(i,4), 1) - P(J(i,5), 1))^2 ...
                 + (P(J(i,4), 2) - P(J(i,5), 2))^2 ... 
                 + (P(J(i,4), 3) - P(J(i,5), 3))^2);
    dist5_3 = sqrt((P(J(i,5), 1) - P(J(i,3), 1))^2 ...
                 + (P(J(i,5), 2) - P(J(i,3), 2))^2 ... 
                 + (P(J(i,5), 3) - P(J(i,3), 3))^2);
             
    % Angulos entre particulas, en el plano perpendicular a la garganta.
    % Uso del teorema del coseno.
    alpha3_4 = acos((dist4_5^2+dist5_3^2-dist3_4^2)/(2*dist4_5*dist5_3));
    alpha4_5 = acos((dist3_4^2+dist5_3^2-dist4_5^2)/(2*dist3_4*dist5_3)); 
    alpha5_3 = acos((dist4_5^2+dist3_4^2-dist5_3^2)/(2*dist4_5*dist3_4));
    
    % Area de los sectores circulares dados los angulos y los radios de
    % cada particula.
    areaSectores = (2*pi^2./[alpha3_4, alpha4_5, alpha5_3])*P(J(i,3:5), 4);
    
    % Area del triangulo. 
    base   = dist3_4;
    altura = sin(alpha4_5)*dist5_3;
    area_triangulo = (1/2)*(base*altura);
    
    Area_seccion = area_triangulo - areaSectores;
    
    % Conductividad en la union entre los dos poros [um^4/(uPa*s)]
    alpha = 1/2; % factor elegido en publicacion
    G(J(i,1), J(i,2)) = alpha*((Area_seccion*RadioHidraulico^2)/(viscosidad));
    G(J(i,2), J(i,1)) = alpha*((Area_seccion*RadioHidraulico^2)/(viscosidad));

end % end if
end % end for

end % end function