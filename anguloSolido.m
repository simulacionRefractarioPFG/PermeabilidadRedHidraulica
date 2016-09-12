function omega = anguloSolido(vert,P1,P2,P3)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Oosterom and Strakee formula
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = [P1(1)-vert(1), P1(2)-vert(2), P1(3)-vert(3)]; 
    b = [P2(1)-vert(1), P2(2)-vert(2), P2(3)-vert(3)]; 
    c = [P3(1)-vert(1), P3(2)-vert(2), P3(3)-vert(3)]; 
    a_norma = sqrt(a(1)^2 + a(2)^2 + a(3)^2);
    b_norma = sqrt(b(1)^2 + b(2)^2 + b(3)^2);
    c_norma = sqrt(c(1)^2 + c(2)^2 + c(3)^2);
    
    numerador   = abs(dot(a,cross(b,c)));
    denominador = a_norma*b_norma*c_norma + dot(a,b)*c_norma ...
                + dot(c,a)*b_norma + dot(b,c)*a_norma;
            
%     a_u = a./a_norma;
%     b_u = b./b_norma;
%     c_u = c./c_norma;
%     numerador   = abs(dot(a_u,cross(b_u,c_u)));
%     denominador = 1 + dot(a_u,b_u) + dot(c_u,a_u) + dot(b_u,c_u);
            
    arcoTan = atan(numerador/denominador);
    if(arcoTan<0)
        arcoTan = arcoTan + pi;
    end
    omega = 2*arcoTan;
    
end