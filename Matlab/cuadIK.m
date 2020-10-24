function [q, k] = cuadIK(Td, q0, DOF, IT_B, ET_F, varargin)
% CUADIK() Funci�n que obtiene la cinem�tica inversa del cuadr�pedo, ya sea
% para las patas traseras (2 DOF) o las patas delanteras (3 DOF). La
% cinem�tica inversa permite mapear del espacio de tarea al espacio de
% configuraci�n, por lo que se le alimenta una pose deseada y su config.
% inicial y este es capaz de computar la siguiente config.
% -------------------------------------------------------------------------
%
% Seg�n el n�mero de inputs adicionales se determina el modo
%   - 0 inputs = IK de posici�n usando Levenberg
%   - 1 input = IK de posici�n / posici�n-orientaci�n usando Levenberg
%   - 2 inputs = IK de posici�n / posici�n-orientaci�n usando "levenberg",
%     "pseudoinversa" o "traspuesta".
%
% -------------------------------------------------------------------------

switch numel(varargin)
    case 0
        solop = 0;
        metnum = "levenberg";
    case 1
        solop = varargin{1};
        metnum = "levenberg";
    case 2
        solop = varargin{1};
        metnum = varargin{2};
end    

tol_p = 1e-06;                          % Tolerancia del error de posici�n
tol_o = 1e-05;                          % Tolerancia del error de orientaci�n

Od = Td(1:3, 4);                        % Se extrae la posici�n y orientaci�n deseada de Td
Rd = Td(1:3, 1:3);
Qd = rot2cuat(Rd);
       
T = cuadFK(q0, DOF, IT_B, ET_F);        % Transformaci�n del espacio de configuraci�n al espacio de tarea   
Ok = T(1:3, 4);                         % Posici�n actual del efector 
Rk = T(1:3, 1:3);                       % Orientaci�n actual del efector
Qk = rot2cuat(Rk);                      % La orientaci�n se pasa a cuaterniones

% Inicializaci�n de variables
q = q0;
K = 500;
k = 0;

% Se obtiene la "diferencia" entre la rotaci�n deseada y la actual
Qe = multcuat(Qd, invcuat(Qk));
eo = Qe(2:4, 1);
ep = Od - Ok;

% C�lculo iterativo de la cinem�tica inversa utilizando m�todos num�ricos
while((norm(ep) > tol_p) && (solop == 1 || norm(eo) > tol_o) && (k < K))

    T = cuadFK(q, DOF, IT_B, ET_F);     % Transformaci�n homogenea de C a T
    Ok = T(1:3, 4);                     % Posici�n actual del efector final
    Rk = T(1:3, 1:3);                   % Orientacion actual del efector final

    Qk = rot2cuat(Rk);                  % Se obtiene la representaci�n de cuaternion de la matriz de rotaci�n
    Qe = multcuat(Qd, invcuat(Qk));     % Se multiplica Qd por Qk^-1

    eo = Qe(2:4, 1);
    ep = Od - Ok;                       % Error = Diferencia entre posici�n deseada y posici�n actual.
    e = [ep ; eo];

    % Jacobiano obtenido num�ricamente analizando la cinem�tica del
    % robot. Recordar que el jacobiano tiene 6 filas: 3 para la
    % posici�n (X, Y y Z) y 3 para los �ngulos (Roll, Pitch y Yaw). De
    % la misma manera, este cuenta con tantas columnas como actuadores
    % en el robot.
    J = cuadJ(q, DOF, IT_B, ET_F); 

    if solop == 1
        J = J(1:3, :);
        e = ep;
    end
    
    % C�lculo de la expresi�n que sustituye a la inversa del jacobiano
    switch metnum
        case 'pseudoinversa'
            if size(J,1) < size(J,2)
                Ji = J' * inv(J * J');
            else
                Ji = inv(J' * J) * J';
            end

        case 'traspuesta'
            alpha = (ep' * (J * J') * e) / (e' * (J * J') * (J * J') * e);
            Ji = alpha * J';

        case 'levenberg'
            Lambda = 0.1;
            I = eye(size(J,1));
            Ji = J' * inv(J * J' + Lambda*I);
    end

    % Matriz que multiplica la retroalimentaci�n por variables de estado
    Kp = eye(size(J,1)); 
    q = q + Ji * Kp * e;                % Se actualiza la soluci�n
    k = k + 1;                          % Se incrementa el n�mero de iteraci�n

end
    
end         