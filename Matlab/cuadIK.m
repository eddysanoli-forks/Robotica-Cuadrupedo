function [q, k] = cuadIK(Td, q0, DOF, IT_B, ET_F, varargin)
% CUADIK() Función que obtiene la cinemática inversa del cuadrúpedo, ya sea
% para las patas traseras (2 DOF) o las patas delanteras (3 DOF). La
% cinemática inversa permite mapear del espacio de tarea al espacio de
% configuración, por lo que se le alimenta una pose deseada y su config.
% inicial y este es capaz de computar la siguiente config.
% -------------------------------------------------------------------------
%
% Según el número de inputs adicionales se determina el modo
%   - 0 inputs = IK de posición usando Levenberg
%   - 1 input = IK de posición / posición-orientación usando Levenberg
%   - 2 inputs = IK de posición / posición-orientación usando "levenberg",
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

tol_p = 1e-06;                          % Tolerancia del error de posición
tol_o = 1e-05;                          % Tolerancia del error de orientación

Od = Td(1:3, 4);                        % Se extrae la posición y orientación deseada de Td
Rd = Td(1:3, 1:3);
Qd = rot2cuat(Rd);
       
T = cuadFK(q0, DOF, IT_B, ET_F);        % Transformación del espacio de configuración al espacio de tarea   
Ok = T(1:3, 4);                         % Posición actual del efector 
Rk = T(1:3, 1:3);                       % Orientación actual del efector
Qk = rot2cuat(Rk);                      % La orientación se pasa a cuaterniones

% Inicialización de variables
q = q0;
K = 500;
k = 0;

% Se obtiene la "diferencia" entre la rotación deseada y la actual
Qe = multcuat(Qd, invcuat(Qk));
eo = Qe(2:4, 1);
ep = Od - Ok;

% Cálculo iterativo de la cinemática inversa utilizando métodos numéricos
while((norm(ep) > tol_p) && (solop == 1 || norm(eo) > tol_o) && (k < K))

    T = cuadFK(q, DOF, IT_B, ET_F);     % Transformación homogenea de C a T
    Ok = T(1:3, 4);                     % Posición actual del efector final
    Rk = T(1:3, 1:3);                   % Orientacion actual del efector final

    Qk = rot2cuat(Rk);                  % Se obtiene la representación de cuaternion de la matriz de rotación
    Qe = multcuat(Qd, invcuat(Qk));     % Se multiplica Qd por Qk^-1

    eo = Qe(2:4, 1);
    ep = Od - Ok;                       % Error = Diferencia entre posición deseada y posición actual.
    e = [ep ; eo];

    % Jacobiano obtenido numéricamente analizando la cinemática del
    % robot. Recordar que el jacobiano tiene 6 filas: 3 para la
    % posición (X, Y y Z) y 3 para los ángulos (Roll, Pitch y Yaw). De
    % la misma manera, este cuenta con tantas columnas como actuadores
    % en el robot.
    J = cuadJ(q, DOF, IT_B, ET_F); 

    if solop == 1
        J = J(1:3, :);
        e = ep;
    end
    
    % Cálculo de la expresión que sustituye a la inversa del jacobiano
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

    % Matriz que multiplica la retroalimentación por variables de estado
    Kp = eye(size(J,1)); 
    q = q + Ji * Kp * e;                % Se actualiza la solución
    k = k + 1;                          % Se incrementa el número de iteración

end
    
end         