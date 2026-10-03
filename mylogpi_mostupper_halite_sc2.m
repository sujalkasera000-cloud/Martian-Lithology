function [logpi, dM] = mylogpi_mostupper_halite_sc2(x, lb, ub, d, s, H)
% MYLOGPI - Calculates the log-posterior probability 
% Integrates with myBerry (Berryman SC + Contact Cement)
% x is now a 7-parameter array.
n = length(x);

% 1. STRICT PRIORS CHECK (SHORT-CIRCUIT)
% Reject proposals outside bounds IMMEDIATELY before running models
if any(x < lb) || any(x > ub)
    logpi = -inf;
    dM = zeros(sum(H), 1); % Dummy output to satisfy function signature
    return; % Exit the function right now, saving computation time and preventing errors
end

% 2. Run the forward model (It is now safe from negative cement values)
dM = myberri_mostupper_halite_sc2(x, H); 
    
isValid = true; % Default flag for physical constraints
    
% 3. Extract values based on H to check physical constraints
% Constraints: Vp > Vs AND 2100 < Rho < 3100
if sum(H) == 3 % H = [1 1 1] -> [Vp, Vs, Rho]
    vp = dM(1); vs = dM(2); rho = dM(3);
    if vp <= vs || rho <= 2100 || rho >= 3100
        isValid = false;
    end
    
elseif sum(H) == 2
    if H(1)==1 && H(2)==1 % H = [1 1 0] -> [Vp, Vs]
        vp = dM(1); vs = dM(2);
        if vp <= vs
            isValid = false;
        end
    elseif H(1)==1 && H(3)==1 % H = [1 0 1] -> [Vp, Rho]
        rho = dM(2);
        if rho <= 2100 || rho >= 3100
            isValid = false;
        end
    elseif H(2)==1 && H(3)==1 % H = [0 1 1] -> [Vs, Rho]
        rho = dM(2);
        if rho <= 2100 || rho >= 3100
            isValid = false;
        end
    end
    
elseif sum(H) == 1
    if H(3)==1 % H = [0 0 1] -> [Rho]
        rho = dM(1);
        if rho <= 2100 || rho >= 3100
            isValid = false;
        end
    end
end
    
% 4. Calculate Log Posterior if all constraints are satisfied
if isValid
    % s.\ means element-wise left division: equivalent to (d - dM) ./ s
    logpi = -0.5 * norm(s.\(d - dM))^2; 
else
    logpi = -inf;
end

end