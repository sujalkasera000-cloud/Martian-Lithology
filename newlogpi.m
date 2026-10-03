function [logpi, dM] = newlogpi(x, lb, ub, d, s, H)
% MYLOGPI - Calculates the log-posterior probability 
% Integrates with myBerry (Berryman SC + Contact Cement)
% x is now a 7-parameter array.

n = length(x);
dM = newmyberri(x, H); % Call the updated forward model

% 1. Check if parameters are within bounds
% Note: using >= and <= so that bounds like P_cem = 0 or P_gas = 0 do not fail.
if sum(x >= lb) == n && sum(x <= ub) == n
    
    isValid = true; % Default flag for physical constraints
    
    % 2. Extract values based on H to check physical constraints
    % Constraints: Vp > Vs AND 2500 < Rho < 3100
    if sum(H) == 3 % H = [1 1 1] -> [Vp, Vs, Rho]
        vp = dM(1); vs = dM(2); rho = dM(3);
        if vp <= vs || rho <= 2500 || rho >= 3100
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
            if rho <= 2500 || rho >= 3100
                isValid = false;
            end
        elseif H(2)==1 && H(3)==1 % H = [0 1 1] -> [Vs, Rho]
            rho = dM(2);
            if rho <= 2500 || rho >= 3100
                isValid = false;
            end
        end
        
    elseif sum(H) == 1
        if H(3)==1 % H = [0 0 1] -> [Rho]
            rho = dM(1);
            if rho <= 2500 || rho >= 3100
                isValid = false;
            end
        end
    end
    
    % 3. Calculate Log Posterior if all constraints are satisfied
    if isValid
        % s.\ means element-wise left division: equivalent to (d - dM) ./ s
        logpi = -0.5 * norm(s.\(d - dM))^2; 
    else
        logpi = -inf;
    end
    
else
    % Fails the upper/lower bounds check
    logpi = -inf; 
end

end