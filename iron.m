function [h_env_rt, H_env] = iron(h, theta, Q, a_id, b_id, F)
% IRON: numerically ironing of a 1D continuous function h
% over the interval [theta(a_id), theta(b_id)], weighted by
% probability density f (corresponding to distribution F)
%
% Inputs:
%   h : a vector representing function values (to be ironed)
%   theta : a vector of original domain
%   Q : range of quantile
%   a_id : left endpoint index of the interval
%   b_id : right endpoint index of the interval
%   F : distribution function on the domain
%
% Outputs:
%   h_env_rt: ironed version of h
%   H_env : concave envelope of the integral function, weighted by f
%
% References: Myerson (1981), Toikka (2011), Kang and Watt (2024a,b)

    % Step 0: Initialization

    if ~issorted(F)
        error('F must be invertible');
    end

    len = b_id - a_id + 1;

   
    F_hd = @(x) interp1(theta, F, x, 'linear', 'extrap'); % F handle

    F_inv = interp1(F, theta, Q); % F^{-1}

    F_inv_hd = @(x) interp1(Q, F_inv, x, 'spline', 'extrap'); % F inverse handle

    h_hd = @(y) interp1(theta(a_id:b_id), h(a_id:b_id), y, 'spline', 'extrap');

    h_comp = @(x) h_hd(F_inv_hd(x)); % h_comp = h(F_inv(q))


    Q_range = F_hd(theta(a_id:b_id)); % range of Q that [theta(a_id), theta(b_id)] maps to
    
    Q_grid = linspace(Q_range(1), Q_range(end),len);

    
   
    % Step 1: find the weighted integration of h
    y = h_comp(Q_grid);    % h(q)
    x = Q_grid;            % q
    stpsz = Q_grid(2)-Q_grid(1);

    integrand =y;
    H = cumtrapz(integrand)*stpsz;
    

    % Step 2: Compute the convex hull of the graph
    K = convhull(x, H);
    xHull = x(K);
    yHull = H(K);

    % Step 3: Sort hull points by x (increasing)
    [xSorted, sortIdx] = sort(xHull);
    ySorted = yHull(sortIdx);

    % Step 4: Exclude points that exhibit concavity
    xLower = [];
    yLower = [];
    for i = 1:length(xSorted)
        while length(xLower) >= 2
            % Calculate slopes of last two segments
            dx1 = xLower(end) - xLower(end-1);
            dy1 = yLower(end) - yLower(end-1);
            dx2 = xSorted(i) - xLower(end);
            dy2 = ySorted(i) - yLower(end);
            if dy1 * dx2 >= dy2 * dx1 % if locally concave, skip the point
                xLower(end) = [];
                yLower(end) = [];
            else
                break; % otherwise, evaluate the point as the conv hull
            end
        end
        xLower(end+1) = xSorted(i);
        yLower(end+1) = ySorted(i);
    end
    
    
    % Step 5: Interpolate the convex envelope over the original x grid
    H_env = interp1(xLower, yLower, Q_grid, 'linear', 'extrap');

    h_env = diff(H_env)/stpsz;
    h_env = [h_env h_env(end)];

    
    h_env_hd = @(y) interp1(Q_grid, h_env, y, 'linear', 'extrap');

    h_env_comp = @(x) h_env_hd(F_hd(x)); % weighted by F

    h_env_rt = h_env_comp(theta(a_id:b_id));

end
