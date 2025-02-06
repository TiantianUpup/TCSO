function runhist = TCSO_WNN_ALM(O, Omega, opts)
   
    % init parameters
    tol = opts.tol;
    max_iter = opts.max_iter;
    beta = opts.beta;
    max_beta = opts.max_beta;
    kappa = opts.kappa;
    rho1 = opts.rho1;
    rho2 = opts.rho2;
    lam = opts.lam;
    max_inner_iter = opts.max_inner_iter;


    [n1, n2] = size(O);

    % Y and Z are auxiliary variables
    % Lam1 and Lam2 are multipliers
    [X,Y,Z]=deal(O);
    [Lam1, Lam2] = deal(zeros(n1, n2));
    res = 1;
    iter = 0;
    opt_tol = 0.998;

    C=opts.C;
    % calculate initial weighted vector
    w=cal_fixed_weight(Y, C);
        
    while (iter < max_iter)
        iter = iter + 1; 
        inner_iter = 0;
        opt = 1;
    
        while opt > opt_tol && inner_iter < max_inner_iter
            inner_iter = inner_iter + 1;

            % Y-subproblem
            Yold = Y;
            [Y, sig]=wshrink((Lam1 + beta * X + rho1 * Y) / (beta + rho1), 1 / (beta + rho1), w);

            % Z-subproblem
            Zold = Z;
            temp = (Lam2 + beta * mirt_dctn(X) + rho1 * Z) / (beta + rho1);
            Z = sign(temp) .* max(abs(temp) - lam / (beta + rho1), 0);

            % X-subproblem
            Xold = X;
            W = (mirt_idctn(beta * Z - Lam2) + rho2 * Xold + beta * Y - Lam1) / (2 * beta + rho2);
            X = proj_polyhedral(W);
            X(Omega) = O(Omega);

            Xtemp = Xold - X;

            if mode(inner_iter, 2) == 1
                opt = sqrt(norm(rho2 * Xtemp) ^ 2 + norm(beta * Xtemp + rho1 * (Yold - Y)) ^ 2 + norm(beta * (mirt_dctn(Xtemp)) + rho1 * (Zold - Z)) ^ 2);
            end
            
        end


        opt_tol = max (0.998 * opt_tol, 0.1);

        % update Lagrange multiplier
        Lam1 = Lam1 + beta * (X - Y);
        Lam2 = Lam2 + beta * (mirt_dctn(X) - Z);

        % update penalty parameter
        beta = min(kappa * beta, max_beta);

        res = norm(X(:) - Xold(:)) / norm(Xold(:));
        
        if iter > 20 && res < tol
            break;
        end

    end

    runhist.Y = X;
    runhist.iter = iter;
end

