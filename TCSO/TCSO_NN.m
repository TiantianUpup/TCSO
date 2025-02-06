 function runhist=TCSO_NN(O, Omega,opts)
    % INPUT
    %====================================
    % O ........... missing matrix according to the missing tensor
    % Omega ....... index set
    % opts
    % OUTPUT
    % ====================================
    % Y ........... recovered matrix
    % iter ........ iteration number 
    
    iter=0;
    
    % init parameters
    tol = opts.tol;
    max_iter = opts.max_iter; 
    beta=opts.beta;
    max_beta = opts.max_beta; 
    kappa=opts.kappa;
    rho1=opts.rho1;
    rho2=opts.rho2;
    lam=opts.lam;
    regular=opts.regular;
    
    [n1, n2] = size(O);

    % Y and Z is auxiliary variable
    % Lam1 and Lam2 is multiplier
    % init optimited variable
    [X,Y,Z,Lam1,Lam2] =deal(zeros(n1, n2));
    
    res = 1;
    
    while (iter < max_iter) 
        iter=iter+1;

        %Y-subproblem
        Y = shrink((Lam1+beta*X+rho1*Y)/(beta+rho1), 1 / (beta+rho1));
    
        %Z-subproblem
        temp=(Lam2+beta*mirt_dctn(X)+rho1*Z)/(beta+rho1);
        Z=sign(temp).*max(abs(temp)-lam/(beta+rho1),0);
     
        %X-subproblem
        Xold=X;
        W=(mirt_idctn(beta*Z-Lam2)+rho2*Xold+beta*Y-Lam1)/(2*beta+rho2);
        X=proj_polyhedral(W);
        X(Omega)=O(Omega);
        
        res=norm(X(:)-Xold(:))/norm(Xold(:));
        
        %update Lagrange multiplier
        Lam1=Lam1+beta*(X-Y);
        Lam2=Lam2+beta*(mirt_dctn(X)-Z);

        % update penalty parameter
        if ~regular
            beta = min(kappa*beta,max_beta); 
        end
        
        % This is because in the first few iterations of the structural test, res is very small
        if iter > 30 && res < tol
            break;
        end  
    end
    
    runhist.Y = X;
    runhist.iter = iter;
end    