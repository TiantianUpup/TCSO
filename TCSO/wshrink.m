function [P,sig] = wshrink(Y, kappa, w)
    % This method calculate the proximal mapping of weighted nuclear norm
    % INPUT
    % ====================================
    % Y ....... matrix
    % w ....... weighted vector
    % opts .... optional parameter       
    % OUTPUT
    % ====================================
    % P ....... proximal mapping
    % sig ..... weighed vector in the next interate
    
    % [U, Sig, V] = svd(X, 'econ');
    % sig=diag(Sig);
    % r=sum(sig~=0);
    % for i = 1:r
    %     eta=sig(i);
    %     c1=eta-epsilon;
    %     c2=(eta+epsilon)^2-2*kappa;
    %     %c2=(eta+epsilon)^2-4;
    %     if c2 < 0
    %         Sig(i, i)=0;
    %     else
    %         Sig(i, i) = (c1+sqrt(c2))/2;
    %     end        
    % end
    %fprintf("rank is %d\n", rank(Y));
    % [U, Sig, V] = lansvd(Y,450,'L');

    m=size(Y,1);
    n=size(Y,2);
  
    [U, Sig, V] = svd(Y, 'econ');
  
    %sig=diag(Sig);
   
    for i = 1:min(m,n)
        wi=w(i);  
        Sig(i, i) = max(Sig(i, i) - kappa*wi, 0);
    end

    % for i=1:150
    %     wi=w(i);
    %     Sig(i, i) = max(Sig(i, i) - kappa*wi, 0);
    % end   

    P = U * Sig * V';
    % [~, Sig, ~] = svd(P, 'econ');
    sig=diag(Sig);
end