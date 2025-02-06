function S = shrink(X, tau)
    % singular value shrinkage operator of matrix
    %
    % INPUT
    % ====================================
    % X ....... matrix
    % tau ......threshold
    % OUTPUT
    % ====================================
    % S ....... shrinkage matrix
    
    [U, Sig, V] = svd(X, 'econ');
    for i = 1:size(Sig, 1)
        Sig(i, i) = max(Sig(i, i) - tau, 0);
    end

    S = U * Sig * V';
end
