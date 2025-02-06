function w = cal_fixed_weight(Y, C)
    % This function calculate the weighted vector (fixed)
    % INPUT
    % ====================================
    % Y ........ initial Y
    % C ..... Constant
    % OUTPUT
    % ====================================
    % w .......... weighted vector
    
    [~, Sig, ~] = svd(Y, 'econ');    % m=n  svd(Y,"econ") equals to svd(Y)
    r=rank(Sig);                     % the rank of Y
    sig=diag(Sig);

    len=length(sig);
    w=zeros(len,1);
    for i=1:r
        w(i)=C/sig(i);
    end  
    
    for i=r+1:len
        w(i)=C/sig(r);
    end    
    % fprintf("weighted is \n");
    % disp(w);
    % fprintf("==============================================\n") 
end