function psnr = PSNR(Xtrue, Xrecover, maxP, omega)
    % INPUT
    %====================================
    % Xtrue ...... true tensor
    % Xrecover ... recovered tensor
    % maxP ....... max pixel
    % OUTPUT
    %====================================
    % psnr

    Xrecover = max(0, Xrecover);
    Xrecover = min(maxP, Xrecover);
    [n1, n2, n3] = size(Xrecover);
    MSE = norm(Xtrue(:) - Xrecover(:))^2 / (n1 * n2 * n3);
    psnr = 10 * log10(maxP^2 / MSE);

    % Xrecover = max(0, Xrecover);
    % Xrecover = min(maxP, Xrecover);

    % dim = size(Xrecover);
    % omegac = setdiff(1:prod(dim), omega);
    % num_missing = length(omegac);

    % Xtemp = Xtrue - Xrecover;
    % %erec = norm(Xtemp(:))^2;
    % MSE = norm(Xtemp(omegac))^2 / num_missing;
    % psnr = 10 * log10(maxP^2 / MSE);

    % d1=max(Xrecover(:));
    % d2=max(Xtrue(:));
    % d=max(d2);
    % sigma=mean2((Xtrue-Xrecover).^2);

    % psnr=10*log10((d.^2)./sigma);
end
