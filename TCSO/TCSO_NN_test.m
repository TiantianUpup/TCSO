function runhist = TCSO_NN_test(T, O, Omega,opts)
    timect = cputime;

    Image = T;
    [n1, n2, n3] = size(T);

    runhist = TCSO_NN(O, Omega, opts);
    Y = runhist.Y;
    % fold matrix to tensor
    Y1 = max(Y(n1 + 1:2 * n1, 2 * n2 + 1:3 * n2), -Y(n1 + 1:2 * n1, 2 * n2 + 1:3 * n2));
    Y2 = max(Y(1:n1, 2 * n2 + 1:3 * n2), -Y(1:n1, 2 * n2 + 1:3 * n2));
    Y3 = max(Y(1:n1, n2 + 1:2 * n2), -Y(1:n1, n2 + 1:2 * n2));
    Xhat = zeros(n1, n2, 3);
    Xhat(:, :, 1) = Y1;
    Xhat(:, :, 2) = Y2;
    Xhat(:, :, 3) = Y3;
    runhist.cput = cputime - timect;
    runhist.X=Xhat;
    maxP = max(abs(Image(:)));
    psnr = PSNR(Image, Xhat, maxP, Omega);
    runhist.psnr = psnr;
    runhist.ssim=ssim(T, Xhat);
    runhist.rse=RSE(Image, Xhat, maxP); 
end
