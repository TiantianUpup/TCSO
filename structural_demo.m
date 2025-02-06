%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used for test subsection 7.3 and subsection 7.4.

% If you want to test subsection 7.3, please set regular = 1, structural = 1, and mixed_sample=0.
% If you want to subsection 7.4, please set regular = 0, structural = 1, and mixed_sample=1.
% You can choose structure mask via mask_id.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
addpath(genpath(cd));

regular = 0;    % test image with regularization
structural = 1; % test structural missing (i.e., subsection 7.3)
mixed_sample=1; % test mixed missing (i.e., subsection 7.4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Our methods
EN_TCSONN = 1;  % TCSO-NN based on nuclear norm
EN_TCSOWNN = 1; % TCSO-WNN based on weighted nuclear norm

image_list = {'airplane.bmp', 'baboon.bmp', 'barbara.bmp', 'sailboat.bmp'};
img_num = 3; % image number, 1-4
Image = im2double(imread(image_list{img_num}));
%% sampling parameters
args.structural=structural;
mask_id= 5; % 1-5 is arranged by test order
args.mask_id=mask_id;

if mixed_sample
    args.sr=1; % 0.3, 0.4, 1
    args.structural = 1; % structural sampling
    [Z1, Omega1, X1, omega1] = sample(Image, args);
    args.structural=0; % random sampling
    [Z2, Omega2, X2, omega2] = sample(Image, args);
    Omega=intersect(Omega1,Omega2);
    [n1,n2,n3]=size(Z1);
    Z=zeros(size(Z1));
    Z(Omega) = Image(Omega);

    X = zeros(3 * n1, 3 * n2);
    X(1:n1, n2 + 1:2 * n2) = Z(:, :, 3);
    X(1:n1, 2 * n2 + 1:3 * n2) = -Z(:, :, 2);
    X(n1 + 1:2 * n1, 1:n2) = -Z(:, :, 3);
    X(n1 + 1:2 * n1, 2 * n2 + 1:3 * n2) = Z(:, :, 1);
    X(2 * n1 + 1:3 * n1, 1:n2) = Z(:, :, 2);
    X(2 * n1 + 1:3 * n1, n2 + 1:2 * n2) = -Z(:, :, 1);
    omega = find(X);
    imshow(Z,'border','tight');
    figure
else    
    [Z, Omega, X, omega] = sample(Image, args);
    imshow(Z,'border','tight');
    figure
end

%%%%%%%%%%%%%% Our methods
opts.max_iter = 150;
opts.max_beta = 1e10;
opts.tol = 1e-4;
opts.rho1 = 0.001;
opts.rho2 = 0.001;
% set iteration number for the PAMM subproblem
opts.max_inner_iter = 1; 
   
% TCSO-NN
if EN_TCSONN
    opts.kappa = 1.15;
    [m, n] = size(X);
    opts.beta = 30 / sqrt(m * n); %[10, 80] choose 30 (ignore the error in paper 30/sqrt(9mn))
    opts.regular = regular;

    % with regularization
    if regular
        opts.lam = 0.2;
        [m, n] = size(X);
        opts.beta = 30 / sqrt(m * n); %[10, 80] choose 30
    else
        % without regularization
        opts.beta = 5e-3;
        opts.lam = 0; %0.1
    end

    runhiststcso = TCSO_NN_test(Image, X, omega, opts);
    cput = runhiststcso.cput;
    iter = runhiststcso.iter;
    psnr = runhiststcso.psnr;
    ssim = runhiststcso.ssim;
    rse = runhiststcso.rse;
    Xtcsonn = runhiststcso.X;
    imshow(Xtcsonn, 'border', 'tight');
    figure
    fprintf("TCSO-NN iter is %d, cost is %3.2f, psnr is %3.4f, rse is %1.4f, ssim is %1.4f\n", iter, cput, psnr, rse, ssim);

end

% TCSO-WNN
if EN_TCSOWNN
    opts.C = 1;
    opts.beta = 1e-2; 
    opts.kappa = 1.2;

    if regular
        opts.lam = 0.1;
    else
        % without regularization
        opts.beta = 1e-4;
        opts.kappa = 1.15;
        opts.lam = 0; %0.1
    end

    runhiststcsownn = TCSO_WNN_ALM_test(Image, X, omega, opts);
    cput = runhiststcsownn.cput;
    iter = runhiststcsownn.iter;
    psnr = runhiststcsownn.psnr;
    ssim = runhiststcsownn.ssim;
    rse = runhiststcsownn.rse;
    Xtcsownn = runhiststcsownn.X;
    imshow(Xtcsownn, 'border', 'tight');
    fprintf("TCSO-WNN-ALM iter is %d, cost is %3.2f, psnr is %3.4f, rse is %1.4f, ssim is %1.4f\n", iter, cput, psnr, rse, ssim);
end

