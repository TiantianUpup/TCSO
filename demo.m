%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used for test subsection 7.2       
%                         
% The test for random missing is easy, you only need to set img_num and sr
% Choose running method via EN_TCSONN and EN_TCSOWNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
addpath(genpath(cd));

img_flag = 0;
structural = 0; % test structural missing (i.e., subsection 7.3)
regular = 1;    % test image with regularization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Our methods
EN_TCSONN = 1;  % TCSO-NN based on nuclear norm
EN_TCSOWNN = 0; % TCSO-WNN based on weighted nuclear norm

image_list = {'airplane.bmp', 'baboon.bmp', 'barbara.bmp', 'sailboat.bmp'};
img_num = 1; % image number, 1-4
Image = im2double(imread(image_list{img_num}));

% 0.05, 0.1, 0.2, 0.3, 0.4
for sr = [0.2]

    %% Sampling
    args.structural = structural;
    args.sr = sr;
    [Z, Omega, X, omega] = sample(Image, args);
    imshow(Z, 'border', 'tight');
    figure


    %%%%%%%%%%%%%% Our methods
    opts.max_iter = 150;
    opts.max_beta = 1e10;
    opts.tol = 1e-4;
    opts.rho1 = 0.001;
    opts.rho2 = 0.001;

    % set iteration number for the PAMM subproblem
    if ismember(sr, [0.05, 0.1, 0.2])
        opts.max_inner_iter = 5; % 5 for low SRs, 1 for high SRs (0.3, 0.4).
    else
        opts.max_inner_iter = 1; % 5 for low SRs (0.05, 0.1, 0.2), 1 for high SRs (0.3, 0.4).
    end

    % TCSO-NN
    if EN_TCSONN
        opts.kappa = 1.15;
        [m, n] = size(X);
        opts.beta = 30 / sqrt(m * n); %[10, 80] choose 30 (ignore the error in paper 30/sqrt(9mn))
        opts.regular = regular;
        opts.lam = 0.1;
            
        runhiststcso = TCSO_NN_test(Image, X, omega, opts);
        cput = runhiststcso.cput;
        iter = runhiststcso.iter;
        psnr = runhiststcso.psnr;
        ssim = runhiststcso.ssim;
        rse = runhiststcso.rse;
        Xtcsonn = runhiststcso.X;
        imshow(Xtcsonn, 'border', 'tight');
        figure
        fprintf("sr is %3.2f, TCSO-NN iter is %d, cost is %3.2f, psnr is %3.4f, rse is %1.4f, ssim is %1.4f\n", sr, iter, cput, psnr, rse, ssim);

    end

    % TCSO-WNN
    if EN_TCSOWNN
        opts.C = 1;
        opts.beta = 1e-2; 
        opts.kappa = 1.2;

        if ismember(sr, [0.05, 0.1])
            opts.lam = 0.1;
        end

        if ismember(sr, [0.2, 0.3, 0.4])
            opts.lam = 0.01;
        end

        runhiststcsownn = TCSO_WNN_ALM_test(Image, X, omega, opts);
        cput = runhiststcsownn.cput;
        iter = runhiststcsownn.iter;
        psnr = runhiststcsownn.psnr;
        ssim = runhiststcsownn.ssim;
        rse = runhiststcsownn.rse;
        Xtcsownn = runhiststcsownn.X;
        imshow(Xtcsownn, 'border', 'tight');
        fprintf("sr is %3.2f, TCSO-WNN-ALM iter is %d, cost is %3.2f, psnr is %3.4f, rse is %1.4f, ssim is %1.4f\n", sr, iter, cput, psnr, rse, ssim);
    end
end
