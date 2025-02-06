function [Z, Omega, X, omega] = sample(T, opts)
    % INPUT
    %====================================
    % T ........... true tensor
    % opts
    % sr .......... sample rate
    % structural .. 1 for structural missing, 0 for random missing
    % mask_id ..... the index of mask image
    % OUTPUT
    % ====================================
    % Z ........... missing tensor
    % Omega ....... known index set of Z
    % X ........... missing tensor after flatting
    % omega ....... known index set of X

    path=fullfile(pwd,'mask');
    structural = opts.structural;
    
    [n1, n2, n3] = size(T);
    if structural
        mask_id=opts.mask_id;
        file_list = dir(fullfile(path, '*.bmp'));
        num_mask = length(file_list);
        mask_list = cell(num_mask, 1);
        for i = 1 : num_mask
            mask_list{i} = file_list(i).name; 
        end

        mask = im2double(imread(mask_list{mask_id}));
        Omega = find(mask);
    else    
        sr=opts.sr;
        Omega = find(rand(n1 * n2 * n3, 1) < sr);
    end    
 
    Z = zeros(n1, n2, n3);
    Z(Omega) = T(Omega);
    X = zeros(3 * n1, 3 * n2);
    X(1:n1, n2 + 1:2 * n2) = Z(:, :, 3);
    X(1:n1, 2 * n2 + 1:3 * n2) = -Z(:, :, 2);
    X(n1 + 1:2 * n1, 1:n2) = -Z(:, :, 3);
    X(n1 + 1:2 * n1, 2 * n2 + 1:3 * n2) = Z(:, :, 1);
    X(2 * n1 + 1:3 * n1, 1:n2) = Z(:, :, 2);
    X(2 * n1 + 1:3 * n1, n2 + 1:2 * n2) = -Z(:, :, 1);
    omega = find(X);
end