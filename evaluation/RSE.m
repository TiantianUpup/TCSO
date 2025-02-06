function rse = RSE(Xtrue, Xrecover, maxP)
    Xrecover = max(0, Xrecover);
    Xrecover = min(maxP, Xrecover);

    rse = norm(Xtrue(:) - Xrecover(:)) / norm(Xtrue(:));
end    