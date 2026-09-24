function H = jointEntropyDiscrete(X)

    % X: N x D matrix of integer-valued observations
    % Each column can have a different range, e.g. 0:K_j
    %
    % H: empirical joint Shannon entropy in bits

    % Find unique joint configurations
    [~, ~, idx] = unique(X, 'rows');

    % Count occurrences of each configuration
    counts = accumarray(idx, 1);

    % Empirical joint probabilities
    p = counts / size(X, 1);

    % Shannon entropy
    H = -sum(p .* log2(p));

end