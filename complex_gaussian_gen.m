function [Z] = complex_gaussian_gen(M,N,mean_matrix,covariance_matrix)

    % M: 矩陣的行數
    % N: 矩陣的列數
    
    size = M * N;

    samples = mvnrnd(mean_matrix , covariance_matrix, size);

    Z = samples(:, 1) + 1i * samples(:, 2);
    
    Z = reshape(Z, M, N);

end