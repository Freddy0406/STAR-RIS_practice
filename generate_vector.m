function[out] = generate_vector(M,N,mean,cov_matrix)

    [V,D] = eig(cov_matrix);

    z = (1/sqrt(2))*(randn(M*N,1)+1i*randn(M*N,1));

    out = mean+V*sqrtm(D)*z;

end