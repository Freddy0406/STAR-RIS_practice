function [sum_SE]=f_function(Rk,Psik,K,SNR_linear)


    sum_SE = 0;
    sum_inter1 = 0;
    sum_inter2 = 0;
    for i = 1:K
        signal = trace(Psik(:,:,i))^2;
        for j = 1:K           
            sum_inter1 = sum_inter1 + trace(Rk(:,:,i)*Psik(:,:,j));
            sum_inter2 = sum_inter2 + trace(Psik(:,:,j));
        end
        interference = sum_inter1 - trace((Psik(:,:,i))^2) + (K/SNR_linear) * sum_inter2;
        sum_SE = sum_SE+log2(1+abs(signal/interference));
    end

end