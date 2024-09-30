function [new_Rk,new_Qk,new_Psik] = new_covariance(pathloss_d,pathloss_hat,R_RIS,R_BS,new_Phi_t,new_Phi_r,noise_var,tao,pilot_P_linear)

    K = length(pathloss_d);
    
    for i = 1:K
        if(i<3)
            new_Rk(:,:,i) =pathloss_d(i)*R_BS + pathloss_hat(i)*...
            trace(R_RIS*new_Phi_t*R_RIS*new_Phi_t')*R_BS;
        else
            new_Rk(:,:,i) = pathloss_d(i)*R_BS + pathloss_hat(i)*...
            trace(R_RIS*new_Phi_r*R_RIS*new_Phi_r')*R_BS;
        end   
    end
        
    for k = 1:K
        new_Qk(:,:,k) = pinv(new_Rk(:,:,k)+(noise_var/tao*pilot_P_linear)*eye(length(new_Rk)));
        new_Psik(:,:,k) = new_Rk(:,:,k)*new_Qk(:,:,k)*new_Rk(:,:,k);
    end



end