function [aging_channel,bessel_coef] = channel_aging(agg_channel,R_channel,velocity,time,carrier_freq)

    coef = zeros(1,time);
    c = 3*power(10,8);
    f_D = (velocity*carrier_freq)/c;
    sampleT = 10^-4;

    for n = 1:time
        coef(n) = 2*pi*f_D*sampleT*n;
    end

    bessel_coef = besselj(0,coef);
    bessel_coef_bar = sqrt(ones(size(bessel_coef))-bessel_coef.^2);

    % Aging channel 只限於R zone users (agg_channel後半部)

    aging_channel = zeros(length(agg_channel),width(agg_channel),time);


    for n = 1:time
        for k = 1:4
            if(k<3)
                aging_channel(:,k,n) = agg_channel(:,k);
            else
                [error] = generate_vector(length(R_channel(:,:,k)),1,0,R_channel(:,:,k));
                aging_channel(:,k,n) = bessel_coef(n) * agg_channel(:,k) + bessel_coef_bar(n) * error; 

            end
        end
    end

end