clc;                                                                                                                                                                                                                                                                  clc;
clear;
close all;


SNR_arr = [-10,0,5,10,15,20,25,30,35,40,45,50];
SE_arr = zeros(1,length(SNR_arr));
iteration  = 20;
% Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);

parfor times = 1:iteration
    temp_SE = zeros(1,length(SNR_arr));
    for i = 1:length(SNR_arr)
        [SE] = pgam_main_SNR_vs_SE(SNR_arr(i));
        temp_SE (i) =   temp_SE (i) + SE;
    end
    SE_arr = SE_arr + temp_SE;
    fprintf('%d complete\n',times);
    %Progressbar
    pause(100/iteration);
    ppm.increment();

end

%Delete Progressbar
delete(ppm);

SE_arr = SE_arr ./ iteration;

figure(1)
plot(SNR_arr,SE_arr,'--o');
grid on;
colororder([1 0 0]);
xlabel("SNR (dB)")
ylabel("Downlink Achievable Sum SE [bit/s/Hz]")
legend('ES protocol (N=64)',"Location","Best")