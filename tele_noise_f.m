clear;
% close all;
clc;
format long
tic;

dt = 0.01;
T_max = 1e2;
t = 0:dt:T_max;
nt = length(t);
pos = 1;
num = 1e3;

dw = 2*pi/T_max;
w = 0:dw:pi/dt;

tele_f_mean = zeros((nt+1)/2,1);

for n = 1:num

    tele = zeros(nt,1);
    tele(1) = 1;
    for i = 2:nt
        if rand < pos*dt
            tele(i) = -tele(i-1);
        else
            tele(i) = tele(i-1);
        end
    end

    tele_f = fft(tele);
    tele_f_abs = abs(tele_f);
    tele_f_half = tele_f_abs(1:(nt+1)/2);
    tele_f_mean = tele_f_mean + tele_f_half;
end

tele_f_mean = tele_f_mean/num;

figure
subplot(2,1,1)
plot(t,tele)
xlabel('t')
ylabel('noise')

subplot(2,1,2)
plot(w,tele_f_mean)
xlabel('w')
ylabel('noise_f')
axis([0,300,-inf,inf])

% figure
% semilogy(w,tele_f_mean)
% xlabel('w')
% ylabel('noise_f')

toc;