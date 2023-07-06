clear;
clc;
format long
tic;

L = 1000;
dt = 0.001;
M = 1;
T_max = 100;
T = 0:M*dt:T_max;
nt = length(T);
nt_real = round(T_max/dt)+1;
num = 2;

k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
E_k = -2*cos(k');
% E_k = -2*(2*rand(length(k),1)-1);
% E_k = ones(L/2,1);

D0 = 3;
Df = 3;

psi0 = D0;

phi0 = ones(L,1);
phi0 = phi0*sqrt(L/(2*sum(abs(phi0).^2)));

phi10 = phi0(1:L/2);
phi20 = phi0(L/2+1:L);
step = 100;
m0 = zeros(step,1);
m0(1) = 2*phi10'*phi20/L;

for i = 2:step
    %     E0 = 0;
    for j = 1:L/2
        H = [E_k(j) -2*m0(i-1)*psi0;-2*m0(i-1)*psi0 -E_k(j)];
        [V,D] = eig(H);
        phi10(j) = V(1,1);
        phi20(j) = V(2,1);
        %         E0 = E0 + V(:,1)'*H*V(:,1);
    end
    m0(i) = 2*phi10'*phi20/L;
    %     E0 = -2*m0(i).^2.*psi0*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
end
% phi = [phi1;phi2];

m_iGS = m0(end);

m_collect = zeros(nt,num);
for n = 1:num
    m = zeros(nt,1);
    phi1 = phi10;
    phi2 = phi20;
    m(1) = m_iGS;
    m_it = m_iGS;

    % Et = zeros(nt,1);
    % Et(1) = E0;

    count = 2;
    t_it = 0;
    
    noise0 = randn(nt_real,1);
    noise0_f = fft(noise0);
    noise0_f_half = noise0_f(2:ceil(nt_real/2));
    dw = 2*pi/T_max;
    w_half = 2*pi/T_max:dw:pi/dt;
    noise_f_half = noise0_f_half./(w_half.^3/2)';
    noise_f = [0;noise_f_half; conj(flip(noise_f_half))];
    noise = ifft(noise_f);
    noise = noise*sqrt(sum(noise0.^2)/sum(noise.^2));

    for i = 2:nt_real
        t_it = t_it + dt;
        psif = Df*noise(i);
%         psif = Df*noise(i)/sqrt(dt);
%         psif = Df*randn;

        [phi1, phi2] = Heun_step(phi1, phi2, m_it, E_k, psif, dt, L);
        m_it = (phi1'*phi2 + phi2'*phi1)/L;
        norm = sqrt(abs(phi1).^2 + abs(phi2).^2);
        phi1 = phi1./norm;
        phi2 = phi2./norm;

        if mod(i-1,M) == 0
            m(count) = real(m_it);
            %         Et(count) = Et(count) -2*m(count).^2*Vf*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
            %         phi_sum_store(count) = phi_sum -1;
            count = count + 1;
        end

    end
    m_collect(:,n) = m;
end

m_abs = sqrt(mean(m_collect.^2,2));

fit_x = log(T(floor(nt/3):end));
fit_y = log(m_abs(floor(nt/3):end));
fit_result = fit(fit_x',fit_y,'poly1')

% m_mean2 = mean(abs(m),2);
% fit_y2 = log(m_mean2(floor(nt/3):end));

toc;

filename = strcat('L = ',num2str(L), ', D0 = ', num2str(D0), ', Df = ', num2str(Df));
figure('Name',filename);
set(gcf, 'position', [100 70 1700 900]);

subplot(1,3,1)
plot(T,m_abs)
xlabel('t')
ylabel('m')

subplot(1,3,2)
plot(log(T(floor(nt*0.05):end)),log(m_abs(floor(nt*0.05):end)))
xlabel('t')
ylabel('m')

subplot(1,3,3)
plot(T(floor(nt*0.05):end),log(m_abs(floor(nt*0.05):end)))
xlabel('t')
ylabel('m')


function [y1, y2] = Heun_step(phi1, phi2, m_it, E_k, psif,dt,L)
    b = 2*m_it*psif;

    fact1t = E_k.*phi1 + b*phi2;
    fact2t = -E_k.*phi2 + b*phi1;
    phi1z = phi1 - 1i*dt*fact1t;
    phi2z = phi2 - 1i*dt*fact2t;
    m_it = (phi1z'*phi2z + phi2z'*phi1z)/L;

    b = 2*m_it*psif;
    fact1 = E_k.*phi1z + b*phi2z;
    fact2 = -E_k.*phi2z + b*phi1z;
    y1 = phi1 - 1i*dt*(fact1t + fact1)/2;
    y2 = phi2 - 1i*dt*(fact2t + fact2)/2;

%     z = x + dt * fact1;
%     fact2 = f(z, field, Jz);
%     y = x + dt * (fact1 + fact2) / 2;
end