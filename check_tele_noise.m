clear;
clc;
format long
tic;

L = 2000;
dt = 0.001;
M = 100;
T_max = 1000;
T = 0:M*dt:T_max;
nt = length(T);
num = 100;
pos = 100;

k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
% E_k = -2*cos(k');
E_k = -2*(2*rand(length(k),1)-1);
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
parfor n = 1:num
    m = zeros(nt,1);
    phi1 = phi10;
    phi2 = phi20;
    m(1) = m_iGS;
    m_it = m_iGS;

    % Et = zeros(nt,1);
    % Et(1) = E0;

    count = 2;
    t_it = 0;
    tele = 1;

    for i = 2:round(T_max/dt)+1
        if rand < pos*dt
            tele = -tele;
        end

        t_it = t_it + dt;
        psif = Df*tele/sqrt(dt);
%         psif = Df*randn;

        b = 2*m_it*psif;

        fact = sqrt(E_k.^2+b^2);
        ft = fact*dt;
        ss = sin(ft);
        ss = ss./fact;
        cc = cos(ft);
        Es = E_k.*ss;
        bs = b*ss;
        phi1n = (cc-1i*Es).*phi1 +1i*bs.*phi2;
        phi2 = (cc+1i*Es).*phi2 +1i*bs.*phi1;
        phi1 = phi1n;
        m_it = (phi1'*phi2 + phi2'*phi1)/L;

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
