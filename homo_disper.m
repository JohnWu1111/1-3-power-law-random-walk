clear;
clc;
format long
tic;

dt = 1e-3;
M = 1;
T = 0:M*dt:100;
nt = length(T);
num = 100;

E_k = 1;

D0 = 3;
Df = 3;

psi0 = D0;

m = zeros(nt,1);
phi10 = 1/sqrt(2);
phi20 = 1/sqrt(2);
step = 100;
m0 = zeros(step,1);
m0(1) = phi10'*phi20;

for i = 2:step
    %     E0 = 0;
    H = [E_k -2*m0(i-1)*psi0;-2*m0(i-1)*psi0 -E_k];
    [V,D] = eig(H);
    phi10 = V(1,1);
    phi20 = V(2,1);
    %         E0 = E0 + V(:,1)'*H*V(:,1);
    m0(i) = phi10'*phi20;
    %     E0 = -2*m0(i).^2.*psi0*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
end
% phi = [phi1;phi2];

m_iGS = m0(end);

for n = 1:num
    phi1 = phi10;
    phi2 = phi20;
    m(1,n) = m_iGS;
    m_it = m_iGS;

    % Et = zeros(nt,1);
    % Et(1) = E0;

    count = 2;
    t_it = 0;
    for i = 2:nt*M
        t_it = t_it + dt;
        psif = Df*randn/sqrt(dt);
%         psif = Df*randn;

        b = 2*m_it*psif;

        fact = sqrt(E_k^2+b^2);
        ft = fact*dt;
        ss = sin(ft);
        ss = ss./fact;
        cc = cos(ft);
        Es = E_k*ss;
        bs = b*ss;
        phi1n = (cc-1i*Es).*phi1 +1i*bs.*phi2;
        phi2 = (cc+1i*Es).*phi2 +1i*bs.*phi1;
        phi1 = phi1n;
        m_it = (phi1'*phi2 + phi2'*phi1)/2;

        if mod(i-1,M) == 0
            m(count,n) = real(m_it);
            %         Et(count) = Et(count) -2*m(count).^2*Vf*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
            %         phi_sum_store(count) = phi_sum -1;
            count = count + 1;
        end

    end
end

m_mean = sqrt(mean(m.^2,2));
m_abs = m_mean;

fit_x = log(T(floor(nt/3):end));
fit_y = log(m_abs(floor(nt/3):end));

m_mean2 = mean(abs(m),2);
fit_y2 = log(m_mean2(floor(nt/3):end));

toc;

filename = strcat(', D0 = ', num2str(D0), ', Df = ', num2str(Df));
figure('Name',filename);
set(gcf, 'position', [250 70 1500 900]);

subplot(1,2,1)
plot(T,m_abs)
xlabel('t')
ylabel('m')

subplot(1,2,2)
plot(log(T),log(m_abs))
xlabel('t')
ylabel('m')
