clear;
clc;
format long
tic;

L = 500;
dt = 1e-3;
M = 1;
T = 0:M*dt:100;
nt = length(T);
num = 100;

k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
E_k = -2*cos(k');
% E_k = -2*(2*rand(length(k),1)-1);

D0 = 3;
Df = 3;

psi0 = D0;

phi0 = ones(L,1);
phi0 = phi0*sqrt(L/(2*sum(abs(phi0).^2)));
m = zeros(nt,3,num);
phi10 = phi0(1:L/2);
phi20 = phi0(L/2+1:L);
step = 100;
m0 = zeros(step,3);
m0(1,1) = 2*phi10'*phi20/L;
m0(1,2) = 0;
m0(1,3) = 1/2;

for i = 2:step
    %     E0 = 0;
    for j = 1:L/2
        H = [E_k(j)+m0(i-1,3) -2*m0(i-1,1)*psi0;-2*m0(i-1,1)*psi0 -E_k(j)-m0(i-1,3)];
        [V,D] = eig(H);
        phi10(j) = V(1,1);
        phi20(j) = V(2,1);
        %         E0 = E0 + V(:,1)'*H*V(:,1);
    end
    m0(i,1) = 2*phi10'*phi20/L;
    m0(i,2) = 0;
    m0(i,3) = (phi10'*phi10 - phi20'*phi20)/L;
    %     E0 = -2*m0(i).^2.*psi0*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
end
% phi = [phi1;phi2];

m_iGS = m0(end,:);

for n = 1:num
    phi1 = phi10;
    phi2 = phi20;
    m(1,:,n) = m_iGS;
    m_it = m_iGS;

    % Et = zeros(nt,1);
    % Et(1) = E0;

    count = 2;
    t_it = 0;
    for i = 2:nt*M
        t_it = t_it + dt;
        psif = Df*randn/sqrt(dt);
%         psif = Df*randn;
        
        a = E_k + m_it(3);
        b = 2*m_it(1)*psif;
        c = 2*m_it(2)*psif;

        fact = sqrt(E_k.^2+b^2+c^2);
        ft = fact*dt;
        ss = sin(ft);
        ss = ss./fact;
        cc = cos(ft);
        Es = E_k.*ss;
        bs = (b-c*1i)*ss;
        phi1n = (cc-1i*Es).*phi1 +1i*bs.*phi2;
        phi2 = (cc+1i*Es).*phi2 +1i*conj(bs).*phi1;
        phi1 = phi1n;
        m_it(1) = (phi1'*phi2 + phi2'*phi1)/L;
        m_it(2) = -1i*(phi1'*phi2 - phi2'*phi1)/L;
        m_it(3) = (phi1'*phi1 - phi2'*phi2)/L;

        if mod(i-1,M) == 0
            m(count,:,n) = real(m_it);
            %         Et(count) = Et(count) -2*m(count).^2*Vf*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
            %         phi_sum_store(count) = phi_sum -1;
            count = count + 1;
        end

    end
end

m_total = sqrt(sum(m.^2,2));

m_mean = mean(m_total,3);

fit_x = log(T(floor(nt/3):end));
fit_y = log(m_mean(floor(nt/3):end));

% m_mean2 = mean(abs(m),2);
% fit_y2 = log(m_mean2(floor(nt/3):end));

toc;

filename = strcat('L = ',num2str(L), ', D0 = ', num2str(D0), ', Df = ', num2str(Df));
figure('Name',filename);
set(gcf, 'position', [250 70 1500 900]);

subplot(1,2,1)
plot(T,m_mean)
xlabel('t')
ylabel('m')

subplot(1,2,2)
plot(log(T),log(m_mean))
xlabel('t')
ylabel('m')
