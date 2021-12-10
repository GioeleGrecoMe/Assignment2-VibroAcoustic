clear all, close all, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           EXERCISE 2                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data & Matrices
M1 = 5;
M2 = 1.25;
M3 = 10;
J1 = 2.5;
J2 = 0.16;

k1 = 1000;
k2 = 100;
k3 = 560;
k4 = 800;

c1 = 0.5;
c2 = 0.5;
c3 = 1;
c4 = 4;

R1 = 1;
R2 = 0.5;

x3_0 = 0.1;
teta1_0 = pi/12;
teta2_0 = -pi/12;

v3_0=1;
om1_0=0.5;
om2_0=2;

%create the matrices
jac_M = [1 0 0;0 0 0; 1 0 0; 0 1 0; 1 0 0; 0 0 1];

M = diag([M3 0 M1 J1 M2 J2]);
M_gen = jac_M'*M*jac_M;


jac_K = [-1 0 0;-1 -R1 R2;-1 R1 0;1 0 0];

K = diag([k1 k2 k3 k4]);
K_gen = jac_K'*K*jac_K;


jac_C = jac_K;

C = diag([c1 c2 c3 c4]);
C_gen = jac_C'*C*jac_C;


%% UNDAMPED - Eigenfrequencies & eigenvectors

% Method: eigenvalues-eigenvectors problem 
% lambda = i*omega;
[V,D] = eig(inv(M_gen)*K_gen); % V are the eigenvectors, D are the eigenvalues
%%
w_nat = sqrt(diag(D)) %natural frequencies
 
%% Normalization
%1st normalization
V_1 = [V]./[V(1,:)];
%2nd normalization
V_2 = [V]./[V(2,:)];
%3rd normalization
V_3 = [V]./[V(3,:)];

%% Sorting the solutions
w_nat_ord = [w_nat(1); -w_nat(1); w_nat(2); -w_nat(2); w_nat(3); -w_nat(3)];
V_ord = [V_2(:,1), V_2(:,1), V_2(:,2), V_2(:,2),V_3(:,3), V_3(:,3)]; %all vibration modes

%% UNDAMPED - Free motion

%set the initial conditions
 pos_t0 = [x3_0;teta1_0;teta2_0]; %position
 vel_t0 = [v3_0;om1_0;om2_0]; %velocity
%%
%Generic initial conditions
% pos_t0 = [2;pi/12]; %position
% vel_t0 = [8;-5]; %velocity

%1st mode
% pos_t0 = [V_1(:,1)]; %position
% vel_t0 = [0;0]; %velocity

%2nd mode
%pos_t0 = [V_1(:,2)]; %position
%vel_t0 = [0;0]; %velocity

%% Analytical solution
const = inv([V_ord;1i*w_nat_ord'.*V_ord])*[pos_t0;vel_t0];

t = 0:2*pi/10/max(w_nat_ord):10;
free_motion = zeros(3,length(t));

for ii = 1:length(t)
    free_motion(:,ii) = (const.'.*V_ord)*exp(1i*w_nat_ord*t(ii)); %const.' = transpose of the non-conjugate
end

%% Integration of Equation Of Motion 
% in [A] the submatrix [C] is [0]
A = - inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*[zeros(size(M_gen)) K_gen; -M_gen zeros(size(M_gen))];
odefun = @(t,y) A*y;
[t_ode,sol] = ode45(odefun,[t(1) t(end)],[vel_t0;pos_t0]);  %solve differential equation (function, time span, initial conditions)

%% Visualization
figure(1)
hold on
plot(t,free_motion(1,:),'r',t,free_motion(2,:)*R1,'b',t,free_motion(3,:)*R2,'g');
title('Undamped Time Responce')
%plot(t_ode,sol(:,4),'or',t_ode,sol(:,5)*R1,'ob',t_ode,sol(:,6)*R1,'--');
grid minor
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('x - analytical solution', '\theta*R1 - analytical solution', '\theta*R2 - analytical solution');% 'x - integrated solution', 'theta*R - integrated solution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DAMPED - Eigenfrequencies & eigenvectors

% in [A_damp] the submatrix [C] is full
% use state form matrix: new variable z = [vel; pos];
A_damp = - inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*[C_gen K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp);
V_damp = V_damp(4:6,:); %consider only the displacements
lambda = diag(D_damp); %complex and conjugate

%% DAMPED - Free motion

% pos_t0 = [2;pi/12];
% vel_t0 = [10;5];

%test initial condition as mode shape
pos_t0 = real(V_damp(:,3));
vel_t0 = real(lambda(3)*V_damp(:,3));

%% Analytical solution
const = [V_damp;lambda.'.*V_damp]\[pos_t0;vel_t0];

t = 0:2*pi/20/max(imag(lambda)):30;
free_motion = zeros(3,length(t));

for ii = 1:length(t)
    free_motion(:,ii) = (const.'.*V_damp)*exp(lambda*t(ii));
end

%% Integration of EOM
odefun = @(t,y) A_damp*y;
[t_ode,sol] = ode45(odefun,[t(1) t(end)],[vel_t0;pos_t0]);
%[t_ode,sol] = ode45(odefun,t,[vel_t0;pos_t0]);

figure(2)
hold on
plot(t,free_motion(1,:),'r',t,free_motion(2,:)*R1,'b',t,free_motion(3,:)*R2,'g');
title('Undamped Time Responce')
%plot(t_ode,sol(:,4),'or',t_ode,sol(:,5)*R1,'ob',t_ode,sol(:,6)*R1,'--');
grid minor
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('x - analytical solution', '\theta*R1 - analytical solution', '\theta*R2 - analytical solution');% 'x - integrated solution', 'theta*R - integrated solution')

% % %% State form of undamped motion
% % % z = {vel; pos};
% % [VV,DD] = eig(A);
% % DD_ord = diag(DD);
% % VV_2 = [VV]./[VV(4,:)];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Raylight damping
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1)C alpha=0.302 beta=0.0014
errorStored=1000000;
c=0;
for i=1:200
    alfa=0.302+i/100000;
    for j=1:200
        c=c+1;
        c/(200*200)
        beta=0.0009+j/100000;
         C_ray=alfa*M_gen+beta*K_gen;
         error=rms(rms(C_gen-C_ray));
         errorPerc=rms(rms((C_gen-C_ray)./(C_gen)))*100;
         if(errorPerc<errorStored)
             ALPBETA=[alfa,beta]
             errorStored=errorPerc;
         end
    end
    
end
%%
% w_quadro=[[1;1;1], w_nat.^2];
% EIG=eig(A_damp);
% %h=[1-(imag(EIG(1))/w_nat(1))^2;1-(imag(EIG(2))/w_nat(2))^2;1-(imag(EIG(3))/w_nat(3))^2];
% h=xi;
% b=h.*(w_nat.^2);
% ALPBETA=pinv(w_quadro)*b; 

Alp=ALPBETA(1);
Beta=ALPBETA(2);
C_ray=Alp.*M_gen+Beta.*K_gen;
error=rms(rms(C_gen-C_ray));
errorPerc=rms(rms((C_gen-C_ray)./(C_gen)))*100
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            Forced motion
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DAMPED - FREQUENCY RESPONSE OF THE SYSTEM

w=0:0.01:20;
H_tot=zeros(1,9);
for ii=1:length(w)
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_ray+K_gen);
    H11(ii)=H(1,1); H12(ii)=H(1,2); H21(ii)=H(2,1); H22(ii)=H(2,2);
    H31(ii)=H(3,1); H32(ii)=H(3,2); H33(ii)=H(3,3); H13(ii)=H(1,3);
    H23(ii)=H(2,3);
    H_tot(ii,:)=reshape(H,1,9);
%     for i=1:3
%         for j=1:3
%             if(i==1 && j==1)
%                 H_tot(ii,1)=H(1,1);
%             else
%                 H_tot(ii,:)=reshape(;
%             end
%         end
%     end
end

% transfer function between x,theta (lines) and q1, q2 (columns)
figure(3)
sgtitle('Frequency Responce Function')
for i=1:size(H_tot,2)
    subplot(3,3,i)

    grid on; 
     xlabel('Frequency [Hz]');
    if (i<4) 
            plot(w/2/pi,abs(H_tot(:,i)),'b'); 
        ylabel('|H_{x_3}| [m/N]'); 
        title(['|H_{x 3}| ',num2str(i)]);
    elseif(i<7)
            plot(w/2/pi,abs(H_tot(:,i))*R1,'b'); 
        ylabel('|H_{\theta_1}| [m/N]'); 
        title(['|H_{\theta 1}| ',num2str(i)]);
    else
            plot(w/2/pi,abs(H_tot(:,i))*R2,'b'); 
        ylabel('|H_{\theta_2}| [m/N]'); 
        title(['|H_{\theta 2}| ',num2str(i)]);
    end
end

figure(4)
sgtitle('Frequency Responce Function')
for i=1:size(H_tot,2)
    subplot(3,3,i)
    plot(w/2/pi,angle(H_tot(:,i)),'b');
    grid on; 
     xlabel('Frequency [Hz]');
    if (i<4) 
        ylabel('\angle H_{x_3} [m/N]'); 
        title(['\angle H_{x 3} ',num2str(i)]);
    elseif(i<7)
        ylabel('\angle H_{\theta_1} [m/N]'); 
        title(['\angle H_{\theta 1} ',num2str(i)]);
    else
        ylabel('\angle H_{\theta_2} [m/N]'); 
        title(['\angle H_{\theta 2} ',num2str(i)]);
    end
end
%% Plot the co-located FRF of point A at the centre of the disk (co-located meaning the FRF
%between the displacement in A and a force applied in A).

%xA=x3+R2*teta2
jac_A=[1;0;R2];

for ii=1:length(w)
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_ray+K_gen);
    FRF(ii)=pinv(jac_A)*H*jac_A;
end
plot(w/2/pi,abs(FRF),'b')

%% DAMPED - Free motion + Forced motion
A1=15;
A2=7;
f1=1.5;
f2=3.5;
O1=f1*2*pi;
O2=f2*2*pi;
%F_x=10;
%F_theta=100;
%F_vec=%[F_x; F_theta];

lambda_f=[1,0,R2];

w=0:0.01:20;
for ii=1:length(w)
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_ray+K_gen);
    F_vec=A1;
    XX=H*lambda_f'*F_vec;
    XX1(ii)=XX(1);
    XX2(ii)=XX(2);
    XX3(ii)=XX(3);
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_ray+K_gen);
    F_vec=A2;
    XX=H*lambda_f'*F_vec;
    XX1(ii)=XX1(ii)+XX(1);
    XX2(ii)=XX2(ii)+XX(2);
    XX3(ii)=XX3(ii)+XX(3);
end

%% FRF between x,theta and F1,F2
figure(4)
subplot(2,2,1); plot(w/2/pi,abs(XX1),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_x| [m/N]'); title('FRF');
subplot(2,2,3); plot(w/2/pi,angle(XX1),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('phase(H_x) [rad]')
subplot(2,2,2); plot(w/2/pi,abs(XX2),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_\theta| [m/N]'); title('FRF');
subplot(2,2,4); plot(w/2/pi,angle(XX2),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('phase(H_\theta) [rad]')

%% Time histories
W=(O1);%imag(lambda(1)));
% W = 3;
index=find(abs(w-W)==min(abs(w-W)));

t=0:2*pi/20/max(imag(lambda)):30;
X=real([XX1(index);XX2(index)]*exp(1i*W*t));

figure(5)
sgtitle('Time Responce')
ax1=subplot(2,1,1); plot(t,X(1,:)+free_motion(1,:),'r--',t,X(1,:),'r'); grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('x - free + forced motion','x - forced motion')
ax2=subplot(2,1,2); plot(t,X(2,:)+free_motion(2,:),'b--',t,X(2,:),'b'); grid on
xlabel('Time [s]')
ylabel('Amplitude [rad]')
legend('\theta - free + forced motion','\theta - forced motion')
linkaxes([ax1 ax2],'x')

