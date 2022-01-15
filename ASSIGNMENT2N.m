clear all, close all, clc


%% Data & Matrices
R1 = 1;
R2 = 0.5;
%Mass parameters
M1 = 5;
J1 = 2.5;
M2 = 1.25;
J2 = 0.16;
M3 = 10;
%Spring parameters
k1 = 1000;
k2 = 100;
k3 = 560;
k4 = 800;
%Damper parameters
c1 = 0.5;
c2 = 0.5;
c3 = 1;
c4 = 4;
%Mass matrix
jac_M =[1    0   R2;
        0    1   0;
        1    0   0;
        0    0   1;
        1    0   0];
M = diag([M1 J1 M2 J2 M3]);
M_gen = jac_M'*M*jac_M;

%Elastic matrix
jac_K = [-1   0       -R2;
         -1   -R1     R2;
         -1   -R1      0;
          1     0      0];
K = diag([k1 k2 k3 k4]);
K_gen = jac_K'*K*jac_K;

%Damping matrix
jac_C = jac_K;
C = diag([c1 c2 c3 c4]);
K_gg=jac_C'*K*jac_C;
C_gen = jac_C'*C*jac_C;

%% 1A DA SCRIVERE A MANO

%% 1B Undamped and damped eigen frequencies and eigenvectors of the system
% UNDAMPED - Eigenfrequencies & eigenvectors

% Method: eigenvalues-eigenvectors problem 
% lambda = i*omega;
[V,D] = eig(inv(M_gen)*K_gen); % V are the eigenvectors, D are the eigenvalues
% DAMPED - Eigenfrequencies & eigenvectors
w0_nat=sqrt(diag(D));
% in [A_damp] the submatrix [C] is full
% use state form matrix: new variable z = [vel; pos];
A_damp = - inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*[C_gen K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp);

%% 1C Rayleigh damping [C]=alpha*[M]+beta*[K]
% ANALYTIC MODE
Cii=reshape(C_gen,size(C_gen,1)*size(C_gen,2),1);
Mii=reshape(M_gen,size(M_gen,1)*size(M_gen,2),1);
Kii=reshape(K_gen,size(K_gen,1)*size(K_gen,2),1);
M_K=[Mii,Kii];
AB=pinv(M_K)*Cii;
alpha=AB(1);
beta=AB(2);
C_R_Analytic=alpha*M_gen+beta*K_gen;
%% %%%%%%%%%%%%%%%%%%% 2 part %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Free motion of the system (considering the Rayleigh damping as in 1.c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 A Plot and comment the free motion of the system starting from the initial conditionsreported

%% DAMPED - Eigenfrequencies & eigenvectors without RAYLEIGH
% Add the matrix [C] to the matrix [A]
A_damp = -inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*...
[C_gen K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp);

% Mode shapes
%Guardando in V_damp noto che gli autovettori più vicini agli undamped sono
%quelli in posizione 1:3, ma dovrebbero essere quelli in seconda posizione.
%come mai?
V_damp = V_damp(4:6,:); %consider only the displacements
V_damp_1 = V_damp./V_damp(1,:);%1st normalization
V_damp_2 = V_damp./V_damp(2,:);%2nd normalization
V_damp_3 = V_damp./V_damp(3,:);%3rd normalization

V_damp_ord = [V_damp_1(:,1), V_damp_1(:,1),...
              V_damp_1(:,2), V_damp_1(:,2),...
              V_damp_1(:,3), V_damp_1(:,3)]; %all vibration modes

%Natural frequencies
lambda = diag(D_damp);
alpha = real(lambda);
w_d = imag(lambda);

% Adimensional damping coefficient
 csi = abs(alpha./w_d); %should be w0 not w_d
 csi_1 = abs(alpha./sqrt(w_d.^2+alpha.^2)); %sotto h=10% simili
 
 pos_t0=[0.1;pi/12;-pi/12];
 vel_t0=[1;0.5;2];


% Analytical
const_damp = [V_damp; lambda.'.*V_damp]\[pos_t0;vel_t0];
t_damp = 0:2*pi/50/max(w_d):20;

free_motion_damp = zeros(3,length(t_damp));

for ii = 1:length(t_damp)
    free_motion_damp(:,ii) = (const_damp.'.*V_damp)*exp(lambda*t_damp(ii));
end

%figure(2)
hold on
plot(t_damp,free_motion_damp(1,:),'b',t_damp,free_motion_damp(2,:)*R2,'k',t_damp,free_motion_damp(3,:)*R1,'r'); 
grid minor
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('x - analytical solution', '\theta_1 - analytical solution', '\theta_2 - analytical solution')
title('Damped Time History')
%% DAMPED - Eigenfrequencies & eigenvectors whith Rayleigh
%C_gen=C_R_Analytic;
% Add the matrix [C] to the matrix [A]
A_damp = -inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*...
[C_R_Analytic K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp);

% Mode shapes
V_damp = V_damp(1:3,:); %consider only the displacements

V_damp_1 = V_damp./V_damp(1,:);%1st normalization
V_damp_2 = V_damp./V_damp(2,:);%2nd normalization
V_damp_3 = V_damp./V_damp(3,:);%3rd normalization

V_damp_ord = [V_damp_1(:,1), V_damp_1(:,1),...
              V_damp_1(:,2), V_damp_1(:,2),...
              V_damp_1(:,3), V_damp_1(:,3)]; %all vibration modes

%Natural frequencies
lambda = diag(D_damp);
alpha = real(lambda);
w_d = imag(lambda);

% Adimensional damping coefficient
 csi = abs(alpha./w_d); %should be w0 not w_d
 csi_1 = abs(alpha./sqrt(w_d.^2+alpha.^2)); %sotto h=10% simili
%% DAMPED - Free motion
% Normalization
%1st normalization
V_1 = [V_damp]./[V_damp(1,:)];
%2nd normalization
V_2 = [V_damp]./[V_damp(2,:)];
%3rd normalization
V_3 = [V_damp]./[V_damp(3,:)];

%pos_t0=V_3(:,6);
%vel_t0=[0,0,0]';

% Analytical
const_damp = [V_damp; lambda.'.*V_damp]\[pos_t0;vel_t0];
%t_damp = 0:2*pi/50/max(w_d):100;

free_motion_damp1 = zeros(3,length(t_damp));

for ii = 1:length(t_damp)
    free_motion_damp1(:,ii) = (const_damp.'.*V_damp)*exp(lambda*t_damp(ii));
end

%Plot the 2 different solutions
%figure(2)
hold on
plot(t_damp,free_motion_damp(1,:),'b',t_damp,free_motion_damp(2,:)*R2,'k',t_damp,free_motion_damp(3,:)*R1,'r'); 

plot(t_damp,free_motion_damp1(1,:),t_damp,free_motion_damp1(2,:)*R2,t_damp,free_motion_damp1(3,:)*R1); 
grid minor
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('x_3', '\theta_1', '\theta_2')
title('Damped Time History')

%% DAMPED - FREQUENCY RESPONSE OF THE SYSTEM
C_ray=C_R_Analytic;
w=0:0.01:30;
H_tot=zeros(1,9);
for ii=1:length(w)
    H = inv((-w(ii)^2)*M_gen+1i*w(ii)*C_ray+K_gen);
    H11(ii)=H(1,1); H12(ii)=H(1,2); H21(ii)=H(2,1); H22(ii)=H(2,2);
    H31(ii)=H(3,1); H32(ii)=H(3,2); H33(ii)=H(3,3); H13(ii)=H(1,3);
    H23(ii)=H(2,3);
    H_tot(ii,:)=reshape(H,1,9);

end

% transfer function between x,theta (lines) and q1, q2 (columns)
figure(3)
sgtitle('Frequency Responce Function')
for i=1:size(H_tot,2)
    subplot(3,3,i)

    grid on
     xlabel('Frequency [Hz]');
    if (i<4) 
            semilogy(w/2/pi,abs(H_tot(:,i)),'b'); 
        ylabel(['|H 1',num2str(i),'| [m/N]']); 
        title(['|H 1',num2str(i),'|']);
    elseif(i<7)
            semilogy(w/2/pi,abs(H_tot(:,i)),'b'); 
        ylabel(['|H 2',num2str(i-3),'| [m/N]']); 
        title(['|H 2',num2str(i-3),'|']);
    else
            semilogy(w/2/pi,abs(H_tot(:,i)),'b'); 
        ylabel(['|H 3',num2str(i-6),'| [m/N]']); 
        title(['|H 3',num2str(i-6),'|']);
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
        ylabel(['\angle H1',num2str(i),' [m/N]']); 
        title(['\angle H1',num2str(i)]);
    elseif(i<7)
        ylabel(['\angle H2',num2str(i-3),' [m/N]']); 
        title(['\angle H2',num2str(i-3)]);
    else
        ylabel(['\angle H3',num2str(i-6),' [m/N]']); 
        title(['\angle H3',num2str(i-6)]);
    end
end

%% DAMPED - Frequency Responce Function (FRF) with normal Damping (NO REYLEIGH

%w = 0:0.01:30;

% [H(w)] = [D(w)]^-1 = [-w^2*[M]+i*w*[C]+[K]]^-1
for ii = 1:length(w)
    FRF = inv(-w(ii)^2*M_gen+1i*w(ii)*C_gen+K_gen);
    FRF11(ii) = FRF(1,1); FRF12(ii) = FRF(1,2); FRF13(ii) = FRF(1,3);
    FRF21(ii) = FRF(2,1); FRF22(ii) = FRF(2,2); FRF23(ii) = FRF(2,3);
    FRF31(ii) = FRF(3,1); FRF32(ii) = FRF(3,2); FRF33(ii) = FRF(3,3);
end

% Plot of the FRF modulus
figure(5)
subplot(3,3,1); plot(w/2/pi,abs(FRF11),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_1')
subplot(3,3,2); plot(w/2/pi,abs(FRF12),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_2| [m/N]'); title('H_1_2')
subplot(3,3,3); plot(w/2/pi,abs(FRF13),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_3| [m/N]'); title('H_1_3')
subplot(3,3,4); plot(w/2/pi,abs(FRF21),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_1| [m/N]'); title('H_2_1')
subplot(3,3,5); plot(w/2/pi,abs(FRF22),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_2| [m/N]'); title('H_2_2')
subplot(3,3,6); plot(w/2/pi,abs(FRF23),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_2_3| [m/N]'); title('H_2_3')
subplot(3,3,7); plot(w/2/pi,abs(FRF31),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_1| [m/N]'); title('H_3_1')
subplot(3,3,8); plot(w/2/pi,abs(FRF32),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_2| [m/N]'); title('H_3_2')
subplot(3,3,9); plot(w/2/pi,abs(FRF33),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [m/N]'); title('H_3_3')
sgtitle('Amplitude of FRF')

% Plot of the FRF phase
figure(6)
subplot(3,3,1); plot(w/2/pi,angle(FRF11)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_1 [deg]'); title('H_1_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,2); plot(w/2/pi,angle(FRF12)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_2 [deg]'); title('H_1_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,3); plot(w/2/pi,angle(FRF13)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_1_3 [deg]'); title('H_1_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,4); plot(w/2/pi,angle(FRF21)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_1 [deg]'); title('H_2_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,5); plot(w/2/pi,angle(FRF22)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_2 [deg]'); title('H_2_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,6); plot(w/2/pi,angle(FRF23)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_2_3 [deg]'); title('H_2_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,7); plot(w/2/pi,angle(FRF31)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_1 [deg]'); title('H_3_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,8); plot(w/2/pi,angle(FRF32)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_2 [deg]'); title('H_3_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,9); plot(w/2/pi,angle(FRF33)*180/pi,'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('H_3_3 [deg]'); title('H_3_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Phase of FRF')
%% 3B Plot the co-located FRF of point A at the centre of the disk (co-located meaning the FRF
%between the displacement in A and a force applied in A).

jac_A=[1;0;R2];

for ii=1:length(w)
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_ray+K_gen);
    FRF(:,ii)=pinv(jac_A)*H*jac_A;
end
figure(17)
plot(w/2/pi,abs(sum(FRF)),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_A| [m/N]'); title('H_A')

%plot(w/2/pi,(abs(FRF(:,:))),'b')
%% 3C Plot the co-located FRF of point 2 at the centre of the disk 2 with a couple

jac_2=[0;0;1];

for ii=1:length(w)
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_ray+K_gen);
    FRF(:,ii)=pinv(jac_2)*H*jac_2;
end
figure(18)
plot(w/2/pi,abs(sum(FRF)),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_1')


%% 3 D DAMPED - Free motion + Forced motion

lambda_f = [1; 0; R2];
%w = 0:0.01:30;

% X = [H(w)]*Q_0 = [H(w)]*[lambda_f]*F_0
f1=1.5;
f2=3.5;
Fdt = 0;
Xdt = 0;
A1=15;
A2=7;
ampli=1000;
F=[A1,A2];
w_f=[f1,f2]*2*pi;
t=t_damp;
for ii = 1:2 %number of harmonics used to reconstruct the square wave
    X(ii,:) = inv(-(w_f(ii))^2*M_gen + 1i*w_f(ii)*C_ray + K_gen)*lambda_f*F(ii)*exp(1i*pi/2);   
    Fdt = Fdt+F(ii)*sin(w_f(ii)*t); %force
    Xdt = Xdt+real(X(ii,:).*exp(1i*w_f(ii)*[t' t' t']));
end
figure(5)
subplot(2,3,1); plot(w/2/pi,ampli.*abs(FRF11),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_1_1|')
stem(w_f/2/pi,F,'or');
xlim([0 5])
legend('FRF','Force')
subplot(2,3,2); plot(w/2/pi,ampli*abs(FRF21),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_2_1|')
stem(w_f/2/pi,F,'or');xlim([0 5]) 
legend('FRF','Force')
subplot(2,3,3); plot(w/2/pi,ampli*abs(FRF31),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_3_1|')
stem(w_f/2/pi,F,'or'); xlim([0 5])
legend('FRF','Force')
subplot(2,3,4); stem(w_f/2/pi,abs(X(:,1)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [m]'); title('|z|'); xlim([0 5])
subplot(2,3,5); stem(w_f/2/pi,abs(X(:,2)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]'); title('|\theta_1|'); xlim([0 5])
subplot(2,3,6); stem(w_f/2/pi,abs(X(:,3)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]'); title('|\theta_2|'); xlim([0 5])

figure(6)
subplot(2,1,1)
plot(t,Xdt(:,1),t,Xdt(:,2),t,Xdt(:,3));
grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('Forced Motion - z - [m]','Forced Motion - \theta_1 - [rad]','Forced Motion - \theta_2 - [rad]','location','southeast')
title('Forced Time Responce')
subplot(2,1,2)
plot(t,Fdt);
grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('Force - [N]','location','southeast')
title('Time History of the Force')

figure(7)
plot(t,free_motion_damp1.'); 
grid minor;   
xlabel('Time [s]')
ylabel('Amplitude')
legend('Free - z','Free - \theta_1','Free - \theta_2')
title('Free Time Responce')


% Pay attention to the initial conditions!!
figure(8)
plot(t,Xdt+free_motion_damp1.'); 
grid minor;   
xlabel('Time [s]')
ylabel('Amplitude')
legend('Free + Forced responce - z - [m]','Free + Forced responce - \theta_1 - [rad]',...
    'Free + Forced responce - \theta_2 - [rad]')
title('Overall Time Responce')
%%
%excite the system with a square wave
% F(t) = 4/pi*sum( (sin(2k-1)*w*t) / (2k-1) ) = 4/pi*(sin(w*t) + 1/3*sin(3*w*t) + 1/5*sin(5*w*t) + ...)
f0 = 0.75; % Frequency of the applied force
Fdt = 0;
Xdt = 0;

for ii = 1:5 %number of harmonics used to reconstruct the square wave
    k=ii-1;
    n = 2*k+1;
    F(ii) = 8/(pi^2)*((-1)^k)*1/(n^2);
    w_f(ii) = 2*pi*n*f0; %omega force
    X(ii,:) = inv(-(w_f(ii))^2*M_gen + 1i*w_f(ii)*C_gen + K_gen)*lambda_f*F(ii)*exp(1i*pi/2);    
    Fdt = Fdt+F(ii)*sin(w_f(ii)*t); %force
    Xdt = Xdt+real(X(ii,:).*exp(1i*w_f(ii)*[t' t' t']));
end

figure(5)
subplot(2,3,1); plot(w/2/pi,abs(FRF11),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_1_1|')
stem(w_f/2/pi,abs(F),'or'); xlim([0 5])
legend('FRF','Force')
subplot(2,3,2); plot(w/2/pi,abs(FRF21),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_2_1|')
stem(w_f/2/pi,abs(F),'or'); xlim([0 5])
legend('FRF','Force')
subplot(2,3,3); plot(w/2/pi,abs(FRF31),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_3_1|')
stem(w_f/2/pi,abs(F),'or'); xlim([0 5])
legend('FRF','Force')
subplot(2,3,4); stem(w_f/2/pi,abs(X(:,1)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [m]'); title('|z|'); xlim([0 5])
subplot(2,3,5); stem(w_f/2/pi,abs(X(:,2)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]'); title('|\theta_1|'); xlim([0 5])
subplot(2,3,6); stem(w_f/2/pi,abs(X(:,3)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]'); title('|\theta_2|'); xlim([0 5])

figure(6)
subplot(2,1,1)
plot(t,Xdt(:,1),t,Xdt(:,2),t,Xdt(:,3));
grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('Forced Motion - z - [m]','Forced Motion - \theta_1 - [rad]','Forced Motion - \theta_2 - [rad]','location','southeast')
title('Forced Time Responce')
subplot(2,1,2)
plot(t,Fdt);
grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('Force - [N]','location','southeast')
title('Time History of the Force')

figure(7)
plot(t,free_motion_damp.'); 
grid minor;   
xlabel('Time [s]')
ylabel('Amplitude')
legend('Free - z','Free - \theta_1','Free - \theta_2')
title('Free Time Responce')


% Pay attention to the initial conditions!!
figure(8)
plot(t,Xdt+free_motion_damp1.'); 
grid minor;   
xlabel('Time [s]')
ylabel('Amplitude')
legend('Free + Forced responce - z - [m]','Free + Forced responce - \theta_1 - [rad]',...
    'Free + Forced responce - \theta_2 - [rad]')
title('Overall Time Responce')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modal approach
[V,D] = eig(inv(M_gen)*K_gen);
% [V,D] = eig(M_gen\K_gen);

%Natural frequencies
w_nat = sqrt(diag(D)); %[rad/s]
%lambda = alpha +- w_d;
%w0_nat = sqrt(w_d.^2+alpha.^2);
w0_nat = w_nat;

% Consider the undamped case
% [V,D]=eig(M_gen\K_gen);
% w_nat=[sqrt(diag(D)); -sqrt(diag(D))];
% V=[V V]./[V(end,:) V(end,:)];
% [w_nat,index]=sort(w_nat);
% V=V(:,index)./max(abs(V(:,index)));

phi = V./max(abs(V)); %matrix of the mode shapes

% Diagonalization of the matrix
M_mod = phi'*M_gen*phi;
K_mod = phi'*K_gen*phi;

h = [csi_1(2); csi_1(4); csi_1(6)]; %adimensional damping coefficient

const_mod = [1./(2*w0_nat)  w0_nat./2]\h; %pseudo-inverse (least mean square error)

alfa = const_mod(1);
beta = const_mod(2);

C_gen_approx = alfa*M_gen + beta*K_gen;
C_mod = phi'*C_gen_approx*phi;

% Put equal to 0 the extra diagonal terms (< eps = 2.2204e-16)
diag_m = diag(M_mod);
MM = [diag(diag_m)];
diag_k = diag(K_mod);
KK = [diag(diag_k)];
diag_c = diag(C_mod);
CC = [diag(diag_c)];
%%
w = 0:0.01:30;

% Diagonal matrix with eps in extra diagonal positions
for ii = 1:length(w)
    FRF_q = inv(-w(ii)^2*M_mod+1i*w(ii)*C_mod+K_mod);
    FRF_q11(ii) = FRF_q(1,1); FRF_q12(ii) = FRF_q(1,2); FRF_q13(ii) = FRF_q(1,3);
    FRF_q21(ii) = FRF_q(2,1); FRF_q22(ii) = FRF_q(2,2); FRF_q23(ii) = FRF_q(2,3);
    FRF_q31(ii) = FRF_q(3,1); FRF_q32(ii) = FRF_q(3,2); FRF_q33(ii) = FRF_q(3,3);
end
figure(9)
subplot(3,3,1); plot(w/2/pi,abs(FRF_q11),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_1_1');
subplot(3,3,2); plot(w/2/pi,abs(FRF_q12),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_1_2');
subplot(3,3,3); plot(w/2/pi,abs(FRF_q13),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_1_3'); 
subplot(3,3,4); plot(w/2/pi,abs(FRF_q21),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_2_1');
subplot(3,3,5); plot(w/2/pi,abs(FRF_q22),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_2_2');
subplot(3,3,6); plot(w/2/pi,abs(FRF_q23),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_2_3'); 
subplot(3,3,7); plot(w/2/pi,abs(FRF_q31),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_3_1');
subplot(3,3,8); plot(w/2/pi,abs(FRF_q32),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_3_2'); 
subplot(3,3,9); plot(w/2/pi,abs(FRF_q33),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H|[m/N]'); title('H_q_3_3');
sgtitle('Magnitude of the Modal FRF')

figure(10)
subplot(3,3,1); plot(w/2/pi,angle(FRF_q11)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_1_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,2); plot(w/2/pi,angle(FRF_q12)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_1_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,3); plot(w/2/pi,angle(FRF_q13)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_1_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,4); plot(w/2/pi,angle(FRF_q21)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_2_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,5); plot(w/2/pi,angle(FRF_q22)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_2_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,6); plot(w/2/pi,angle(FRF_q23)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_2_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,7); plot(w/2/pi,angle(FRF_q31)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_3_1'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,8); plot(w/2/pi,angle(FRF_q32)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_3_2'); yticks([-180 -90 0 90 180]); ylim([-200,200])
subplot(3,3,9); plot(w/2/pi,angle(FRF_q33)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angleH[deg]'); title('H_q_3_3'); yticks([-180 -90 0 90 180]); ylim([-200,200])
sgtitle('Phase of the Modal FRF')
%% 4.B Reconstruct the co-located FRF of point 𝐴 employing a modal approach and compare with
%the one obtained using physical coordinates in (3.b)
jac_A=[1;0;R2];
jac_qA=phi'*jac_A;
for ii=1:length(w)
    H = inv(-w(ii)^2*M_mod+1i*w(ii)*C_mod+K_mod);
    FRF(:,ii)=pinv(jac_qA)*H*jac_qA;
end
figure(19)
plot(w/2/pi,abs(sum(FRF)),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_A| [m/N]'); title('H_A')

%% 4 C
jac_2=[0;0;1];
jac_q2=phi'*jac_2;
for ii=1:length(w)
    H = inv(-w(ii)^2*M_mod+1i*w(ii)*C_mod+K_mod);
    FRF(:,ii)=pinv(jac_q2)*H*jac_q2;
end
figure(20)
plot(w/2/pi,abs(sum(FRF)),'b'); grid minor; xlabel('Frequency [Hz]'); ylabel('|H_c2| [m/N]'); title('H_c2')

%% Compute the steady state amplitude responce of the system with 2 forces

%% Modal Approach - time histories reconstruction

% Set the initial conditions
pos_t0 = [0.1;pi/12;-pi/12];
vel_t0 = [1;0.5;2];

% Test initial conditions as mode shape
% pos_t0 = V_ord(:,4);
% vel_t0 = 0*V_ord(:,4);

qp_t0 = inv(phi)*vel_t0;
q_t0 = inv(phi)*pos_t0;

% Integration of EOM
A_mod = -inv([M_mod zeros(size(M_mod)); zeros(size(M_mod)) M_mod])*[C_mod K_mod; -M_mod zeros(size(M_mod))];
% A_mod = -inv([MM zeros(size(MM)); zeros(size(MM)) MM])*[CC KK; -MM zeros(size(MM))];
odefun_mod = @(t,y) A_mod*y;
[t_ode_mod,sol_mod] = ode45(odefun_mod,[t(1) t(end)],[qp_t0;q_t0]);

figure(11)
plot(t_ode_mod,sol_mod(:,4),'ob-',t_ode_mod,sol_mod(:,5),'ok-',t_ode_mod,sol_mod(:,6),'or-'); grid minor
xlabel('Time [s]')
legend('q_1','q_2','q_3')
title('Modal Coordinates')

X_rec = (phi*sol_mod(:,4:end)')';

figure(12)
plot(t_ode_mod,X_rec(:,1),'ob-',t_ode_mod,X_rec(:,2),'ok-',t_ode_mod,X_rec(:,3),'or-'); grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('z','\theta_1','\theta_2')
title('Reconstructed time history through modal coordinates')


%% Comparison btw the direct and modal approach
% Solve the direct problem using the approximated [C] matrix
A_damp_mod = -inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*...
    [C_gen_approx K_gen; -M_gen zeros(size(M_gen))]; 
% Integration of EOM
odefun_damp_mod = @(t_damp,y) A_damp_mod*y;
[t_ode_damp_mod,sol_damp_mod] = ode45(odefun_damp_mod,[t_damp(1) t_damp(end)],[vel_t0;pos_t0]);

figure(13)
plot(t_ode_mod,X_rec(:,1),'ob-.',t_ode_mod,X_rec(:,2),'ok-.',t_ode_mod,X_rec(:,3),'or-.');
hold on
plot(t_ode_damp_mod,sol_damp_mod(:,4),'*b',t_ode_damp_mod,sol_damp_mod(:,5),'*k',t_ode_damp_mod,sol_damp_mod(:,6),'*r');
grid minor;
xlabel('Time [s]')
ylabel('Amplitude')
legend('z - modal approach','\theta_1 - modal approach','\theta_2 - modal approach','z - direct approach', '\theta_1 - direct approach', '\theta_2 - direct approach')
title('Damped Time Histoty')

%% Reconstruction of FRF as a linear combination of the mode shapes

figure(14)
ii = 1;
AA = phi(3,ii)*phi(3,ii)*FRF_q11;
ii = 2;
BB = phi(3,ii)*phi(3,ii)*FRF_q22;
ii = 3;
CC = phi(3,ii)*phi(3,ii)*FRF_q33;
TT = AA + BB + CC;
plot(w/2/pi,abs(FRF33),'k','linewidth',1)
hold on
ii = 1;
plot(w/2/pi,abs(AA),'r')
ii = 2;
plot(w/2/pi,abs(BB),'g')
ii = 3;
plot(w/2/pi,abs(CC),'b')
xlabel('Frequency [Hz]')
ylabel('Amplitude')
plot(w/2/pi,abs(TT),'m')
grid minor
legend('H_3_3','H_q_1_1','H_q_2_2','H_q_3_3','Sum')


figure(15)
plot(w/2/pi,angle(FRF33)*180/pi,'k','linewidth',1)
hold on
ii = 1;
plot(w/2/pi,angle(phi(3,ii)*phi(3,ii)*FRF_q11)*180/pi,'r')
ii = 2;
plot(w/2/pi,angle(phi(3,ii)*phi(3,ii)*FRF_q22)*180/pi,'g')
ii = 3;
plot(w/2/pi,angle(phi(3,ii)*phi(3,ii)*FRF_q33)*180/pi,'b')
xlabel('Frequency [Hz]')
ylabel('\angle [deg]')
plot(w/2/pi,angle(TT)*180/pi,'m')
grid minor
legend('H_3_3','H_q_1_1','H_q_2_2','H_q_3_3','Sum')

%% One more example

figure(16)
plot(w/2/pi,abs(FRF13),'k','linewidth',1)
hold on
ii = 1;
plot(w/2/pi,abs(phi(1,ii)*phi(3,ii)*FRF_q11),'r')
ii = 2;
plot(w/2/pi,abs(phi(1,ii)*phi(3,ii)*FRF_q22),'g')
ii = 3;
plot(w/2/pi,abs(phi(1,ii)*phi(3,ii)*FRF_q33),'b')
xlabel('Frequency [Hz]')
ylabel('Amplitude')
grid minor
legend('H_1_3','H_q_1_1','H_q_2_2','H_q_3_3')
 

figure(17)
plot(w/2/pi,angle(FRF13)*180/pi,'k','linewidth',1)
hold on
ii = 1;
plot(w/2/pi,angle(phi(1,ii)*phi(3,ii)*FRF_q11)*180/pi,'r')
ii = 2;
plot(w/2/pi,angle(phi(1,ii)*phi(3,ii)*FRF_q22)*180/pi,'g')
ii = 3;
plot(w/2/pi,angle(phi(1,ii)*phi(3,ii)*FRF_q33)*180/pi,'b')
xlabel('Frequency [Hz]')
ylabel('\angle [deg]')
grid minor
legend('H_1_3','H_q_1_1','H_q_2_2','H_q_3_3')

%% 4 d. Compute the steady state amplitude of response for the three degrees of freedom when
%excited by a horizontal force applied in A. Compare the complete system response with the
%one obtained considering only the first mode of vibration, for the two following cases:
%i. a harmonic force 𝐹(𝑡) = 𝐴1 𝑐𝑜𝑠(2𝜋𝑓1 𝑡);
%ii. a harmonic force 𝐹(𝑡) = 𝐴2 𝑐𝑜𝑠(2𝜋𝑓2 𝑡).

lambda_f = [1; 0; R2];
%phi1=[phi(:,1),zeros(3,2)];

lambda_qf=phi'*lambda_f;
f1=1.5;
f2=3.5;
Fdt = 0;
Xdt = 0;
A1=15;
A2=7;
ampli=A1;
F=A1;
w_f=f1*2*pi;
t=t_damp;
 
    X(1,:) = phi(1,ii)*phi(1,ii)*inv(-(w_f)^2*M_mod + 1i*w_f*C_mod + K_mod)*lambda_qf*F*exp(1i*pi/2);   
    Fdt = Fdt+F*sin(w_f*t); %force
    Xdt = Xdt+real(X(1,:).*exp(1i*w_f*[t' t' t']));

% figure(5)
% subplot(2,3,1); plot(w/2/pi,ampli.*abs(FRF_q11),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_1_1|')
% stem(w_f/2/pi,F,'or');
% xlim([0 5])
% legend('FRF','Force')
% subplot(2,3,2); plot(w/2/pi,ampli*abs(FRF_q22),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_2_1|')
% stem(w_f/2/pi,F,'or');xlim([0 5]) 
% legend('FRF','Force')
% subplot(2,3,3); plot(w/2/pi,ampli*abs(FRF_q33),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude'); title('|H_3_1|')
% stem(w_f/2/pi,F,'or'); xlim([0 5])
% legend('FRF','Force')
% subplot(2,3,4); stem(w_f/2/pi,abs(X(:,1)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [m]'); title('|z|'); xlim([0 5])
% subplot(2,3,5); stem(w_f/2/pi,abs(X(:,2)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]'); title('|\theta_1|'); xlim([0 5])
% subplot(2,3,6); stem(w_f/2/pi,abs(X(:,3)),'b'); grid on; hold on; xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]'); title('|\theta_2|'); xlim([0 5])

figure(6)
subplot(2,1,1)
plot(t,Xdt(:,1),t,Xdt(:,2),t,Xdt(:,3));
grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('Forced Motion - q1','Forced Motion - q2','Forced Motion - q3','location','southeast')
title('Forced Time Responce')
subplot(2,1,2)
plot(t,Fdt);
grid minor
xlabel('Time [s]')
ylabel('Amplitude')
legend('Force - [N]','location','southeast')
title('Time History of the Force')

% figure(7)
% plot(t,free_motion_damp1.'); 
% grid minor;   
% xlabel('Time [s]')
% ylabel('Amplitude')
% legend('Free - z','Free - \theta_1','Free - \theta_2')
% title('Free Time Responce')
% 
% 
% % Pay attention to the initial conditions!!
% figure(8)
% plot(t,Xdt+free_motion_damp1.'); 
% grid minor;   
% xlabel('Time [s]')
% ylabel('Amplitude')
% legend('Free + Forced responce - z - [m]','Free + Forced responce - \theta_1 - [rad]',...
%     'Free + Forced responce - \theta_2 - [rad]')
% title('Overall Time Responce')
 


