clear all

%% load data
load('Data')
%n_ext=10;
n_z=100000;
t=[Data(:,1);zeros(n_z,1)];
%t_ext=interp(t,n_ext);
%t=t_ext;
F=[Data(:,2);zeros(n_z,1)];
%F_ext=interp(F,n_ext);
%F=F_ext;
xf=Data(:,3:end);
x(:,1)=[xf(:,1);zeros(n_z,1)];
x(:,2)=[xf(:,2);zeros(n_z,1)];
x(:,3)=[xf(:,3);zeros(n_z,1)];
x(:,4)=[xf(:,4);zeros(n_z,1)];
%x=[interp(x(:,1),n_ext),interp(x(:,2),n_ext),interp(x(:,3),n_ext),interp(x(:,4),n_ext)];
np=size(x,1);
nj=size(x,2);

dt=t(2)-t(1);
fsamp=1/dt;
disp(' ')
disp(['Number of sensors: ' num2str(nj)])

% figure(1)
% sf1(1)=subplot(3,1,1);
% plot(t,F)
% grid on
% ylabel('F [N]')
% title('Force')
% sf1(2)=subplot(3,1,2);
% plot(t,x(:,1),'b',t,x(:,2),'r',t,x(:,3),'b',t,x(:,4),'r')
% grid on
% ylabel('x_1,x_2,x_3,x_4 [m]')
% title('Displacement')
% sf1(3)=subplot(3,1,3);
% plot(t,x(:,3),'b',t,x(:,4),'r')
% grid on
% ylabel('x_3,x_4 [m]')
% title('Displacement')
% linkaxes(sf1,'x')
% sgtitle('Measurements')


%% 1) PLOT EXPERIMENTAL FRF diagrams Definition of experimental Transfer Function

% force and measurement spectra
[xfft,frq]=ffg(x,np,dt);
[Ffft,frq]=ffg(F,np,dt);

% transfer function
Hjkexp=xfft./Ffft;

figure(2)
sf2(1)=subplot(2,1,1);
for jj=1:nj
    plot(frq,abs(Hjkexp(:,jj)))
    hold on
    legenda{jj}=['H_' num2str(jj) '_k'];
end
grid on
title('Magnitude')
sgtitle('Experimental Transfer Functions H_j_k')
ylabel(['|H_j_k| [m/N]'])
legend(legenda)
sf2(2)=subplot(2,1,2);
for jj=1:nj
    plot(frq,angle(Hjkexp(:,jj))*180/pi)
    hold on
end
grid on
ylabel('\angleH_j_k [deg]')
xlabel('Freq [Hz]')
yticks([-180 -90 0 90 180])

title('phase')
linkaxes(sf2,'x')
xlim([0 5])
%% 2) Estimate the natural frequencies, damping ratios and mode shapes of the resonating modes in the
%range 0 − 5 𝐻𝑧 employing simplified methods (e.g. half power point method). Comment the
%obtained results.
FRF=Hjkexp;
ind=find(frq<5);
for i=1:4
    [pks(:,i),idx(:,i)]=findpeaks(abs(FRF(1:(ind(end)-15),i)));
    pks_im(:,i)=FRF(idx(:,i),i);
    w_nat(:,i)=frq(idx(:,1))*2*pi;%idx(:,i)/duration*2*pi;
    %with half power method
    h_vec(:,i)=HalfPowerMethod(frq(1:ind(end)),FRF(1:(ind(end)-15),i));
end


%with mqii=1 normalized we can calculate cqii passing through the h_vec
%hii=cii/(2mii*w0i)=>cii=hii*2*w0i
for i=1:4
   cq(:,i)=h_vec(:,i)*2.*w_nat(:,i);
   X(:,i)=cq(:,i).*w_nat(:,i).*pks_im(:,i);
end
X_norm=X./(X(1,:));

%% identification parameters single range
nrange=4;
[pks,idx_min]=(findpeaks(1./sum(abs(FRF(1:(ind(end)),:))')'));
range_min=[frq(80);frq(idx_min)];
range_max=[frq(idx_min-2);frq(ind(end)-300)];
%%
for irangle=1:nrange
    fini=range_min(irangle);
    ffin=range_max(irangle);
    iini=min(find(round(frq*1000)/1000>=fini));
    ifin=max(find(round(frq*1000)/1000<=ffin));
    disp(' ')
    disp('Frequency range')
    disp(['- fmin [Hz]: ' num2str(frq(iini))])
    disp(['- fmax [Hz]: ' num2str(frq(ifin))])
    npid=ifin-iini+1;
    disp(['Number of points for the identification: ' num2str(npid)])
    
    % Function reduced to the range of interest
    rfHjki=frq(iini:ifin);
    Hjkiexp=Hjkexp(iini:ifin,:);
    
    % First guess parameters: SIMPLIFIED METHODS
    for nx=1:nj

        disp(' ')
        jj=nx;%input('Which diagram for identification? ');
        %modal parameters guessing
        [vmax,iwmax]=max(abs(Hjkiexp(:,jj)));
        f0i=rfHjki(iwmax);
        w0i0=2*pi*f0i;
        derFIjki=(angle(Hjkiexp(iwmax+1,jj))-angle(Hjkiexp(iwmax-1,jj)))/(2*pi*(rfHjki(iwmax-1)));
        csii0=-1/(w0i0*derFIjki);
        r0i=2*w0i0*csii0;
        % Constant parameters guessing
        Aj0=-imag(Hjkiexp(iwmax,jj))*w0i0*r0i;
        % all the other constants are 0
        
        disp(' ')
        disp(['Displacement x_' num2str(jj)])
        disp(['Init guess f0 [Hz]: ' num2str(f0i)])
        disp(['Init guess csi [-]: ' num2str(csii0)])
        disp(['Init guess A0 [Hz]: ' num2str(Aj0)])
        
        %filling of vector xpar
        
        xpar0=[csii0;w0i0;Aj0;zeros(5,1)];
        
        %identification: single channel
        options=optimset('fminsearch');
        options=optimset(options,'TolFun',1e-8,'TolX',1e-8);
        xpar=fminsearch(@(xpar) errHjki_cw(xpar,rfHjki,Hjkiexp(:,jj)),xpar0,options);
        %plot resuts of identification
        vpar=[1;2*xpar(1)*xpar(2);xpar(2)^2;xpar(3:8)];
        %    [m; c=2 m w0 csi; k=w0^2 m; A B C D E F]
        
        Hjkiid=funHjki(vpar,rfHjki);
        
        %
        
        figure(3)
        sf3(1)=subplot(2,1,1);
        plot(rfHjki,abs(Hjkiexp(:,jj)),'*b-',rfHjki,abs(Hjkiid),'r--','linewidth',2);
        grid on
        xlabel('Frequency [Hz]')
        ylabel(['|H' num2str(jj) 'k| [m/N]'])
        title('Magnitude')
        legend('Esperimental','Identified')
        sf3(2)=subplot(2,1,2);
        plot(rfHjki,angle(Hjkiexp(:,jj)),'*b-',rfHjki,angle(Hjkiid),'r--','linewidth',2);
        grid on
        ylabel(['\abgleH' num2str(jj) 'k [rad]'])
        xlabel('Frequency [Hz]')
        title('Phase')
        legend('Experimental','Identified')
        linkaxes(sf3,'x')
        sgtitle(['Transfer functions H' num2str(jj)])
    end
end
    
%% 4) Compare and comment the identified modal parameters with ones obtained with the abovementioned methods