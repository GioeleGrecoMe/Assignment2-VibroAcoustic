function h_vec=HalfPowerMethod(f,H) 
    %f=frequency domain and H=transfer function
    %1 normalize H
    H_norm=abs(H./(H(1)));
    %2 find peaks and valley
    [pks,idx]=findpeaks(abs(H_norm));
    [valley,vidx]=findpeaks(abs(1./H_norm));
    vidx=[1;vidx;length(H)];
    %3 find w0_i
    w0=f(idx)*2*pi;
    n=length(w0);
    %4 find w1 and w2 as half power H
    for i=1:n
            start=(vidx(i));
            H_filt=abs(H_norm(vidx(i):vidx(i+1)));
     
        H_pow=H_filt.*(H_filt>=(H_norm(idx(i))/sqrt(2)));
        HPowInd=find(H_pow>0);
        w1=f(HPowInd(1)+start)*2*pi;
        w2=f(HPowInd(end)+start)*2*pi;
        h_vec(i)=(w2^2-w1^2)/(4*w0(i)^2);
    end
    h_vec=h_vec';
end