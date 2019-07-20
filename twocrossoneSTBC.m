%%%%%%%%%%%%% two cross one STBC %%%%%%%%%%%%%%%%
N_frame=130; N_packet=4000;NT=2; NR=1; b=2;
SNRdBs=[0:2:30]; sq_NT=sqrt(NT); sq2=sqrt(2);
for i_SNR=1:length(SNRdBs)
SNRdB=SNRdBs(i_SNR); sigma=sqrt(0.5/(10^(SNRdB/10)));
for i_packet=1:N_packet
msg_symbol=randint(N_frame*b,NT);

%%%%%%%%%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx_bits=msg_symbol.'; tmp=[]; tmp1=[];

for i=1:NT
[tmp1,sym_tab,P]=modulator(tx_bits(i,:),b); tmp=[tmp; tmp1];
end
X=tmp.'; X1=X; X2=[-conj(X(:,2)) conj(X(:,1))];

%%%%%%%%%%%%%%%%%%%%% noise channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:NT
    hr(n,:,:)=(randn(N_frame,NT)+j*randn(N_frame,NT))/sq2;       %noise
end

h=reshape(hr(n,:,:),N_frame,NT); Habs(:,n)=sum(abs(h).^2,2);

%%%%%%%%%%%%%%%%%%%% receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = sum(h.*X1,2)/sq_NT+sigma*(randn(N_frame,1)+j*randn(N_frame,1));   % received signal + noise matrix for time T1
r2 = sum(h.*X2,2)/sq_NT+sigma*(randn(N_frame,1)+j*randn(N_frame,1));   % received signal + noise matrix for time T2

Z1 = r1.*conj(h(:,1)) + conj(r2).*h(:,2);
Z2 = r1.*conj(h(:,2)) - conj(r2).*h(:,1);

%%%%%%%%%%%%% estimation and distance minimization %%%%%%%%%%%%%%%%%%%%%%%
for m=1:P
xi = (-1+sum(Habs,2))*abs(sym_tab(m))^2;
e1(:,m) = abs(sum(Z1,2)-sym_tab(m)).^2 + xi;    % estimated signal for T1
e2(:,m) = abs(sum(Z2,2)-sym_tab(m)).^2 + xi;    % estimated signal for T2
end
[y1,i1]=min(e1,[],2); S1d=sym_tab(i1).'; clear e1    % row distance minimization for time T1
[y2,i2]=min(e2,[],2); S2d=sym_tab(i2).'; clear e2    % row distance minimization for time T2
Xd = [S1d S2d];                                      % decoded signal
tmp1=X>0 ; tmp2=Xd>0;
noeb_p(i_packet) = sum(sum(tmp1~=tmp2));% for coded
end % end of FOR loop for i_packet
BER(i_SNR) = sum(noeb_p)/(N_packet*N_frame*b);
end % end of FOR loop for i_SNR
semilogy(SNRdBs,BER), axis([SNRdBs([1 end]) 1e-6 1e0]);
grid on; xlabel('SNR[dB]'); ylabel('BER');