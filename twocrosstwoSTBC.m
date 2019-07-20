% two cross two%
clear
clc
N_frame=130; 
N_packet=4000;
NT=2; 
NR=2; 
b=2;
SNRdBs=[0:2:30]; 
sq_NT=sqrt(NT); 
sq2=sqrt(2);
for i_SNR=1:length(SNRdBs)
    SNRdB=SNRdBs(i_SNR); 
    sigma=sqrt(0.5/(10^(SNRdB/10)));
    for i_packet=1:N_packet
        msg_symbol=randint(N_frame*b,NT);
        
        %%%%%%%%%%%%% transmitter side %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tx_bits=msg_symbol.'; tmp=[]; tmp1=[];
        for i=1:NT
            [tmp1,const,M]=modulator(tx_bits(i,:),b); tmp=[tmp; tmp1];     %
        end
        X=tmp.'; X1=X; X2=[-conj(X(:,2)) conj(X(:,1))];
        
        %%%%%% noise channel %%%%%%%%%%%%%%%%
        
        for n=1:NT
            n1(n,:,:)=(randn(N_frame,NT)+j*randn(N_frame,NT))/sq2;
        end
        for n=1:NT
            n2(n,:,:)=(randn(N_frame,NT)+j*randn(N_frame,NT))/sq2;
        end
        h1=reshape(n1(n,:,:),N_frame,NT); Habs(:,n)=sum(abs(h1).^2,2);
        h2=reshape(n2(n,:,:),N_frame,NT); Habs1(:,n)=sum(abs(h2).^2,2);
        
        %%%%%%%%%%%%%%%%%%%%% receiver side %%%%%%%%%%%%%%%%%%%%%
        
        r11 = sum(h1.*X1,2)/sq_NT+sigma*(randn(N_frame,1)+j*randn(N_frame,1)); % receiver matrix at rx1 at time T1
        r21 = sum(h1.*X2,2)/sq_NT+sigma*(randn(N_frame,1)+j*randn(N_frame,1));  % receiver matrix at rx2 at time T1
        r12 = sum(h2.*X1,2)/sq_NT+sigma*(randn(N_frame,1)+j*randn(N_frame,1));  % receiver matrix at rx1 at time T2
        r22 = sum(h2.*X2,2)/sq_NT+sigma*(randn(N_frame,1)+j*randn(N_frame,1)); % receiver matrix at rx2 at time T2
       
        Z1 = r11.*conj(h1(:,1)) + conj(r21).*h1(:,2) + r12.*conj(h2(:,1)) + conj(r22).*h2(:,2);
        Z2 = r11.*conj(h1(:,2)) - conj(r21).*h1(:,1) + r12.*conj(h2(:,2)) - conj(r22).*h2(:,1);
                       %%%%%%%%%%%%% estimation and distance minimization %%%%%%%%%%%%%%%%%%%%%%%
          for m=1:M
           xi = (-1+sum(Habs,2))*abs(const(m))^2;
           e1(:,m) = abs(sum(Z1,2)-const(m)).^2 + xi;    % estimated signal for T1
           e2(:,m) = abs(sum(Z2,2)-const(m)).^2 + xi;    % estimated signal for T2
          end
              [y1,i1]=min(e1,[],2); S1d=const(i1).'; clear e1    % row distance minimization for time T1
              [y2,i2]=min(e2,[],2); S2d=const(i2).'; clear e2    % row distance minimization for time T2
              Xd = [S1d S2d];                                      % decoded signal
               tmp1=X>0 ; tmp2=Xd>0;
              noeb_p(i_packet) = sum(sum(tmp1~=tmp2));% for coded
          end % end of FOR loop for i_packet
   BER(i_SNR) = sum(noeb_p)/(N_packet*N_frame*b);
      end % end of FOR loop for i_SNR
   semilogy(SNRdBs,BER), axis([SNRdBs([1 end]) 1e-6 1e0]);
   grid on; xlabel('SNR[dB]'); ylabel('BER');