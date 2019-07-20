%%%%%two cross two %%%%%%%%%%
clear all
n_frame=1000;
n_packet=100;
b=2;
m=2^b;
mod_obj=modem.qammod('m',m,'symbolorder','gray','inputtype','bit');
demod_obj=modem.qamdemod(mod_obj);
%MIMO parameters
t_tx=4;
code_length=64;
nt=2;
nr=2;
n_pbits=nt*b*n_frame;n_tbits=n_pbits*n_packet;
code_book=codebook_generator;
fprintf('===============\n');
fprintf('precoding transmission');
fprintf('n %d x %d MIMO\n %d QAM',nt,nr,m);
fprintf('\n simulation bits : %d',n_tbits);
fprintf('\n=============\n');
SNRdbs = [0:2:10];
sq2=sqrt(2);
for i_SNR=1:length(SNRdbs)
    SNRdb=SNRdbs(i_SNR)
    noise_var=nt*0.5*10^(-SNRdb/10);sigma=sqrt(noise_var);
    rand('seed',1);randn('seed',1);n_ebits=0;
    for i_packet=1:n_packet
        msg_bit=randint(n_pbits,1);
        
        %%%%%%%transmitter%%%%%%%%
        
        s=modulate(mod_obj,msg_bit);
        scale=modnorm(s,'avpow',1);
        S=reshape(scale*s,nt,1,n_frame);
        tx_symbol=[S(1,1,:) -conj(S(2,1,:));S(2,1,:) conj(S(1,1,:))];
        
        
        %%%%%%%%%%%%%%%%%%%channel noise%%%%%%%%%%%%%%%%%%%%
        
       hr=(randn(nr,t_tx)+j*randn(nr,t_tx))/sq2;
         
            %hr1=(randn(nr,t_tx)+j*randn(nr,t_tx))/sq2;
        
        %h=[hr;hr1];
        for i=1:code_length
            cal1(i)=norm(hr*code_book(:,:,i),'fro');
        end
        %for k=1:code_length
           % cal2(k)=norm(hr1*code_book(:,:,k),'fro');
        %end
        [val,index]=max(cal1);
        %[val1,index1]=max(cal2);
        he1=hr*code_book(:,:,index);
        %he2=hr1*code_book(:,:,index1);
        %he=[he1;he2];
        norm_h1=norm(he1)^2;
         %norm_h2=norm(he2)^2;
        for i=1:n_frame
         r1(:,:,i)=he1*tx_symbol(:,:,i)+sigma*(randn(nr,2)+j*randn(nr,2));
       % r2(:,:,i)=he2*tx_symbol(:,:,i)+sigma*(randn(nr,2)+j*randn(nr,2));
        end
        
        %%%%%%%%%%%%%%%%%receiver%%%%%%%%%%%%%%%%%%
      
        for i=1:n_frame
            y(1,i,:)=((he1(1,1)'*r1(1,1,i)'+he1(1,2)*r1(2,1,i)')+(he1(2,1)'*r1(1,1,i)'+he1(2,2)*r1(2,2,i)'))/norm_h1;   %yT1
            y(2,i,:)=((he1(1,2)'*r1(2,2,i)'-he1(2,1)*r1(1,2,i)')+(he1(2,2)'*r1(2,2,i)'-he1(1,2)*r1(2,2,i)'))/norm_h1;   %yT2
            
        end
        s_hat=reshape(y/scale,nt*n_frame,1);
        msg_hat=demodulate(demod_obj,s_hat);
        n_ebits=n_ebits+sum(msg_hat~=msg_bit);
    end
    BER(i_SNR)=n_ebits/n_tbits;
end
semilogy(SNRdbs,BER,'-k^','linewidth',2); hold on;grid on;
xlabel('SNR[db]'),ylabel('BER');legend('precoded');