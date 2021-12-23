clear all;
close all;

%%%%%%% ��������
fd=8000; % �ź�Դ�����ٶ�

%%%%%%% �źŲ���

m1=mseq([0 0 1 0 1]); % ������һ��31λ�����ź�
m2=mseq([1 0 1 1 1]); % �����ڶ���31λ�����ź�

tm=linspace(0,1/fd*30,31);
subplot(211)
stairs(tm,m1)
axis([0 0.0038 -0.1 1.1])
title('�ź�Դ1')
subplot(212)
stairs(tm,m2)
axis([0 0.0038 -0.1 1.1])
title('�ź�Դ2')
pause;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Hamming�ŵ�����
m1=[m1 0];
m2=[m2 0];
n=7;%�볤
k=4;%��Ϣλ��
P=[1 1 1; 1 1 0;1 0 1;0 1 1];
%P�����˼ල��ϵ
G=[eye(k) P];%���ɾ���
H=[P' eye(n-k)];%У�����
%����
th=linspace(0,1/fd*63,64);
ham1=zeros(1,8*length(m1)/4); % 64
ham2=zeros(1,8*length(m2)/4); % length=64
%�����λ0�����8λ
for i=1:length(m1)/4
    ham1(1,(8*(i-1)+1):8*i)=[mod(m1(1,4*(i-1)+1:4*i)*G,2) 0];
    ham2(1,(8*(i-1)+1):8*i)=[mod(m2(1,4*(i-1)+1:4*i)*G,2) 0];    
end
subplot(311)
stairs(th,ham1)
axis([0 0.0079 -0.1 1.1])
title('�ź�Դ1�ı��������')
subplot(312)
stairs(th,ham2)
axis([0 0.0079 -0.1 1.1])
title('�ź�Դ2�ı��������')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ʱ�ָ��� 
m=[];
for i=1:length(m1)/4
    m=[m ham1(1,(8*(i-1)+1):8*i) ham2(1,(8*(i-1)+1):8*i)];
end
subplot(313)
tt=linspace(0,(length(m)-1)*1/fd,length(m));
stairs(tt,m)
axis([0 0.0160 -0.1 1.1])
title('ʱ�ָ���')
pause;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%FSK����
fc1=32000;
fc2=64000;
Nds2=10;
len=length(m);
tc=linspace(0,1/fc2/10*(len*10*Nds2-1),len*10*Nds2);   %�����ز��Ķ������������
ycc1=cos(2*pi*fc1*tc);
ycc2=cos(2*pi*fc2*tc);
yout=[];
%����
for i=1:len
    if(m(i)>0.5)
        yout=[yout ycc2(1,((i-1)*Nds2*10+1):(i*10*Nds2))];
    else
        yout=[yout ycc1(1,((i-1)*Nds2*10+1):(i*10*Nds2))];
    end
end
subplot(211)
plot(tc,yout)
title('FSK��������')
axis([3.5e-3 5.5e-3 -1.1 1.1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
yout = awgn(yout,3) % ����5dB����
subplot(212)
plot(tc,yout)
title('��������')
axis([3.5e-3 5.5e-3 -1.1 1.1])
pause;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FSK���
Fs = 640000;  % ����Ƶ��
N    = 40;       % Order
Fc1  = 31000;     % First Cutoff Frequency
Fc2  = 33000;     % Second Cutoff Frequency
Fc11=63000;
Fc22=65000;
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);
% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd1 = dfilt.dffir(b);
b  = fir1(N, [Fc11 Fc22]/(Fs/2), 'bandpass', win, flag);
Hd2 = dfilt.dffir(b);

yc1=filter(Hd1,yout); % �˲�
yc2=filter(Hd2,yout);

Nds=10;
len=length(yc1);
num=Nds*10;
out2=[];
for i=1:len/num
    if(sum(abs(yc2(1,((i-1)*num+1):(i*num))))>30)
        out2=[out2 1];
    else
        out2=[out2 0];
    end
end
subplot(211)
plot(tc,yc2)
title('����ź�')
subplot(212)
stairs(tt,out2)
axis([0 0.0160 -0.1 1.1])
title('����ź�')
pause;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%�⸴��
reham1=[];
reham2=[];
for i=1:length(out2)/16
    reham1=[reham1 out2(1,((i-1)*16+1):((i-1)*16+1+7))];
    reham2=[reham2 out2(1,((i-1)*16+9):((i-1)*16+9+7))];
end
subplot(411)
stairs(ham1)
axis([0 64 -0.1 1.1])
title('������1')
subplot(412)
stairs(reham1)
axis([0 64 -0.1 1.1])
title('�⸴�ú�����1')
subplot(413)
stairs(ham2)
axis([0 64 -0.1 1.1])
title('������2')
subplot(414)
stairs(reham2)
axis([0 64 -0.1 1.1])
title('�⸴�ú�����2')
pause;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%����
for i=1:length(reham1)/8
    s=mod(reham1(1,(8*(i-1)+1:(8*i-1)))*H',2);
    %Ѱ�Ҵ������ֵ�λ��
    pos=bin2dec(num2str(s));
    if(pos~=0)
        if(pos==4)
            pos=7-2;
        elif(pos==3)
            pos=7-3;
        else
            pos=7-pos+1;
        end
    end
    if(pos)
        reham1(1,8*(i-1)+pos)=~reham1(1,8*(i-1)+pos);
    end
    s=mod(reham2(1,(8*(i-1)+1:(8*i-1)))*H',2);
    %Ѱ�Ҵ������ֵ�λ��
    pos=bin2dec(num2str(s));
    if(pos~=0)
        if(pos==4)
            pos=7-2;
        elseif(pos==3)
            pos=7-3;
        else
            pos=7-pos+1;
        end
    end
    if(pos)
        reham2(1,8*(i-1)+pos)=~reham2(1,8*(i-1)+pos);
    end
end
subplot(411)
stairs(ham1)
axis([0 64 -0.1 1.1])
title('������1')
subplot(412)
stairs(reham1)
axis([0 64 -0.1 1.1])
title('���������1')
subplot(413)
stairs(ham2)
axis([0 64 -0.1 1.1])
title('������2')
subplot(414)
stairs(reham2)
axis([0 64 -0.1 1.1])
title('���������2')
pause;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%����
rem1=[];
rem2=[];
for i=1:length(reham1)/8
    rem1=[rem1 reham1(1,((i-1)*8+1):((i-1)*8+1+3))];
    rem2=[rem2 reham2(1,((i-1)*8+1):((i-1)*8+1+3))];
end
rem1=rem1(1,1:length(rem1)-1);
rem2 =rem2(1,1:length(rem2)-1);    
subplot(411)
stairs(tm,m1(1,1:31))
axis([0 0.0038 -0.1 1.1])
title('ԭ��1')
subplot(412)
stairs(tm,rem1)
axis([0 0.0038 -0.1 1.1])
title('����������1')
subplot(413)
stairs(tm,m2(1,1:31))
axis([0 0.0038 -0.1 1.1])
title('ԭ��2')
subplot(414)
stairs(tm,rem2)
axis([0 0.0038 -0.1 1.1])
title('����������2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����������
error_bit1=sum(rem1~=m1(1,1:31)); 




