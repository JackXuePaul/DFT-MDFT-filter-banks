%%
clear all;
close all;
load initial_coefficients.mat;
load input.mat;
fs = 10e5;                   %sampling frequency
fc = fs/2;
M = 64;                     %the factor of downsampling(upsampling)
I=2;                        %oversample ratio
K=M*I;                      %number of channels
T = 1;                     %Sampling time

% B=fs/8/K;                   %band width
% k=B/T;%调频斜率
% n=round(T*fs);%采样点个数
t = 1/fs:(1/fs):(T+(2*M-mod(fs*T,2*M))*1/fs);
x=[zeros(1,6) xx zeros(1,6)];
% xx=exp(1j*fs/K*12*tt).*exp(1j*pi*k*tt);%LFM信号
%xx = 2*cos(2*pi*fc/M*3.2*tt)+cos(2*pi*fc/M*3.25*tt);   %input signal
lx=length(x);
XX=zeros(K,lx);
X=zeros(K,lx/K);
X=reshape(x,K,lx/K);

%%
%滤波器设计
disp('design prototype filter ...');
m = 8;                      %m detemines the length of filter
E = 1e-8;                   
a = 1e5;                    %the energy ratio of passband and stopband               
t = 0.6;
N = 2*m*M;                  %length of filter
F=10; L = 2*F*M;
freq_vector = linspace (0, pi, L);
Ut = zeros(F,N/2);
Us = zeros(F,N/2);
for k = 1:F
    for n = 1:N/2
        Ut(k,n) = 2*cos(((n-1)-(N-1)/2)*freq_vector(k));
        Us(k,n) = 2*cos(((n-1)-(N-1)/2)*(freq_vector(k)-pi/M));
    end
end
mU = zeros(L-2*F,N/2);
for k = 2*F+1:L
      for n = 1:N/2
        mU(k,n) = 2*cos(((n-1)-(N-1)/2)*freq_vector(k));
      end
end
Q = mU'*mU;

b_vector = fcoe;
h_n_m = b_vector(1,1:N/2)';
Error = 0.01;
c = 0;          %number of iterations
while (Error>E)
    U = diag(Ut*h_n_m)*Ut+diag(Us*h_n_m)*Us;
    g_n_m = (U'*U+a*Q)\U'*ones (F, 1);
    h_n_m = (1-t)*h_n_m+t*g_n_m;
    Error = norm (h_n_m-g_n_m)/norm (h_n_m);
    c = c+1;
end
h = [h_n_m;fliplr(eye(N/2))*h_n_m];
h=h';
h=K*h;
lh=length(h);

%%
%multiphase filters
% s=floor(lh/M/2);
% d=2*M-1;
% D=2*s*M+d;     %delay time
% for k=1:2*M
%     for n=1:lh
%         hk(k,n)=sqrt(2)*h(1,n)*exp(sqrt(-1)*pi*(k-1)/M*(n-D/2));     %multiphase filter banks
%     end
% end
%  tic;
P=zeros(K,lh);
p=zeros(K,lh/M);       %analysis polyphase fillters
pk=zeros(2*K,lh/K);    %polyphase of polyphase
px=zeros(2*K,lx/K);
Y=zeros(K,lx/M);

for i = 1:K
    P(i,i:end)=h(1,1:lh-i+1);
    p(i,:)=P(i,1:M:end);
end

for i =1:K
    pk(2*i-1:2*i,:)=reshape(p(i,:),I,lh/K);
    px(2*i-1,:)=filter(pk(2*i-1,:),1,X(i,:));
    px(2*i,:)=filter(pk(2*i,:),1,X(i,:));
end
for i =1:K
    Y(i,:)=reshape(px(2*i-1:2*i,:),1,lx/M);
end
y=ifft(Y);
%%
XK=zeros(4*M,lx/2/M);
YK=zeros(K,lx/M);
yk=zeros(K,lx/M);
YY=zeros(2*K,lx/K);
RX=zeros(K,lx/K);

for i = 1:K    
     XK(2*i-1:2*i,:)=reshape(y(i,:),2,lx/K);
end
for i = 1:M
    XK(4*i-3,:)=j*imag(XK(4*i-3,:));
    XK(4*i,:)=j*imag(XK(4*i,:));
    XK(4*i-2,:)=real(XK(4*i-2,:));
    XK(4*i-1,:)=real(XK(4*i-1,:));
end

for i =1:K
    YK(i,:)=reshape(XK(2*i-1:2*i,:),1,lx/M);
end

yk=fft(YK);
%%
%综合
Q=zeros(K,lh);
q=zeros(K,lh/M);
qk=zeros(2*K,lh/K);       %polyphase of polyphase
xq=zeros(2*K,lx/K);

for i = 1:K
    Q(i,1:lh-i+1)=h(1,i:lh);
    q(i,:)=Q(i,1:M:end);
end
for i = 1:K
    qk(2*i-1:2*i,:)=reshape(q(i,:),I,lh/K);
    xq(2*i-1,:)=downsample(yk(i,:),I);
    xq(2*i,:)=downsample(circshift(yk(i,:),1),I);
end

for i =1:2*K
    YY(i,:)=filter(qk(i,:),1,xq(i,:));
end

for i = 1:K
    RX(i,:)=YY(2*i-1,:)+YY(2*i,:);
end

r_x=reshape(RX,1,lx); %reconstruct signal;



[td,lags]=xcorr(x,r_x);
figure();
stem(lags,abs(td));

rx=x(1,1:end-1024)-r_x(1,1025:end);
figure();
plot(rx);
figure();
plot(abs(fft(xx)));
figure();
plot(abs(fft(r_x)));

% N=fs/M*T;
% n=nst channel
% L=the length of nst channel
% if L>=N/2
%    op=(N-L)*fs/M/N;
% end
% if L<N/2
%     op=(N/4-L)*fs/M/N+(n-1.5)*fs/K;
% end
