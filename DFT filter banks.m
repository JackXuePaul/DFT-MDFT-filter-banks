%%
%clear all;
close all;
fs = 10e5;                   %sampling frequency
fc = fs/2;
M = 8;                     % the number of channels & the factor of downsampling(upsampling)
%Kds=6;                     
T = 1;                        %Sampling time
t = T/fs:(T/fs):T;

xx = 2*cos(2*pi*(fc/M*3)*t)+cos(2*pi*fc/M*6.2*t);  %input signal
lx=length(xx);

XX=zeros(M,lx/M);

% all above could be changed if needed;
XX=flipud(reshape(xx,M,lx/M));   %the signal matrix
%%
%design of prototype filter
[n,fo,ao,w] = firpmord([fs/(2*M) fs/(2*M)*1.05],[1 0],[0.001 0.0001],fs);
n=n-mod(n,M)+M-1;
h = firpm(n,fo,ao,w);
%fvtool(h);
lh = length(h);

%%

E=reshape(h,M,lh/M);     %analysis filters

T=zeros(M,lx/M);
[rt,ct]=size(T);
for i = 1:1:M
    T(i,:)=filter(E(i,:),1,XX(i,:));
end

%V=(1./dftmtx(M))*T;
V=ifft(T);
[rv,cv]=size(V);

%%
%synthesis part
Q=fft(V);
R=flipud(reshape(h,M,lh/M));      %synthesis filters
Y=zeros(M,lx/M);
X=zeros(M,lx);    
[ry,cy]=size(Y);
for i = 1:1:M
    Y(i,:)=filter(R(i,:),1,Q(i,:));
    X(i,:)=upsample(Y(i,:),M);
    X(i,:)=circshift(X(i,:),1-i);
end
r_x = sum(X);    %reconstruction signal, comparing with the input xx;

figure();
plot(abs(fft(xx)));
figure();
plot(abs(fft(r_x)));
  