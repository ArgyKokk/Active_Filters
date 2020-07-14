% Argyrios Kokkinis
% 8459
% Low pass filter
% Butterworth

a1 = 8;
a2 = 4;
a3 = 5;
a4 = 9;

%specs
fp = 0.6 *(3+2);
fs = 2 * fp;
fp = fp*1000;
fs = fs*1000; 
a_min = 17.5 + (9-5)*0.5;
a_max = 0.6 + (5-3)/10;

% pass band from 0 to 2ðfp (0:6000ð) , a < 0.8 dB
% stop band from 2ðfs to inf (12000ð:inf), a > 19.5 dB
% transition band from 2ðfp to 2ðfs (6000ð:12000ð)

ws = 2*pi*fs;
wp = 2*pi*fp;
n = log((10^(a_min/10) - 1) / (10^(a_max/10) - 1)) / (2*log(ws / wp));
n = ceil(n)

w0 = wp / ((10^(a_max/10) - 1)^(1/(2*n)));

%poles : s1,2 = -0.8090 +- 0.5877j , s2,3 = -0.3090 +- 0.9510j , s5 = -1
%Q : 0.5 , 0.62 , 1.62
%angles: 0, +-36 ,+- 72
Q1 = 0.5;
Q2 = 0.62;
Q3 = 1.62;

%strategy 1
% gains
k1 = 1;
k2 = 3 - 1/Q2;
k3 = 3 - 1/Q3;

%Transfer functions
s = tf('s');
H1 = w0/(s+w0);
H2 = k2 * (w0^2)/(s^2+w0/Q2*s + w0^2);
H3 = k3 * (w0^2)/(s^2+w0/Q3*s + w0^2);

%gain specification K=1 , 0dB
K = 10^(0);
k = K/(k1*k2*k3);

%Transfer function
H = k*H1*H2*H3;

%create triangular periodic signal for 10 periods
T = 1/2000;
Total_time = 10*T;
f_s = 200000;
t=0:1/f_s:Total_time-1/f_s;
sig = sawtooth(2*pi*2000*t,1/2);

%plot transient analysis
out=lsim(H,sig,t);
plot(t,sig,'red')
hold on;
plot(t,out,'green')
title('Transient Analysis');
ylabel('Magnitude');
hold off;

%input signal fourier
N = 1000;
n=2^nextpow2(N);
xfourier= fft(sig,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
plot(f,p1,'red')
title('input spectrum');
xlabel('Frequency [Hz]');

%output signal fourier
N = 1000;
n=2^nextpow2(N);
xfourier= fft(out,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
plot(f,p1,'red')
title('output spectrum');
xlabel('Frequency [Hz]');

