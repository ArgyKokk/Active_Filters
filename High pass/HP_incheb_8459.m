% Argyrios Kokkinis
% 8459
% High pass filter
% Inverse Chebyshev

a1=8;
a2=4;
a3=5;
a4=9;

%specs
m=0;
fp= (3+m)*1000;
fs= fp/1.8;
amin=25 + a3*(3/9);
amax=0.4+ (a4/36);


ws = 2*pi*fs;
wp = 2*pi*fp;


%stop band from 0 to 2ðfs
%pass band from 2ðfp to inf

Wp=1;
Ws= wp/ws;

%transform frequencies for ICH filter 
Wp=Wp/Ws;
Ws=1;

%filter's order
n = acosh((((10^(amin/10)-1))/((10^(amax/10)-1)))^(1/2))/(acosh(1/Wp));
n1 = n;
n = ceil(n);

e = 1/((10^(amin/10)-1)^(1/2));
a = (1/n)*asinh(1/e);

%half-power frequency
whp = 1/(cosh((1/n)*acosh(1/e)));

%butterworth angles for n=4
%a = +- 67.5 , +- 22.5

a1 = 22.5;
a2 = -22.5;
a3 = 67.5;
a4 = -67.5;

%poles chebyshev
p1 = -sinh(a)*cosd(a1) + cosh(a)*sind(a1)*1i;
p2 = -sinh(a)*cosd(a2) + cosh(a)*sind(a2)*1i;
p3 = -sinh(a)*cosd(a3) + cosh(a)*sind(a3)*1i;
p4 = -sinh(a)*cosd(a4) + cosh(a)*sind(a4)*1i;

%find Ùïê
W01 = sqrt(real(p1)^2+imag(p1)^2);
W02 = W01;
W03 = sqrt(real(p3)^2+imag(p3)^2);
W04 = W03;

%find Qs
Q1 = sqrt(real(p1)^2 + imag(p1)^2)/(2*abs(real(p1)));
Q2 = Q1;
Q3 = sqrt(real(p3)^2 + imag(p3)^2)/(2*abs(real(p3)));
Q4 = Q3;
Q12 = Q1;
Q34 = Q3;

%reverse chebyshev poles for ICH
p1_r = 1/p1;
p2_r = 1/p2;
p3_r = 1/p3;
p4_r = 1/p4;

W01r = sqrt(real(p1_r)^2+imag(p1_r)^2);
W02r = W01r;
W03r = sqrt(real(p3_r)^2+imag(p3_r)^2);
W04r = W03r;
W012r= W01r;
W034r = W03r;

%find zeros
Wz1 = sec(1*pi/(2*n));
Wz2 = sec(3*pi/(2*n));
Wz12 = Wz1;
Wz34 = Wz2;

%reverse poles
w012 = wp/(2*W012r);
w034 = wp/(2*W034r);
wz12 = wp/(2*Wz12);
wz34 = wp/(2*Wz34);

% module one : pole w012 , zero wz34
% module two : pole w034 , zero wz12 

Qm1 = 1/(1-(wz34/w012)^2);
Qm2 = 1/(1-(wz12/w034)^2);

%transfer functions
s = tf('s');
% I wannt 10dB gain at high frequencies
% I set 5db gain for each module
gain = 1.78;
T1 = gain *(s^2+(wz34)^2)/(s^2+w012/Q12*s+w012^2);
T2 = gain *(s^2+(wz12)^2)/(s^2+w034/Q34*s+w034^2);
TF = T1*T2;

%for 0 gain
gain1 = 1;
T1_0 = gain1 *(s^2+(wz34)^2)/(s^2+w012/Q12*s+w012^2);
T2_0 = gain1 *(s^2+(wz12)^2)/(s^2+w034/Q34*s+w034^2);
T_0 = T1_0*T2_0;

% set R4 to 100Ohm
R4 = 100;

%module 1
BoctorHighPass(wz34,w012,Q12,R4,10^(-6));
module1 = circuit;

%module 2
BoctorHighPass(wz12,w034,Q34,R4,10^(-6));
module2 = circuit;

%input signal vin=
%cos(0.5*ws*t)+0.6*cos(0.8*ws*t)+cos(1.2*wp*t)+0.8*cos(2.4*wp*t)+0.4*cos(3.5*wp*t)

ws1 = 0.5*ws;
ws2 = 0.8*ws;
ws3 = 1.2*wp;
ws4 = 2.4*wp;
ws5 = 3.5*wp;

%signal's duration
T = 1/100;
Total_time = 10*T;
f_s = 20000;
t=0:1/f_s:Total_time-1/f_s;

%input signal
vin = cos(ws1*t)+0.6*cos(ws2*t)+cos(ws3*t)+0.8*cos(ws4*t)+0.4*cos(ws5*t);

%inputs frequencies
fs1 = ws1/(2*pi);
fs2 = ws2/(2*pi);
fs3 = ws3/(2*pi);
fs4 = ws4/(2*pi);
fs5 = ws5/(2*pi);

%plot transient analysis
out=lsim(TF,vin,t);
plot(t,vin,'red')
hold on;
plot(t,out,'green')
title('Transient Analysis');
ylabel('Magnitude');
hold off;

%input signal fourier
N = 100;
n=2^nextpow2(N);
xfourier= fft(vin,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
plot(f,p1,'red')
title('input spectrum');
xlabel('Frequency [Hz]');
hold on

%output signal fourier
N = 100;
n=2^nextpow2(N);
xfourier= fft(out,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
plot(f,p1,'green')
title('output spectrum');
xlabel('Frequency [Hz]');