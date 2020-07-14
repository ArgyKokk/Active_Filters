% Argyrios Kokkinis
% 8459
% Band elimination filter
% Inverse Chebyshev

a1 = 8;
a2 = 4;
a3 = 5;
a4 = 9;

%specs
f0 = 1.8*1000;
f1 = 1200 + 25*(9-a4);
f2 = f0^2/f1;
D = (f0^2-f1^2)/(1.8*f1);
f3 = (-D+sqrt(D^2+4*f0^2))/2;
f4 = f0^2/f3;
amin = 30 - a3;
amax = 0.5 + a4/18;

%pass band from 0 to 2ðf1 and 2ðf2 to inf
%stop band from 2ðf3 to 2pf4
%pass band attenuation less than amin
%stop band attenuation bigger than amax

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

bw = w2 - w1;
Wp = 1;
Ws = bw/(w4-w3);
qc = w0/bw;

%filter's order
n = acosh((((10^(amin/10)-1))/((10^(amax/10)-1)))^(1/2))/(acosh(Ws));
n1 = n;
n = ceil(n);
e = 1/((10^(amin/10)-1)^(1/2));
a = (1/n)*asinh(1/e);

%half-power frequency
whp = 1/(cosh((1/n)*acosh(1/e)));

%Butterworth angles for n=4
% a = +- 67.5 , +- 22.5
a1 = 22.5;
a2 = -22.5;
a3 = 67.5;
a4 = -67.5;

%poles
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

%find inverse chebyshev poles
W012_ICH = 1/W01;
W034_ICH = 1/W04;

W012_n = W012_ICH*Ws;
W034_n = W034_ICH*Ws;

%find zeros
z1 = sec(1*pi/(2*n));
z3 = sec(3*pi/(2*n));

z1n = z1*Ws;
z3n = z3*Ws;

%reverse poles
W0_1_2 = 1/W012_n; %for Q12
W0_3_4 = 1/W034_n; %for Q34

%reverse zeros
Wz12 = 1/z1n;
Wz34 = 1/z3n;

%poles for high-pass
S12 = - (W0_1_2)/(2*Q12);
W12 = sqrt(W0_1_2^2-S12^2);

S34 = - (W0_3_4)/(2*Q34);
W34 = sqrt(W0_3_4^2-S34^2);

%Geffe transform for poles and zeros

%pole S12 +- W12j
S_12 = abs(S12);
W_12 = abs(W12);
C = S_12^2 + W_12^2;
D = 2*S_12/qc;
E = 4 + C/(qc^2);
G = (E^2-4*D^2)^(1/2);
Q_12 = (1/D)*sqrt((1/2)*(E+G));
k12 = S_12*Q_12/qc;
w = k12 + sqrt(k12^2-1);
w02 = w*w0;
w01= w0/w;

%pole S34 +- W34;
S_34 = abs(S34);
W_34 = abs(W34);
C = S_34^2 + W_34^2;
D = 2*S_34/qc;
E = 4 + C/(qc^2);
G = (E^2-4*D^2)^(1/2);
Q_34 = (1/D)*sqrt((1/2)*(E+G));
k34 = S_34*Q_34/qc;
w = k34 + sqrt(k34^2-1);
w04 = w*w0;
w03= w0/w;

%zero Wz12
K = 2 + Wz12^2/qc^2;
x = (K + sqrt(K^2-4))/2;
wz1 = w0*sqrt(x);
wz2 = w0/sqrt(x);

%zero Wz34
K = 2 + Wz34^2/qc^2;
x = (K + sqrt(K^2-4))/2;
wz3 = w0*sqrt(x);
wz4 = w0/sqrt(x);

%transfer functions using 7,21 - 7,23 for LPN and HPN

%unit 1 LPN
omega_z1 = wz1/w01;
R1 = 1;
R4 = 1;
C = 1/(2*Q_12);
R2 = 4*Q_12^2;
R3 = omega_z1^2/(2*Q_12^2);
R5 = (4*Q_12^2)/(omega_z1^2-1);

kf = w01;
Capac = 10^(-8);
km1 = C/(kf*Capac);
R11 = R1*km1;
R14 = R4*km1;
R13 = R3*km1;
R12 = R2*km1;
R15 = R5*km1;
gain1 = 1/(R3+1);

%unit 2 HPN
k21 = (w02/wz2)^2-1;
omega_z2 = wz2/w02;
R1 = 1;
R3 = 1;
R2 = Q_12^2*(k21+2)^2;
R4 = Q_12^2*(k21+2);
C = 1/(Q_12*(2+k21));
C1 = k21*C;

kf=w02;
Capac = 10^(-8);
km2 = C/(kf*Capac);
R21 = km2*R1;
R22 = km2*R2;
R23 = km2*R3;
R24 = km2*R4;
C21 = C1/(km2*kf);
gain2 = ((2+k21)*Q_12^2)/(1+(2+k21)*Q_12^2);
g2 = gain2*(w02/wz2)^2;

%unit 3 LPN
omega_z3 = wz3/w03;
R1 = 1;
R4 = 1;
C = 1/(2*Q_34);
R2 = 4*Q_34^2;
R3 = omega_z3^2/(2*Q_34^2);
R5 = (4*Q_34^2)/(omega_z3^2-1);

kf = w03;
Capac = 10^(-8);
km3 = C/(kf*Capac);
R31 = R1*km3;
R34 = R4*km3;
R33 = R3*km3;
R32 = R2*km3;
R35 = R5*km3;
gain3 = 1/(R3+1); 

%unit 4 HPN
k41 = (w04/wz4)^2-1;
omega_z4 = wz4/w04;
R1 = 1;
R3 = 1;
R2 = Q_34^2*(k41+2)^2;
R4 = Q_34^2*(k41+2);
C = 1/(Q_34*(2+k41));
C1 = k41*C;

kf=w04;
Capac = 10^(-8);
km4 = C/(kf*Capac);
R41 = km4*R1;
R42 = km4*R2;
R43 = km4*R3;
R44 = km4*R4;
C41 = C1/(km4*kf);
gain4 = ((2+k41)*Q_34^2)/(1+(2+k41)*Q_34^2);
g4 = gain4*(w04/wz4)^2;

%transfer functions
s = tf('s');
T1 = gain1*(s^2+wz1^2)/(s^2+(w01/Q_12)*s+w01^2);
T2 = g2*(s^2+wz2^2)/(s^2+(w02/Q_12)*s+w02^2);
T3 = gain3*(s^2+wz3^2)/(s^2+(w03/Q_34)*s+w03^2);
T4 = g4*(s^2+wz4^2)/(s^2+(w04/Q_34)*s+w04^2);

TF=T1*T2*T3*T4;

%gain 0 db
gain = 1;
gain_now = 10^0.5;
K = gain/gain_now;
H = K*TF;

%The gain is close to 10 db (10.35dB) 

%input signal vin = 0.8*(cos(w0-(w0-w3)/2)*t) +
%cos((w0+(w0+w3)/2)*t)+cos(0.5*w1*t)+ 0.8*cos(2.4*w2*t)+ 0.4*cos(3.5*w2*t)

ws1 = w0 - (w0-w3)/2;
ws2 = w0 + (w0+w3)/2;
ws3 = 0.5*w1;
ws4 = 2.4*w2;
ws5 = 3.5*w2;

%signal's duration
T = 1/200;
Total_time = 10*T;
f_s = 20000;
t=0:1/f_s:Total_time-1/f_s;

%input signal
vin = 0.8*cos(ws1*t)+cos(ws2*t)+cos(ws3*t)+0.8*cos(ws4*t)+0.4*cos(ws5*t);

%inputs frequencies
fs1 = ws1/(2*pi);
fs2 = ws2/(2*pi);
fs3 = ws3/(2*pi);
fs4 = ws4/(2*pi);
fs5 = ws5/(2*pi);

%plot transient analysis
out=lsim(TF,vin,t);
plot(t,vin,'green')
hold on;
plot(t,out,'red')
title('Transient Analysis');
ylabel('Magnitude');

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



