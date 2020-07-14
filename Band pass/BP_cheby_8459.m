% Argyrios Kokkinis
% 8459
% Band pass filter
% Chebyshev

a1=8;
a2=4;
a3=5;
a4=9;

%specs
f0 = 0.65*1000;
f1 = 400 + 25*a3;
f2 = f0^2/f1;
D = 2.3 * (f0^2-f1^2)/f1;
f3 = (-D + sqrt(D^2+4*f0^2))/2;
f4 = f0^2/f3;
amin = 27.5+a4;
amax = 0.5 + (a3-5)/10;

%stop band from 0 to 2ðf3 , and 2ðf4 to inf 
%(attenuation greater than 36.5dB)
%pass band from 2ðf1 to 2ðf2
%(attenuation less than 0.5dB)

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;
bw = w2 - w1; 

Wp = 1;
Ws = (w4-w3)/(w2-w1);
%filter's order
n = acosh(((10^(amin/10)-1)/(10^(amax/10)-1))^(1/2))/acosh(Ws);
n1=n;
n = ceil(n);
e = (10^(amax/10)-1)^(1/2);
a = (1/n)*asinh(1/e);

%half-power frequency
whp = cosh((1/n)*acosh(1/e));

%Butterworth angles for n=5
%0 , +- 36 , +- 72
a1 = 0;
a2 = 36;
a3 = -36;
a4 = 72;
a5 = -72;

%poles
p1 = -sinh(a)*cosd(a1) + cosh(a)*sind(a1)*1i;
p2 = -sinh(a)*cosd(a2) + cosh(a)*sind(a2)*1i;
p3 = -sinh(a)*cosd(a3) + cosh(a)*sind(a3)*1i;
p4 = -sinh(a)*cosd(a4) + cosh(a)*sind(a4)*1i;
p5 = -sinh(a)*cosd(a5) + cosh(a)*sind(a5)*1i;

%Qs
Q1 = sqrt(real(p1)^2 + imag(p1)^2)/(2*abs(real(p1)));
Q2 = sqrt(real(p2)^2 + imag(p2)^2)/(2*abs(real(p2)));
Q3 = sqrt(real(p4)^2 + imag(p4)^2)/(2*abs(real(p4)));

%transform pole p1
S1 = abs(real(p1));
qc = w0/bw;
Q_1 = qc/S1;
psi1 = acos(1/(2*Q_1));

%transform complex poles p2,3
S2 = abs(real(p2));
Omega2 = imag(p2);
C = S2^2 + Omega2^2;
D = (2*S2)/qc;
E = 4 + C/(qc^2);
G = (E^2-4*D^2)^(1/2);
Q_2 = (1/D)*sqrt((1/2)*(E+G));
k = (S2*Q_2)/qc;
w = k + sqrt(k^2-1);
w02 = w*w0;
w01 = (1/w)*w0;

%transform complex poles p4,5
S4 = abs(real(p4));
Omega4 = imag(p4);
C = S4^2 + Omega4^2;
D = (2*S4)/qc;
E = 4 + C/(qc^2);
G = (E^2-4*D^2)^(1/2);
Q_4 = (1/D)*sqrt((1/2)*(E+G));
k = (S4*Q_4)/qc;
w = k + sqrt(k^2-1);
w04 = w*w0;
w03 = (1/w)*w0;

%Unit 1 (strategy 1), 1ìF capacitances
Capacitance = 10^(-6);
R_01 = 1;
R_02 = 4*Q_1^2;
C_01 = 1/(2*Q_1);
kf1 = w0;
km1 = C_01/(Capacitance*kf1);
R11 = km1 * R_01;
R12 = km1 * R_02;

% Transfer function unit1
s = tf('s');
H1 = - (2*Q_1*w0*s)/(s^2+(w0/Q_1)*s+w0^2);

%Unit 2 (strategy 1), 1ìF capacitances
Capacitance = 10^(-6);
R_11 = 1;
R_12 = 4*Q_2^2;
C_11 = 1/(2*Q_2);
kf2 = w01;
km2 = C_11/(Capacitance*kf2);
R21 = km2 * R_11;
R22 = km2 * R_12;

% Transfer function unit2
s = tf('s');
H2 = - (2*Q_2*w01*s)/(s^2+(w01/Q_2)*s+w01^2);

%Unit 3 (strategy 1), 1ìF capacitances
Capacitance = 10^(-6);
R_21 = 1;
R_22 = 4*Q_2^2;
C_21 = 1/(2*Q_2);
kf3 = w02;
km3 = C_21/(Capacitance*kf3);
R31 = km3 * R_21;
R32 = km3 * R_22;

% Transfer function unit3
s = tf('s');
H3 = - (2*Q_2*w02*s)/(s^2+(w02/Q_2)*s+w02^2);

%Unit 4 (strategy 1), 1ìF capacitances
Capacitance = 10^(-6);
R_31 = 1;
R_32 = 4*Q_4^2;
C_31 = 1/(2*Q_4);
kf4 = w03;
km4 = C_31/(Capacitance*kf4);
R41 = km4 * R_31;
R42 = km4 * R_32;

% Transfer function unit4
s = tf('s');
H4 = - (2*Q_4*w03*s)/(s^2+(w03/Q_4)*s+w03^2);

%Unit 5 (strategy 1), 1ìF capacitances
Capacitance = 10^(-6);
R_41 = 1;
R_42 = 4*Q_4^2;
C_41 = 1/(2*Q_4);
kf5 = w04;
km5 = C_41/(Capacitance*kf5);
R51 = km5 * R_41;
R52 = km5 * R_42;

% Transfer function unit5
s = tf('s');
H5 = - (2*Q_4*w04*s)/(s^2+(w04/Q_4)*s+w04^2);

H = H1*H2*H3*H4*H5;

% I want 5 db gain at w0 for Hbp
gain = 10 ^(0.25);
cur_gain = 10 ^(9.35);
k = gain/cur_gain;

%regulate gain for module1 to 1
%regulate gain for module2 to 0.5
%regulate gain for module3 to 0.5
%regulate gain for module4 to 2.666
%regulate gain for module5 to 2.666

%module1 - gain regulation
gain1 = 10^0;
cur_gain1 = 10^(1.915);
k1 = gain1/cur_gain1;

%attenuate the signal for gain1
Za1 = R11/k1;
Zb1 = R11/(1-k1);

%module2 - gain regulation
gain2 = 0.5;
cur_gain2 = 53.8;
k2 = gain2/cur_gain2;

%attenuate the signal for gain2
Za2 = R21/k2;
Zb2 = R21/(1-k2);

%module3 - gain regulation
gain3 = 0.5;
cur_gain3 = 53.8;
k3 = gain3/cur_gain3;

%attenuate the signal for gain3
Za3 = R31/k3;
Zb3 = R31/(1-k3);

%module4 - gain regulation
gain4 = 2.6667;
cur_gain4 = 96.73;
k4 = gain4/cur_gain4;

%attenuate the signal for gain4
Za4 = R41/k4;
Zb4 = R41/(1-k4);

%module5 - gain regulation
gain5 = 2.6667;
cur_gain5 = 96.73;
k5 = gain5/cur_gain5;

%attenuate the signal for gain4
Za5 = R51/k5;
Zb5 = R51/(1-k5);

% for 0dB gain
gain0 = 1;
k0 = gain0/cur_gain;

%attenuate the signal before the first module
%Za = R11/k;
%Zb = R11/(1-k);
k=k1*k2*k3*k4*k5;
H1 = k1*H1;
H2 = k2*H2;
H3 = k3*H3;
H4 = k4*H4;
H5 = k5*H5;
Hnew = k*H;

%input signal
%vin = cos((w0-(w0-w1)/2)*t)+
%0.8*cos((w0+(w0+w1)/3)*t)+0.8*cos(0.4*w3*t)+0.6*cos(2.5*w4*t)+0.5*cos(3*w4*t)

ws1 = w0-(w0-w1)/2;
ws2 = w0+(w0+w1)/3;
ws3 = 0.4*w3;
ws4 = 2.5*w4;
ws5 = 3*w4;

%input signal frequencies
fs1 = ws1/(2*pi);
fs2 = ws2/(2*pi);
fs3 = ws3/(2*pi);
fs4 = ws4/(2*pi);
fs5 = ws5/(2*pi);

%signal's duration
T = 1/200;
Total_time = 10*T;
f_s = 200000;
t=0:1/f_s:Total_time-1/f_s;

%input signal
vin = cos(ws1*t)+0.8*cos(ws2*t)+0.8*cos(ws3*t)+0.6*cos(ws4*t)+0.5*cos(ws5*t);

%plot transient analysis
out=lsim(Hnew,vin,t);
plot(t,vin,'red')
hold on;
plot(t,out,'green')
title('Transient Analysis');
ylabel('Magnitude');
hold off

%input signal fourier
N = 5000;
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
N = 5000;
n=2^nextpow2(N);
xfourier= fft(out,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
plot(f,p1,'green')
title('output spectrum');
xlabel('Frequency [Hz]');
