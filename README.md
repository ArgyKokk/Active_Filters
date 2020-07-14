# Active_Filters
Design of four different analog filters with different specifications each.
Low Pass filter, Band Pass filter, Band Elimination filter , High Pass filter
All filters were designed in MATLAB and implemented on Multisim.<br>
course: Design of active filters, AUTH, 2020

## Low Pass filter
Design of a Butterworth low pass  filter with specifications regarding the frequency and attenuation
```MATLAB
fp = 3 kHz , fs = 6 kHz
amin = 0.8 dB , amax = 19.5 dB
```
## Band Pass filter
Design of a Chebyshev band pass filter with specifications regarding the frequency and attenuation
```MATLAB
fo = 650 Hz , f1 = 525 Hz , f2 = 804.7 Hz , f3 = 403.5 Hz , f4 = 1.04 kHz
amin = 0.5 dB , amax = 36.5 dB
```
## Band elimination filter
Design of an Inverse Chebyshev band elimination filter with specifications regarding the frequency and attenuation
```MATLAB
fo = 1.8 kHz , f1 = 1.2 kHz , f2 = 2.7 kHz , f3 = 1.43 kHz , f4 = 2.264 kHz 
amin = 1 dB , amax = 25 dB
```
## High Pass filter
Design of an Inverse Chebyshev high pass with specifications regarding the frequency and attenuation
```MATLAB
fp = 3 kHz , fs = 1.666 kHz
amin = 0.65 dB , amax = 26.66 dB
```
