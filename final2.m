clc; close all; clear all;

%generate message signal
[m1, fs]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_BBCArabic2.wav');
[m2, fs]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_FM9090.wav');
[m3, fs]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_QuranPalestine.wav');
[m4, fs]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_RussianVoice.wav');
[m5, fs]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_SkyNewsArabia.wav');
[m6, fs]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_WRNArabic.wav');

L = 36;
fs = L*fs;

m1 = m1(:,1)+m1(:,2);
m2 = m2(:,1)+m2(:,2);
m3 = m3(:,1)+m3(:,2);
m4 = m4(:,1)+m4(:,2);
m5 = m5(:,1)+m5(:,2);
m6 = m6(:,1)+m6(:,2);

m1 = interp(m1, L);
m2 = interp(m2, L);
m3 = interp(m3, L);
m4 = interp(m4, L);
m5 = interp(m5, L);
m6 = interp(m6, L);

n1=length(m1);
n2=length(m2);
n3=length(m3);
n4=length(m4);
n5=length(m5);
n6=length(m6);

N_arr=[n1 n2 n3 n4 n5 n6];
NN=max(N_arr);

m1_padding = zeros(NN,1,'double');
m2_padding = zeros(NN,1,'double');
m3_padding = zeros(NN,1,'double');
m4_padding = zeros(NN,1,'double');
m5_padding = zeros(NN,1,'double');
m6_padding = zeros(NN,1,'double');

for i = 1: n1
    m1_padding(i)= m1(i); 
end
for i = 1: n2
    m2_padding(i)= m2(i); 
end
for i = 1: n3
    m3_padding(i)= m3(i); 
end
for i = 1: n4
    m4_padding(i)= m4(i); 
end
for i = 1: n5
    m5_padding(i)= m5(i); 
end
for i = 1: n6
    m6_padding(i)= m6(i); 
end

Ts = 1/fs;
t = (0:Ts:(NN-1)*Ts).';
f= (-NN/2: (NN/2)-1)*(fs/NN).';

M1 = abs(fftshift(fft(m1_padding)));
M2 = abs(fftshift(fft(m2_padding)));
M3 = abs(fftshift(fft(m3_padding)));
M4 = abs(fftshift(fft(m4_padding)));
M5 = abs(fftshift(fft(m5_padding)));
M6 = abs(fftshift(fft(m6_padding)));

%generate carrier signal
wn1 = 2*pi*100*10^3;
ct1 = cos(wn1*t);
cf1 = abs(fftshift(fft(ct1)));

wn2 = 2*pi*(100 + 1*50)*10^3;
ct2 = cos(wn2*t);
cf2 = abs(fftshift(fft(ct2)));

wn3 = 2*pi*(100 + 2*50)*10^3;
ct3= cos(wn3*t);

wn4 = 2*pi*(100 + 3*50)*10^3;
ct4= cos(wn4*t);

wn5 = 2*pi*(100 + 4*50)*10^3;
ct5= cos(wn5*t);

wn6 = 2*pi*(100 + 5*50)*10^3;
ct6= cos(wn6*t);

%modulation signal
st1 = m1_padding.*ct1;
sf1 = abs(fftshift(fft(st1)));

st2 = m2_padding.*ct2;
sf2 = abs(fftshift(fft(st2)));

st3 = m3_padding.*ct3;
sf3 = abs(fftshift(fft(st3)));

st4 = m4_padding.*ct4;
sf4 = abs(fftshift(fft(st4)));

st5 = m5_padding.*ct5;
sf5 = abs(fftshift(fft(st5)));

st6 = m6_padding.*ct6;
sf6 = abs(fftshift(fft(st6)));

st = st1 + st2 + st3 + st4 + st5 + st6;
sf = abs(fftshift(fft(st)));

%% Mixer and oscillator
wif = 2*pi*25*10^3;
wbw = 2*pi*10000;
 
yt1_mixer11= st.*cos((wn1+wif)*t);
yf1_mixer11= abs(fftshift(fft(yt1_mixer11)));

yt2_mixer11= st.*cos((wn2+wif)*t);
yf2_mixer11= abs(fftshift(fft(yt2_mixer11)));

yt3_mixer11= st.*cos((wn3+wif)*t); 
yf3_mixer11= abs(fftshift(fft(yt3_mixer11)));

yt4_mixer11= st.*cos((wn4+wif)*t); 
yf4_mixer11= abs(fftshift(fft(yt4_mixer11)));

yt5_mixer11= st.*cos((wn5+wif)*t); 
yf5_mixer11= abs(fftshift(fft(yt5_mixer11)));

yt6_mixer11= st.*cos((wn6+wif)*t); 
yf6_mixer11= abs(fftshift(fft(yt6_mixer11)));

%IF stage
bnd9 =[wif-wbw wif+wbw]./(pi*fs);
num9 = fir1(256,bnd9);

yt1_if2= filter(num9,1,yt1_mixer11);
yf1_if2= abs(fftshift(fft(yt1_if2)));

yt2_if2= filter(num9,1,yt2_mixer11);
yf2_if2= abs(fftshift(fft(yt2_if2)));

yt3_if2= filter(num9,1,yt3_mixer11);
yf3_if2= abs(fftshift(fft(yt3_if2)));

yt4_if2= filter(num9,1,yt4_mixer11);
yf4_if2= abs(fftshift(fft(yt4_if2)));

yt5_if2= filter(num9,1,yt5_mixer11);
yf5_if2= abs(fftshift(fft(yt5_if2)));

yt6_if2= filter(num9,1,yt6_mixer11);
yf6_if2= abs(fftshift(fft(yt6_if2)));

%Baseband detection
yt1_mixer22=yt1_if2.*cos(wif*t);
yf1_mixer22= abs(fftshift(fft(yt1_mixer22)));

yt2_mixer22=yt2_if2.*cos(wif*t); 
yf2_mixer22= abs(fftshift(fft(yt2_mixer22)));

yt3_mixer22=yt3_if2.*cos(wif*t); 
yf3_mixer22= abs(fftshift(fft(yt3_mixer22)));

yt4_mixer22=yt4_if2.*cos(wif*t); 
yf4_mixer22= abs(fftshift(fft(yt4_mixer22)));

yt5_mixer22=yt5_if2.*cos(wif*t); 
yf5_mixer22= abs(fftshift(fft(yt5_mixer22)));

yt6_mixer22=yt6_if2.*cos(wif*t); 
yf6_mixer22= abs(fftshift(fft(yt6_mixer22)));

bnd10 = wbw/(pi*fs);
num10 = fir1(256,bnd10,'low');
den = [1];

yt1_lpf2 = filter(num10,den,yt1_mixer22);
yf1_lpf2= abs(fftshift(fft(yt1_lpf2)));

yt2_lpf2 = filter(num10,den,yt2_mixer22);
yf2_lpf2= abs(fftshift(fft(yt2_lpf2)));

yt3_lpf2 = filter(num10,den,yt3_mixer22);
yf3_lpf2= abs(fftshift(fft(yt3_lpf2)));

yt4_lpf2 = filter(num10,den,yt4_mixer22);
yf4_lpf2= abs(fftshift(fft(yt4_lpf2)));

yt5_lpf2 = filter(num10,den,yt5_mixer22);
yf5_lpf2= abs(fftshift(fft(yt5_lpf2)));

yt6_lpf2 = filter(num10,den,yt6_mixer22);
yf6_lpf2= abs(fftshift(fft(yt6_lpf2)));

audiowrite('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_BBCArabic2_recv2.wav',yt1_lpf2,fs);
audiowrite('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_FM9090_recv2.wav',yt2_lpf2,fs);
audiowrite('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_QuranPalestine_recv2.wav',yt3_lpf2,fs);
audiowrite('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_RussianVoice_recv2.wav',yt4_lpf2,fs);
audiowrite('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_SkyNewsArabia_recv2.wav',yt5_lpf2,fs);
audiowrite('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_WRNArabic_recv2.wav',yt6_lpf2,fs);

[recv1,fs_new1]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_BBCArabic2_recv2.wav');
[recv2,fs_new2]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_FM9090_recv2.wav');
[recv3,fs_new3]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_QuranPalestine_recv2.wav');
[recv4,fs_new4]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_RussianVoice_recv2.wav');
[recv5,fs_new5]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_SkyNewsArabia_recv2.wav');
[recv6,fs_new6]= audioread('C:\Users\nessm\OneDrive\Documents\MATLAB/Short_WRNArabic_recv2.wav');

recv1 = downsample(recv1,L);
recv2 = downsample(recv2,L);
recv3 = downsample(recv3,L);
recv4 = downsample(recv4,L);
recv5 = downsample(recv5,L);
recv6 = downsample(recv6,L);

duration = 20;
pause(duration)                                 
sound(recv1,fs_new1/L)
pause(duration)                                 
sound(recv2,fs_new2/L)
pause(duration)                                 
sound(recv3,fs_new3/L)
pause(duration)                                 
sound(recv4,fs_new4/L)
pause(duration)                                 
sound(recv5,fs_new5/L)
pause(duration)                                 
sound(recv6,fs_new6/L)
%% drawing 2
figure(7)
subplot(3,2,1)
plot(f,yf1_mixer11);
xlabel('frequency');ylabel('freq Amplitude');title('RF Mixer signal 1');

subplot(3,2,2)
plot(f,yf2_mixer11);
xlabel('frequency');ylabel('freq Amplitude');title('RF Mixer signal 2');

subplot(3,2,3)
plot(f,yf3_mixer11);
xlabel('frequency');ylabel('freq Amplitude');title('RF Mixer signal 3');

subplot(3,2,4)
plot(f,yf4_mixer11);
xlabel('frequency');ylabel('freq Amplitude');title('RF Mixer signal 4');

subplot(3,2,5)
plot(f,yf5_mixer11);
xlabel('frequency');ylabel('freq Amplitude');title('RF Mixer signal 5');

subplot(3,2,6)
plot(f,yf6_mixer11);
xlabel('frequency');ylabel('freq Amplitude');title('RF Mixer signal 6');

figure(8)
subplot(3,2,1)
plot(f,yf1_if2);
xlabel('frequency');ylabel('freq Amplitude');title('IF signal 1');

subplot(3,2,2)
plot(f,yf2_if2);
xlabel('frequency');ylabel('freq Amplitude');title('IF signal 2');

subplot(3,2,3)
plot(f,yf3_if2);
xlabel('frequency');ylabel('freq Amplitude');title('IF signal 3');

subplot(3,2,4)
plot(f,yf4_if2);
xlabel('frequency');ylabel('freq Amplitude');title('IF signal 4');

subplot(3,2,5)
plot(f,yf5_if2);
xlabel('frequency');ylabel('freq Amplitude');title('IF signal 5');

subplot(3,2,6)
plot(f,yf6_if2);
xlabel('frequency');ylabel('freq Amplitude');title('IF signal 6');

figure(9)
subplot(3,2,1)
plot(f,yf1_mixer22);
xlabel('frequency');ylabel('freq Amplitude');title('IF Mixer signal 1');

subplot(3,2,2)
plot(f,yf2_mixer22);
xlabel('frequency');ylabel('freq Amplitude');title('IF Mixer signal 2');

subplot(3,2,3)
plot(f,yf3_mixer22);
xlabel('frequency');ylabel('freq Amplitude');title('IF Mixer signal 3');

subplot(3,2,4)
plot(f,yf4_mixer22);
xlabel('frequency');ylabel('freq Amplitude');title('IF Mixer signal 4');

subplot(3,2,5)
plot(f,yf5_mixer22);
xlabel('frequency');ylabel('freq Amplitude');title('IF Mixer signal 5');

subplot(3,2,6)
plot(f,yf6_mixer22);
xlabel('frequency');ylabel('freq Amplitude');title('IF Mixer signal 6');


figure(10)
subplot(3,2,1)
plot(f,yf1_lpf2);
xlabel('frequency');ylabel('freq Amplitude');title('LPF signal 1');

subplot(3,2,2)
plot(f,yf2_lpf2);
xlabel('frequency');ylabel('freq Amplitude');title('LPF signal 2');

subplot(3,2,3)
plot(f,yf3_lpf2);
xlabel('frequency');ylabel('freq Amplitude');title('LPF signal 3');

subplot(3,2,4)
plot(f,yf4_lpf2);
xlabel('frequency');ylabel('freq Amplitude');title('LPF signal 4');

subplot(3,2,5)
plot(f,yf5_lpf2);
xlabel('frequency');ylabel('freq Amplitude');title('LPF signal 5');

subplot(3,2,6)
plot(f,yf6_lpf2);
xlabel('frequency');ylabel('freq Amplitude');title('LPF signal 6');
