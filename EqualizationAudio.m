% Trabalho 2 - Eletrônica Aplicada
% Lorenzo Pereira Piccoli Xavier
% Filtagem de sinal de audio ruidoso baseado na
% Equalização de Sinais em Bandas de Freqüência
% Filtros passa bandas elípticos
% Transformação bilinear
% Coeficientes IIR 
% Input: 'shrek.wav'
% Output: 'shrek.wav' equalizado

clear all;
close all;
format longEng

function ybp = bp(sinal,fs,Rp,Rs,fl,fh)

Wp = [fl fh*0.9]/fs;
Ws = [fl*0.9 fh]/fs;

[n,Wp1] = ellipord(Wp,Ws,Rp,Rs);

ordem = n
[z,p,k] = ellipap(ordem,Rp,Rs);
[b0,a0] = zp2tf(z,p,k);

%if fl>=2200
%    fineBw = 1.1
%end

Wo = 2*pi*sqrt(fl*fh); % center frequency
Bw = 2*pi*(fh-fl); % bandwidth

[num,den] = lp2bp(b0,a0,Wo,Bw);

fs = 44.1e3; % para evitar o aliasing

%[num1,den1] = bilinear(num,den,fs); %non optimum
% better reponse
[num1,den1] = impinvar(num,den,fs);

%figure
%title('FREQS')
%freqs(num,den)

% figure
% title('FREQZ')
% freqz(tf2sos(num1,den1),20048,fs)

%similar to freqz
omega = linspace(0,fs/2,20e3)
H =  filt(num1,den1,1/fs)

%{
figure
H_omega = squeeze(freqresp(H,omega))
plot(omega/(2*pi),abs(H_omega))
grid on
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
%}


b = cell2mat(H.Numerator);
a = cell2mat(H.Denominator);
% O último índice é o de menor ordem z^-6
% O primeiro índice é o de maior ordem z^0

% Resposta ao sinal senoidal de 500 Hz
xn = zeros(1,length(a));
yn = zeros(1,length(b));
y500 = zeros(1,length(sinal));
for n=1:length(sinal)
     xn(1)=sinal(n);
     
     y500(n)=xn*transpose(b)-yn(2:length(a))*transpose(a(2:length(a)));
     yn(1)=y500(n);
     xn(2:length(b))=xn(1:length(b)-1);
     yn(2:length(a))=yn(1:length(a)-1);
end

zz=20*log10(abs(fft(y500))); 
zzi=20*log10(abs(fft(sinal))); 
samples = length(sinal);
w=0:fs/samples:fs-fs/samples;
% figure
% plot(w,zzi,':k',w,zz,'b');
% title('Mag resposta filtrada');
% figure
% plot(w,angle(fft(y500)));
% title('Phase');

ybp=y500;

end

function ylp = lp(sinal,fs,Rp,Rs,fl)

Wp = fl*0.9/fs;
Ws = fl/fs;

[n,Wp] = ellipord(Wp,Ws,Rp,Rs);

ordem = n
[z,p,k] = ellip(ordem,Rp,Rs,fl*2*pi,'s');
[num,den] = zp2tf(z,p,k);

%figure
%freqs(num,den)

fs = 44.1e3; % para evitar o aliasing

[num1,den1] = bilinear(num,den,fs);

% figure
% freqz(tf2sos(num1,den1),20048,fs)

%similar to freqz
H =  filt(num1,den1,1/fs);

%{
figure
H_omega = squeeze(freqresp(H,omega))
plot(omega/(2*pi),abs(H_omega))
grid on
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
%}


b = cell2mat(H.Numerator);
a = cell2mat(H.Denominator);
% O último índice é o de menor ordem z^-6
% O primeiro índice é o de maior ordem z^0

% Resposta ao sinal senoidal de 500 Hz
xn = zeros(1,length(a));
yn = zeros(1,length(b));
y500 = zeros(1,length(sinal));
for n=1:length(sinal)
     xn(1)=sinal(n);
     
     y500(n)=xn*transpose(b)-yn(2:length(a))*transpose(a(2:length(a)));
     yn(1)=y500(n);
     xn(2:length(b))=xn(1:length(b)-1);
     yn(2:length(a))=yn(1:length(a)-1);
end

zz=20*log10(abs(fft(y500))); 
zzi=20*log10(abs(fft(sinal))); 
samples = length(sinal);
w=0:fs/samples:fs-fs/samples;
% figure
% plot(w,zzi,':k',w,zz,'b');
% title('Mag resposta filtrada');
% figure
% plot(w,angle(fft(y500)));
% title('Phase');

ylp=y500;

end


% Run only one
[signal, fs] = audioread('shrek.wav');
%fs = 44.1e3;
%signal = eye(1,length(signal));

% Período de amostragem
T=1/fs;

%Declara as faixas dos filtros passa bandas elípticos
faixas = [ 
    500 800;
    900 2000;
    2200 4000;
    4200 7000;
    7200 10000;
    ];
%Declara os pesos dos filtros passa bandas elípticos
pesos = [ 
    1;
    1;
    1;
    0.75;
    0.5;
    ];

%Inicia vetor que acumulará as respostas ponderadas dos filtros
sinal_filtrado = zeros(1,length(signal));

%Filtro passa baixa elíptico IIR
Rp = 0.5;
Rs = 40;
y = lp(signal,fs,Rp,Rs,400);
sinal_filtrado = sinal_filtrado + y*2;

%Filtros passa banda elíticos IIR
for n=1:length(pesos)
    fl = faixas(n,1);
    fh = faixas(n,2);
    Rp = 1;
    Rs = 40;
    y = bp(signal,fs,Rp,Rs,fl,fh);
    sinal_filtrado = sinal_filtrado + y*pesos(n);
end

%Toca o som original
%sound(signal,fs)

%Toca o som filtrado
sound(sinal_filtrado,fs)

%Plot the sound original and filtred ones
figure
tempo=0:T:(length(signal)-1)*T;
plot(tempo,signal,'-r',tempo,sinal_filtrado,':k')


zz=20*log10(abs(fft(sinal_filtrado))); 
zzi=20*log10(abs(fft(signal))); 
samples = length(signal);
w=0:fs/samples:fs-fs/samples;
figure
plot(w,zzi,':k',w,zz,'b');
title('Mag resposta filtrada');