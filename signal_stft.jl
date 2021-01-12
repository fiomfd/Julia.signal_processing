# shignal_stft.jl
##########################
fs=Int(1000);
Ts=1/fs;
N=Int(12000);
T=N/fs;
#t=(0:Int(N-1))*Ts;

x=zeros(Complex, N)
for m=1:N
    x[m]=sin(2*Ts*(m-1)*(Ts*(m-1)-3)*(Ts*(m-1)-6)*(Ts*(m-1)-9)*(Ts*(m-1)-12))
end

# setting integers
L=256;
M1=128;
M=2*M1;
N0=250;
n0=1+floor(Int, (N-L)/(L-N0));
s=L-N0;

using DSP, FFTW

# window
h=hanning(L);

using LinearAlgebra

# the N-th root of unity
omega=exp(2*pi*im/M);

# sqrt(N) times the Fourier matrix
F=zeros(Complex, M, M);
for m=1:M
    for n=1:M
        F[m,n]=omega^(-(m-1)*(n-1))
    end
end

# sampling and locarization
Y=zeros(Complex, L,n0);
for l=1:L
    for j=1:n0
        Y[l,j]=h[l]*x[s*(j-1)+l]
    end
end

# discrete short-time Fourier transform
Z=F*Y;

# spectrogram
X=zeros(Float64, M1+1,n0);
for k=1:M1+1
    for j=1:n0
        X[k,j]=abs.(Y[M1+2-k,j])
    end
end

using Plots


# plotting
plot(X, 
    title="Spectrogram of x(t)=sin(2t(t-3)(t-6)(t-9)(t-12))", 
    st=:heatmap, 
    color=:cool, 
    xaxis="n = 1, 1, 2, ..., n0-1", 
    yaxis="k = 0, 1, 2, ..., M/2")
plot!(xticks = ([0 1000 1957;], [0 1000 1957]))
plot!(yticks = ([0 50 100;], [0 50 100]))
savefig("signal_stft.png")






