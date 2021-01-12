# shignal_dft_shifted.jl
##########################
fs=Int(1000);
Ts=1/fs;
N=Int(12000);
T=N/fs;
#t=(0:Int(N-1))*Ts;

x=zeros(Float64, N)
for m=1:N
    x[m]=sin(2*Ts*(m-1)*(Ts*(m-1)-3)*(Ts*(m-1)-6)*(Ts*(m-1)-9)*(Ts*(m-1)-12))
end

using FFTW, DSP 

hatx=fft(x);
k = (0:length(hatx)-1)*50/length(hatx);
G=broadcast(abs, hatx);
kshift = (-length(x)/2:length(x)/2-1)*(50/length(x));
Gshift = fftshift(G);

using Plots

plot(kshift, Gshift, 
    xaxis="k = -N/2, ..., 0, ..., N/2", 
    yaxis="|F[x](k)|", 
    title="DFT of x(t)=sin(2t(t-3)(t-6)(t-9)(t-12))", 
    grid=false, 
    legend=false)
plot!(xticks = ([-25 0 25;], [-6000 0 6000]))
plot!(yticks = ([0 200 400;], [0 200 400]))
savefig("signal_dft_shifted.png")
