# signal.jl
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

using Plots


plot(x, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude x[n]", 
    title="Discrete time signal x(t)=sin(2t(t-3)(t-6)(t-9)(t-12))", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 5000 10000;], [0, 5000, 10000]))
plot!(yticks = ([-1 1;], [-1, 1]))

plot(x, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Audio Signal x(t)=sin(2t(t-3)(t-6)(t-9)(t-12))", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 5000 10000;], [0, 5, 10]))
plot!(yticks = ([-1 1;], [-1, 1]))
savefig("signal.png") 
