# signal_haar.jl
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

using DSP, Wavelets

y1=dwt(x, wavelet(WT.haar), 1);
y2=dwt(x, wavelet(WT.haar), 2);
y3=dwt(x, wavelet(WT.haar), 3);
y4=dwt(x, wavelet(WT.haar), 4);
y5=dwt(x, wavelet(WT.haar), 5);

y1a=zeros(N);
y2a=zeros(N);
y3a=zeros(N);
y4a=zeros(N);
y5a=zeros(N);
y1a[1:Int(N/2)]=y1[1:Int(N/2)];
y2a[1:Int(N/4)]=y2[1:Int(N/4)];
y3a[1:Int(N/8)]=y3[1:Int(N/8)];
y4a[1:Int(N/16)]=y4[1:Int(N/16)];
y5a[1:Int(N/32)]=y5[1:Int(N/32)];
y1d=y1-y1a;
y2d=y2-y2a;
y3d=y3-y3a;
y4d=y4-y4a;
y5d=y5-y5a;

x1a=idwt(y1a, wavelet(WT.haar), 1);
x1d=idwt(y1d, wavelet(WT.haar), 1);
x2a=idwt(y2a, wavelet(WT.haar), 2);
x2d=idwt(y2d, wavelet(WT.haar), 2);
x3a=idwt(y3a, wavelet(WT.haar), 3);
x3d=idwt(y3d, wavelet(WT.haar), 3);
x4a=idwt(y4a, wavelet(WT.haar), 4);
x4d=idwt(y4d, wavelet(WT.haar), 4);
x5a=idwt(y5a, wavelet(WT.haar), 5);
x5d=idwt(y5d, wavelet(WT.haar), 5);

using Plots

# x1a
p1a=plot(x1a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 1 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1d
p1d=plot(x1d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 1 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1
plot(p1a, p1d, layout = (1,2))
savefig("signal_haar_level_1.png") 

# x2a
p2a=plot(x2a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 2 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x2d
p2d=plot(x2d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 2 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1
plot(p2a, p2d, layout = (1,2))
savefig("signal_haar_level_2.png") 

# x3a
p3a=plot(x3a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 3 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x3d
p3d=plot(x3d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 3 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1
plot(p3a, p3d, layout = (1,2))
savefig("signal_haar_level_3.png") 

# x4a
p4a=plot(x4a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 4 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x4d
p4d=plot(x4d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 4 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1
plot(p4a, p4d, layout = (1,2))
savefig("signal_haar_level_4.png") 

# x5a
p5a=plot(x5a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 5 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x5d
p5d=plot(x5d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 5 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000 12000;], [0 6000 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1
plot(p5a, p5d, layout = (1,2))
savefig("signal_haar_level_5.png") 


