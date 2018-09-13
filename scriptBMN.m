xstart=[ 0.75, 0.65, 0.0, 0.6, 0.25,    0.75, 0.35, 0.5,    0.75, 0.9, 0.5,    0.75, 0.1, 0.0,    0.15, 0.5, 0.25,    0.5, 0.75, 0.9,    0.0, 0.75, 0.1,    0.5, 0.75, 0.35,    0.5, 0.25, 0.6,    0.5, 0.25, 0.15,    0.4, 0.0, 0.25,    0.25, 0.15, 0.5,    0.25, 0.85, 0.0,    0.1, 0.0, 0.75,    0.25, 0.4, 0.0,    0.85, 0.0, 0.25,    0.25, 0.6, 0.5,    0.65, 0.0, 0.75,    0.0, 0.25, 0.4,    0.0, 0.25, 0.85,    0.0, 0.75, 0.65,    1,1]

m=0*[1:380]
ax = 0*[1:380];
bx=0*[100];
A=[];
b=[];
deltastart = 1+0.002;
energ = @(x) energy(x,deltastart);
con  = @(x) cons1(x,deltastart);
xb=fmincon(energ,xstart,A,b,A,b,-Inf,Inf,con);

outPoissonY=0*[1:100];
outPoissonZ=0*[1:100];
outPoissonInstY=0*[1:100];
outPoissonInstZ=0*[1:100];
outMultipleY=0*[1:388];
outMultipleZ=0*[1:388];


total=1;
for j=1:total
    x0 = xb;%[ 0.75      ,  0.73000002,  0.  ,        0.51999998,  0.25,        0.75,  0.27000001,  0.5 ,        0.75,        0.98000002,  0.5 ,        0.75,        0.02,  0.        ,  0.23,        0.5 ,        0.25      ,  0.5 ,        0.75,  0.98000002,  0.  ,        0.75,        0.02      ,  0.5 ,        0.75,  0.27000001,  0.5 ,        0.25,        0.51999998,  0.5 ,        0.25,        0.23,  0.47999999,  0.  ,        0.25,        0.25      ,  0.23,        0.5 ,        0.25,  0.76999998,  0.  ,        0.02,        0.        ,  0.75,        0.25,  0.47999999,  0.  ,        0.76999998,  0.        ,  0.25,        0.25,  0.51999998,  0.5 ,        0.73000002,  0.        ,  0.75,        0.   ,       0.25,  0.47999999,  0.  ,        0.25,        0.76999998, 0.   ,       0.75,  0.73000002,  1.  ,        1.        ];
    xa=xb;
for i = 1:100
    delta = 1+0.5*i*10^(-2)
    ax(i)=0.5*i*10^(-2)/deltastart;
    bx(i)=ax(i);
    energ = @(x) energy(x,delta);
    con  = @(x) cons1(x,delta);
    rd = randn(1,67);
    dummy1=x0(66);
    dummy2=x0(67);
    x0=(1+0*10^(-3))*x0+.5*10^(-3)*rd;
    %options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(i)=x0(66);
    
    outMultipleY(i)=outMultipleY(i)-(log((x0(66)-xa(66))/xa(66)+1))/log(ax(i)+1);
    outMultipleZ(i) = outMultipleZ(i)-(log((x0(67)-xa(67))/xa(67)+1))/log(ax(i)+1);
    outPoissonY(i)=outPoissonY(i)-(log((x0(66)-xa(66))/xa(66)+1))/log(ax(i)+1);
    outPoissonZ(i) = outPoissonZ(i)-(log((x0(67)-xa(67))/xa(67)+1))/log(ax(i)+1);
    outPoissonInstY(i)=outPoissonInstY(i)+-(log((x0(66)-dummy1)/dummy1+1))/log(ax(i)+1);
    outPoissonInstZ(i) = outPoissonInstZ(i)-(log((x0(67)-dummy2)/dummy2+1))/log(ax(i)+1);
end
    temp1PoissonY = 0*[1:100];
    temp2PoissonZ = 0*[1:100];
    xa=x0;
for i = 1:70
    delta = 1.5-0.5*i*10^(-2)
    ax(100+i)=-0.5*(i)*10^(-2)/1.5;
    energ = @(x) energy(x,delta);
    con  = @(x) cons1(x,delta);
    rd = randn(1,67);
    x0=(1-0*10^(-3))*x0+.5*10^(-3)*rd;
    %options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(100+i)=x0(66);
    [Ab,Bb]=con(x0);
    Ab
    Bb
    temp1PoissonY(i)=-(log((x0(66)-xa(66))/(xa(66))+1))/log(ax(i+100)+1);
    outMultipleY(100+i)=outMultipleY(100+i)+temp1PoissonY(i);
    temp2PoissonZ(i) =  -(log((x0(67)-xa(67))/(xa(67))+1))/log(ax(i+100)+1);
    outMultipleZ(100+i) = outMultipleZ(100+i)+temp2PoissonZ(i);
end
    start=delta;
    temp1PoissonY = 0*[1:100];
    temp2PoissonZ = 0*[1:100];
    xa=x0;
for i = 1:70
    delta = start+0.5*i*10^(-2)
    ax(100+70+i)=(0.5*i*10^(-2)/start);
    energ = @(x) energy(x,delta);
    con  = @(x) cons1(x,delta);
    rd = randn(1,67);
    x0=(1+0*10^(-3))*x0+.5*10^(-3)*rd;
    %options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(100+70+i)=x0(66);
    temp1PoissonY(i)=-(log((x0(66)-xa(66))/xa(66)+1))/log(ax(i+100+70)+1);
    outMultipleY(100+70+i)=outMultipleY(100+70+i)+temp1PoissonY(i);
    temp2PoissonZ(i) =  -(log((x0(67)-xa(67))/xa(67)+1))/log(ax(i+100+70)+1);
    outMultipleZ(100+70+i) = outMultipleZ(100+70+i)+temp2PoissonZ(i);
end
    start=delta
    temp1PoissonY = 0*[1:100];
    temp2PoissonZ = 0*[1:100];
    xa=x0;
for i = 1:70
    delta = start-0.5*i*10^(-2)
    ax(100+2*70+i)=-0.5*(i)*10^(-2)/start;
    energ = @(x) energy(x,delta);
    con  = @(x) cons1(x,delta);
    rd = randn(1,67);
    x0=(1-0*10^(-3))*x0+.5*10^(-3)*rd;
    %options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(100+2*70+i)=x0(66);
    temp1PoissonY(i)=-(log((x0(66)-xa(66))/(xa(66))+1))/log(ax(i+100+2*70)+1);
    outMultipleY(100+2*70+i)=outMultipleY(100+2*70+i)+temp1PoissonY(i);
    temp2PoissonZ(i) =  -(log((x0(67)-xa(67))/(xa(67))+1))/log(ax(i+100+2*70)+1);
    outMultipleZ(100+2*70+i) = outMultipleZ(100+2*70+i)+temp2PoissonZ(i);
end
    start=delta;
    temp1PoissonY = 0*[1:100];
    temp2PoissonZ = 0*[1:100];
    xa=x0;
for i = 1:70
    delta = start+0.5*i*10^(-2)
    ax(100+3*70+i)=(0.5*i*10^(-2)/start);
    energ = @(x) energy(x,delta);
    con  = @(x) cons1(x,delta);
    rd = randn(1,67);
    x0=(1+0*10^(-3))*x0+.5*10^(-3)*rd;
    %options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(100+3*70+i)=x0(66);
    temp1PoissonY(i)=-(log((x0(66)-xa(66))/xa(66)+1))/log(ax(i+100+3*70)+1);
    outMultipleY(100+3*70+i)=outMultipleY(100+3*70+i)+temp1PoissonY(i);
    temp2PoissonZ(i) =  -(log((x0(67)-xa(67))/xa(67)+1))/log(ax(i+100+3*70)+1);
    outMultipleZ(100+3*70+i) = outMultipleZ(100+3*70+i)+temp2PoissonZ(i);
end
end


in = [1:380];
plot(in,m)
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(bx,1/total*outPoissonY,bx,1/total*outPoissonZ)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonMedium-A.png')
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(bx,1/total*outPoissonInstY,bx,1/total*outPoissonInstZ)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonMedium-B.png')
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(in,1/total*outMultipleY,in,1/total*outMultipleZ)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonMedium-C.png')
