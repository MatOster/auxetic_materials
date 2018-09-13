%script to run big sized connfiguration

%starting configuration
xstart=[.82852399349212646, 0.68926197290420532, 0.0, 0.71778500080108643, 0.32852301001548767, 0.82852399349212646, 0.46778500080108643, 0.65704697370529175, 0.82852399349212646, 1.0177899599075317, 0.65704697370529175, 0.82852399349212646, 0.13926200568675995, 0.0, -0.046308800578117371, 0.57852399349212646, 0.32852301001548767, 0.5, 0.90704697370529175, 1.096310019493103, -0.15704700350761414, 0.90704697370529175, 0.21778500080108643, 0.5, 0.90704697370529175, 0.54630899429321289, 0.5, 0.25, 0.87483197450637817, 0.5, 0.25, 0.11073800176382065, 0.28221499919891357, -0.078523501753807068, 0.32852301001548767, 0.17147700488567352, 0.032214701175689697, 0.65704697370529175, 0.17147700488567352, 1.1248300075531006, 0.0, 0.060738299041986465, -0.078523501753807068, 0.98557001352310181, 0.17147700488567352, 0.36073800921440125, 0.0, 1.0463099479675293, -0.078523501753807068, 0.32852301001548767, 0.17147700488567352, 0.79630899429321289, 0.65704697370529175, 0.61073797941207886, -0.078523501753807068, 0.98557001352310181, -0.15704700350761414, 0.25, 0.43926200270652771, -0.15704700350761414, 0.25, 1.2033599615097046, -0.15704700350761414, 0.90704697370529175, 0.76778501272201538, 1.3140934556722641, 1.3140900135040283]


m=[1:380]       %periodicity in y-direciton
ax = [1:380];   %strain in x-direction
bx=0*[100];     %strain in x-direciton for plot
A=[];
b=[];

%calculating new starting configuration slightly perturbed from old one
deltastart = 1.314593455672264+0.0005;
energ = @(x) energyBig(x,deltastart);
con  = @(x) consBig(x,deltastart);
xa=fmincon(energ,xstart,A,b,A,b,-Inf,Inf,con);
xb=xa;

%output Poisson ratio
outPoissonY=0*[1:100];
outPoissonZ=0*[1:100];
outPoissonInstY=0*[1:100];
outPoissonInstZ=0*[1:100];
outMultipleY=0*[1:388];
outMultipleZ=0*[1:388];


total=1;        %number of runs with same start to be averaged over
for j=1:total
    x0 = xb;%[ 0.75      ,  0.73000002,  0.  ,        0.51999998,  0.25,        0.75,  0.27000001,  0.5 ,        0.75,        0.98000002,  0.5 ,        0.75,        0.02,  0.        ,  0.23,        0.5 ,        0.25      ,  0.5 ,        0.75,  0.98000002,  0.  ,        0.75,        0.02      ,  0.5 ,        0.75,  0.27000001,  0.5 ,        0.25,        0.51999998,  0.5 ,        0.25,        0.23,  0.47999999,  0.  ,        0.25,        0.25      ,  0.23,        0.5 ,        0.25,  0.76999998,  0.  ,        0.02,        0.        ,  0.75,        0.25,  0.47999999,  0.  ,        0.76999998,  0.        ,  0.25,        0.25,  0.51999998,  0.5 ,        0.73000002,  0.        ,  0.75,        0.   ,       0.25,  0.47999999,  0.  ,        0.25,        0.76999998, 0.   ,       0.75,  0.73000002,  1.  ,        1.        ];
    xa=xb;
for i = 1:100
    delta = deltastart+0.5*i*10^(-2)
    ax(i)=0.5*i*10^(-2)/deltastart;
    bx(i)=ax(i);
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBig(x,delta);
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
    start=delta;
for i = 1:70
    delta = start-0.5*i*10^(-2)
    ax(100+i)=-0.5*(i)*10^(-2)/start;
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBig(x,delta);
    rd = randn(1,67);
    x0=(1-0*10^(-3))*x0+.5*10^(-3)*rd;
    %options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(100+i)=x0(66);
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
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBig(x,delta);
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
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBig(x,delta);
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
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBig(x,delta);
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

length(m)

in = [1:380];
plot(in,m)
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(bx,1/total*outPoissonY,bx,1/total*outPoissonZ)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonBig-A.png')
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(bx,1/total*outPoissonInstY,bx,1/total*outPoissonInstZ)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonBig-B.png')
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(in,1/total*outMultipleY,in,1/total*outMultipleZ)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonBig-C.png')
