%script for commpressing big sized structure

%starting configuration
xstart=[.82852399349212646, 0.68926197290420532, 0.0, 0.71778500080108643, 0.32852301001548767, 0.82852399349212646, 0.46778500080108643, 0.65704697370529175, 0.82852399349212646, 1.0177899599075317, 0.65704697370529175, 0.82852399349212646, 0.13926200568675995, 0.0, -0.046308800578117371, 0.57852399349212646, 0.32852301001548767, 0.5, 0.90704697370529175, 1.096310019493103, -0.15704700350761414, 0.90704697370529175, 0.21778500080108643, 0.5, 0.90704697370529175, 0.54630899429321289, 0.5, 0.25, 0.87483197450637817, 0.5, 0.25, 0.11073800176382065, 0.28221499919891357, -0.078523501753807068, 0.32852301001548767, 0.17147700488567352, 0.032214701175689697, 0.65704697370529175, 0.17147700488567352, 1.1248300075531006, 0.0, 0.060738299041986465, -0.078523501753807068, 0.98557001352310181, 0.17147700488567352, 0.36073800921440125, 0.0, 1.0463099479675293, -0.078523501753807068, 0.32852301001548767, 0.17147700488567352, 0.79630899429321289, 0.65704697370529175, 0.61073797941207886, -0.078523501753807068, 0.98557001352310181, -0.15704700350761414, 0.25, 0.43926200270652771, -0.15704700350761414, 0.25, 1.2033599615097046, -0.15704700350761414, 0.90704697370529175, 0.76778501272201538, 1.3140934556722641, 1.3140900135040283]

out=0*[1:200];  %output strain in y-direction
m=0*[1:200]     %output periodicity in y-direction
ax = 0*[1:200]; %strain in x-direction
A=[];
b=[];

%calculating new starting configuration with same boundary conditions but changing minimal spring length
deltastart = 1.314593455672264;
energ = @(x) energyBig(x,deltastart);
con  = @(x) consBigComp(x,deltastart,.130);
xa=fmincon(energ,xa,A,b,A,b,-Inf,Inf,con);
con=@(x) consBigComp(x,deltastart,.122);
xa=fmincon(energ,xa,A,b,A,b,-Inf,Inf,con);
con=@(x) consBigComp(x,deltastart,.118);
xa=fmincon(energ,xa,A,b,A,b,-Inf,Inf,con);
con=@(x) consBigComp(x,deltastart,.112);
xa=fmincon(energ,xa,A,b,A,b,-Inf,Inf,con);
% con=@(x) consBigComp(x,deltastart,.115);
% xa=fmincon(energ,xa,A,b,A,b,-Inf,Inf,con);
% con=@(x) consBigComp(x,deltastart,.112);
% xa=fmincon(energ,xa,A,b,A,b,-Inf,Inf,con);
xb=xa;
out2=0*[1:200]; %output strain in z-direciton
total=1;
for j=1:total
    temp1 = 0*[1:200];
    temp2 = 0*[1:200];
    x0=xb;%[.82852399349212646, 0.68926197290420532, 0.0, 0.71778500080108643, 0.32852301001548767, 0.82852399349212646, 0.46778500080108643, 0.65704697370529175, 0.82852399349212646, 1.0177899599075317, 0.65704697370529175, 0.82852399349212646, 0.13926200568675995, 0.0, -0.046308800578117371, 0.57852399349212646, 0.32852301001548767, 0.5, 0.90704697370529175, 1.096310019493103, -0.15704700350761414, 0.90704697370529175, 0.21778500080108643, 0.5, 0.90704697370529175, 0.54630899429321289, 0.5, 0.25, 0.87483197450637817, 0.5, 0.25, 0.11073800176382065, 0.28221499919891357, -0.078523501753807068, 0.32852301001548767, 0.17147700488567352, 0.032214701175689697, 0.65704697370529175, 0.17147700488567352, 1.1248300075531006, 0.0, 0.060738299041986465, -0.078523501753807068, 0.98557001352310181, 0.17147700488567352, 0.36073800921440125, 0.0, 1.0463099479675293, -0.078523501753807068, 0.32852301001548767, 0.17147700488567352, 0.79630899429321289, 0.65704697370529175, 0.61073797941207886, -0.078523501753807068, 0.98557001352310181, -0.15704700350761414, 0.25, 0.43926200270652771, -0.15704700350761414, 0.25, 1.2033599615097046, -0.15704700350761414, 0.90704697370529175, 0.76778501272201538, 1.3140934556722641, 1.3140900135040283]

for i = 1:100
    delta = deltastart-0.5*i*10^(-2)
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBigComp(x,delta,.112);
    rd = randn(1,67);
    ax(i)=-0.5*i*10^(-2)/deltastart;    
    dummy1=x0(66);
    dummy2=x0(66);
    x0=(1+0*10^(-3))*x0+.5*10^(-3)*rd;
    options = optimoptions(@fmincon,'GradObj','on', 'GradConstr','on');%'Algorithm', 'interior-point',
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(i)=x0(66);
    %xa(66)=dummy1;
    %xa(67)=dummy2;
    temp1(i)=-(log((x0(66)-xa(66))/xa(66)+1))/log(ax(i)+1);
    out(i)=out(i)+temp1(i);
    temp2(i) =  -(log((x0(67)-xa(67))/xa(67)+1))/log(ax(i)+1);
    out2(i) = out2(i)+temp2(i);
end
    deltastart=delta;
    xa=x0;
for i = 1:100
    delta = deltastart+0.5*i*10^(-2)
    energ = @(x) energyBig(x,delta);
    con  = @(x) consBigComp(x,delta,.112);
    rd = randn(1,67);
    ax(i)=0.5*i*10^(-2)/deltastart;    
    dummy1=x0(66);
    dummy2=x0(66);
    x0=(1+0*10^(-3))*x0+.8*10^(-3)*rd;
    options = optimoptions(@fmincon,'GradObj','on', 'GradConstr','on');%'Algorithm', 'interior-point',
    x0=fmincon(energ,x0,A,b,A,b,-Inf,Inf,con);
    m(100+i)=x0(66);
    %xa(66)=dummy1;
    %xa(67)=dummy2;
    temp1(i)=-(log((x0(66)-xa(66))/xa(66)+1))/log(ax(i)+1);
    out(100+i)=out(100+i)+temp1(i);
    temp2(i) =  -(log((x0(67)-xa(67))/xa(67)+1))/log(ax(i)+1);
    out2(100+i) = out2(100+i)+temp2(i);
end

end
in = [1:200];
plot(in,m)
figure 
title('Poisson ratio in y and z direction')
xlabel('(engineers) strain in x direction')
ylabel('(engineers) Poisson ratio')
plot(in,1/total*out,in,1/total*out2)
legend('Poisson ration in y-direction','Poisson ratio in z-direction')
saveas(gcf,'PoissonBigComp.png')
