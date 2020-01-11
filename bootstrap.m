function se=bootstrap(Y,X,Z,B);

se=1;

N=length(Y);
res=zeros(B,1);
for b=1:B,
    u=floor(N*rand(N,1)+1);
    
    Yb=Y(u,:);
    Xb=X(u,:);
    Zb=Z(u,:);
    ZZ=[Xb,Zb];
    beta=inv(ZZ'*ZZ)*(ZZ'*Yb);
    res(b,1)=beta(1,1);
    end
    
res=sort(res);
std(res);
se=(res(round(B*0.975),1)-res(round(B*0.025),1))/(2*1.96);
    

