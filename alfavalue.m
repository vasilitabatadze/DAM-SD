function [alfaalfa] = alfavalue(M,P0)



J=length(P0);

x=linspace(1,J,J);




A=zeros(M+1,M+1);
B=zeros(M+1,1);

A(1,1)=J;
nal=100;
alfa1=linspace(0.01,1,nal);
min=1000000000000;
for al=1:nal
    
    alfa=alfa1(al);
    
for i=1:M
  sum=0;
  for k=1:J
     sum=sum+CC(x(k),i,alfa) ;
  end
  A(1,i+1)=sum;
end

for i=1:M
  sum=0;
  for k=1:J
     sum=sum+CC(x(k),i,alfa) ;
  end
  A(i+1,1)=sum;
end


   sum=0;
  for k=1:J
     sum=sum+P0(k);
  end
    
  B(1)=sum;
for i=1:M
    for j=1:M
 sum=0;
  for k=1:J
     sum=sum+CC(x(k),j,alfa)*CC(x(k),i,alfa) ;
  end
  A(i+1,j+1)=sum;
        
    end
    
    sum=0;
  for k=1:J
     sum=sum+P0(k)*CC(x(k),i,alfa);
  end
    
  B(i+1)=sum;
  
end

X=A\B;


F=zeros(J,1);
eps1=zeros(J,1);
for i=1:J
   F(i)=X(1);
  
  for n=1:M
    F(i)=F(i)+CC(x(i),n,alfa)*X(n+1);
  end
end

dat=0;

eps=0;

for i=1:J

eps=eps+abs((P0(i)-F(i))/abs(P0(i)+0.00001));
end
eps=eps*100/(J);

if(eps<min)
    min=eps;
    alfamin=alfa;
    minerror=min;
    finalF=F;
    Xfinal=X;
end



end
alfamin
for i=1:J
   F1(i)=Xfinal(1);
  
  for n=1:M
    F1(i)=F1(i)+CC(x(i),n,alfamin)*Xfinal(n+1);
  end
end

alfaalfa=alfamin;





% figure(11)
% plot(P0)
% figure(11)
% hold on;
% plot(F1)
end

