clc;
clear all;
P_new =[3 4 3 12 3 11 33 26 43 41 62 39 56 70 117 189 310 228 137 369 613 583 647 938 603 877 1293 1316 1929 2614 2307 2204 2373 2726 3917 3845 4032 3384 5348 3444 3293 4974 3935 7860 4741 4791 3933 4758 4170 4183 5073 5006 5300 4236 3896 4033 4152 4879 4451 4044 3903 3620 3693 5465 5618 4354 3931 3610 3992 5536 5086 4212 3529 3554 3512 3083 2938 3122 3225 3126 3201 2456 2185 1372 2369 2978 2681 2179 1472 1815 1823 1710 1898 1453 1754 1422 1497 1695 1636 1495 1410 1201 1092 1578 909 1147 1396 1291 1372 877 1158 1010 1103 1220 1173 1106 868 834 591 1013 1251 806 816 737 624 751 578 521 626 518 353 583 632 644 513 822 651 530 398 538 641 687 827 726 580 445 560 769 768 767 745 685 581 763 846 880 771 743 928 670 891 950 871 758 1062 816 1148 1009 1129 1440 1077 1040 713 1089 812 1182 1033 1288 1160 853 1184 1048 1522 1276 1108 1715 1406 1295 1508 1735 1940 1813 2988 2948 2460 2659 2919 3539 3497 3330 2619 3105 3991 3395 4322 4516 4156 4437 5079 6899 7357 7518 7579 7713 4662 9116 10157 11047 11754 7070 7981 12594 14542 14162 17540 13863 15166 12872 13970 17234 19722 18965 15650 16170 16981 18803 21330 26684 21243 20530 23012 19790 20890 22884 24701 23065 24405 21915 23254 18950 20018 25174 24141 23287 24957 20572 21350 20412 22950 33470 27301 26860 24962 21363 20051 19609 22915 20252 19875 18662 15450 11299 18213 16272 16022 15871 12155 12329 13430 16169 14879 16298 15538 17272 14717 12282 16578 20964 21671 21502 18447 20263 18450 25161 35383 28507 27052 35928 33363 36804 39237 39036 32725 36532 30501 41385 53135 50023 55892 53284 57725 54990 58784 60916 62322 52618 68053 59937 54940 46169 45533 47525 48682 55757 41346 38598 37535 33355 38905 37892 40261 33552 30004 22194 20089 25308 28680 29079 23275 21088 18607 16840 19202 20634 19114 18262 15845 14104 12364 13013 13494 15143 13308 10972 9765 10624 12718 12057 12027 10406 9834 10641 8489 9938 9985 8523 7434 6035 5455 6391 6385 6572 5946 6040 5177 4712 5766 5926];
P0=P_new;
m=length(P0);
 
P=zeros(m,1);
NPS=1;
optimal_error=1000000;
M=3;
[alfa] = alfavalue(M,P0)
t=5
 
%l=10;
l=t-1;
 
m=length(P0);
%m=59;
P=zeros(m,1);
 
 
x=linspace(1,m,m);

for extr=1:1
%     if (-1)^m==1
%     t=m/2;
%     else
%       t=(m-1)/2;  
%     end
  



if extr==1

 for i=1:m
     P(i)=P0(i);
 end





finalF=zeros(m,1);

else
    for iq=1:m-1
  %P(iq)= x(iq)+2*sin(x(iq)*pi/5);
  
    end
    P(m)=F11(m);
   %P(m)= 32903.53;
 % x=linspace(1,NOD,num1);
    %P(iq)= Psm(iq); 
     %P(iq)=x(iq)+2*sin(x(iq)*pi/5);
finalF=zeros(m,1);
end







MD=1+l*M+l*M+l*M;

NONIU=160;
eps=zeros(NONIU,1);
niu1=linspace(0.01,1,NONIU);
F=zeros(NONIU,m);
B=zeros(NONIU,MD);
     

parfor i1=1:NONIU
    
   
 
    
  [eps(i1),F(i1,:),B(i1,:)]= MATRIX(MD,niu1(i1),alfa,m,M,l,t,P)
   





end


min=10000000000;
for i1=1:NONIU
if(eps(i1)<min)
    min=eps(i1);
    niumin=niu1(i1);
    minerror=min;
    finalF=F(i1,:);
    Bfinal=B(i1,:);
end
end

eps
minerror
niumin
%fclose(fileID);
%e%ps1


 %figure(1)
% plot(finalF)
% 
%hold on;
% plot(P)
 %grid on
 %set(gca,'xtick',[1:1:21]);
%format long
%B


%for id=1:NOD
  
 
  
  % m=m+1;
   
  % num1=num;
  % xsm=linspace(1,num,num1);
   x=linspace(1,m,m);
    Psm=zeros(m,1);
   F11=zeros(m,1);
   for i=1:m
    
   
  Psm(i)=P0(i);
   end

   for i=1:t-1
   %F11(i)=Psm(i);
    F11(i)=P(i);
   end

   
   
   for i=t:m
     
   F11(i)=Bfinal(1);
   irow=1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F11(i)=F11(i)+gamma(alfa*n+1)/gamma(alfa*n+niumin)*(i-k).^(niumin+alfa*n-1)*Bfinal(irow);
 
       end
      end
   
      
        irow=l*M+1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
  
    F11(i)=F11(i)+gamma(alfa*n+1)/gamma(alfa*n+niumin-1)*(i-k).^(niumin+alfa*n-2)*Bfinal(irow);
 
       end
      end
   
      
     irow=l*M+l*M+1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
  
    F11(i)=F11(i)+gamma(alfa*n+1)/gamma(alfa*n+niumin-2)*(i-k).^(niumin+alfa*n-3)*Bfinal(irow);
 
       end
      end
      
      
   end

end
 
alfa
%format long;
%for w=59:m
%F11(w)
%end
   total_error = 0;
  for i=(t):length(Psm)
  total_error = total_error + abs((Psm(i)-F11(i))/(Psm(i)));
  end
  
  total_error = (total_error/(length(Psm)-t+1))*100
