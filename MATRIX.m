function [epss,FF,BB] = MATRIX(MD,niu,alfa,m,M,l,t,P)

 
    A=zeros(MD,MD);
C=zeros(MD,1);
    %%%%%First ROW%%%%%%%%%%%%%%%%
    A(1,1)=m-l+1;
   
    icol=1;
      for k=1:l
    for n=1:M
      
            icol=icol+1;
            A(1,icol)=0;
            for i=t:m
               
              A(1,icol)=A(1,icol)+(i-k).^(niu+alfa*n-1)  ;
               
            end
            A(1,icol)=A(1,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu);
        end
      end
    
      
      icol=1+l*M;
       for k=1:l
    for n=1:M
      
            icol=icol+1;
            A(1,icol)=0;
            for i=t:m
               
              A(1,icol)=A(1,icol)+(i-k).^(niu+alfa*n-2)  ;
               
            end
            A(1,icol)=A(1,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-1);
        end
       end
      
      
       
       
       
       icol=1+l*M+l*M;
       for k=1:l
    for n=1:M
      
            icol=icol+1;
            A(1,icol)=0;
            for i=t:m
               
              A(1,icol)=A(1,icol)+(i-k).^(niu+alfa*n-3)  ;
               
            end
            A(1,icol)=A(1,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-2);
        end
      end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    
    %%%%%%%%% first scolomn%%%%%%%%%%%%%%%%
       irow=1;
        for k=1:l
    for n=1:M
       
            irow=irow+1;
            A(irow,1)=0;
            for i=t:m
              
              A(irow,1)=A(irow,1)+(i-k).^(niu+alfa*n-1)  ;
               
            end
            A(irow,1)=A(irow,1)*gamma(alfa*n+1)/gamma(alfa*n+niu);
        end
    end
    
       irow=1+l*M;
        for k=1:l
    for n=1:M
       
            irow=irow+1;
            A(irow,1)=0;
            for i=t:m
              
              A(irow,1)=A(irow,1)+(i-k).^(niu+alfa*n-2)  ;
               
            end
            A(irow,1)=A(irow,1)*gamma(alfa*n+1)/gamma(alfa*n+niu-1);
        end
        end
    
     irow=1+l*M+l*M;
        for k=1:l
    for n=1:M
       
            irow=irow+1;
            A(irow,1)=0;
            for i=t:m
              
              A(irow,1)=A(irow,1)+(i-k).^(niu+alfa*n-3)  ;
               
            end
            A(irow,1)=A(irow,1)*gamma(alfa*n+1)/gamma(alfa*n+niu-2);
        end
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
  
    
    %%%%%%%%%%%%%%%%%%%%% 1X1 %%%%%%%%%%
    
        irow=1;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-1)*(i-kpr).^(niu+alfa*npr-1) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu)*gamma(alfa*npr+1)/gamma(alfa*npr+niu);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
    
    
    
     %%%%%%%%%%%%%%%%%%%%%  1X2 %%%%%%%%%%
    
        irow=1;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-2)*(i-kpr).^(niu+alfa*npr-1) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-1)*gamma(alfa*npr+1)/gamma(alfa*npr+niu);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
      %%%%%%%%%%%%%%%%%%%%% 1X3 %%%%%%%%%%
    
        irow=1;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-3)*(i-kpr).^(niu+alfa*npr-1) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-2)*gamma(alfa*npr+1)/gamma(alfa*npr+niu);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%% 2X1 %%%%%%%%%%
    
        irow=1+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-1)*(i-kpr).^(niu+alfa*npr-2) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu)*gamma(alfa*npr+1)/gamma(alfa*npr+niu-1);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
     %%%%%%%%%%%%%%%%%%%%% 2X2 %%%%%%%%%%
    
        irow=1+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-2)*(i-kpr).^(niu+alfa*npr-2) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-1)*gamma(alfa*npr+1)/gamma(alfa*npr+niu-1);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%% 2X3 %%%%%%%%%%
    
        irow=1+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-3)*(i-kpr).^(niu+alfa*npr-2) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-2)*gamma(alfa*npr+1)/gamma(alfa*npr+niu-1);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%% 3X1 %%%%%%%%%%
    
        irow=1+l*M+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-1)*(i-kpr).^(niu+alfa*npr-3) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu)*gamma(alfa*npr+1)/gamma(alfa*npr+niu-2);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%% 3X2 %%%%%%%%%%
    
        irow=1+l*M+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-2)*(i-kpr).^(niu+alfa*npr-3) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-1)*gamma(alfa*npr+1)/gamma(alfa*npr+niu-2);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%% 3X3 %%%%%%%%%%
    
        irow=1+l*M+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+alfa*n-3)*(i-kpr).^(niu+alfa*npr-3) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-2)*gamma(alfa*npr+1)/gamma(alfa*npr+niu-2);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
   %%%%%%%%%%%%%%%%%% C matrica%%%%%%%%%%%%%%%%%%%%%
   C(1)=0;
   for i=t:m
   C(1)=C(1)+P(i);    
   end
   
  
   
   icol=1;
    for k=1:l
   for n=1:M
      
           icol=icol+1;
          C(icol)=0;
   for i=t:m
      
   C(icol)=C(icol)+P(i)*(i-k).^(niu+alfa*n-1);  
     
   end  
       C(icol)=C(icol)*gamma(alfa*n+1)/gamma(alfa*n+niu)  ;
       end
    end
   
   
      icol=1+l*M;
    for k=1:l
   for n=1:M
      
           icol=icol+1;
          C(icol)=0;
   for i=t:m
      
   C(icol)=C(icol)+P(i)*(i-k).^(niu+alfa*n-2);  
     
   end  
        C(icol)=C(icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-1);     
       end
    end
    
    
     icol=1+l*M+l*M;
    for k=1:l
   for n=1:M
      
           icol=icol+1;
          C(icol)=0;
   for i=t:m
      
   C(icol)=C(icol)+P(i)*(i-k).^(niu+alfa*n-3);  
     
   end  
        C(icol)=C(icol)*gamma(alfa*n+1)/gamma(alfa*n+niu-2);     
       end
    end
    
    
   B=A\C; 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 F=zeros(m,1);


for i=t:m
   F(i)=B(1);
  
   irow=1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F(i)=F(i)+gamma(alfa*n+1)/gamma(alfa*n+niu)*(i-k).^(niu+alfa*n-1)*B(irow);
 
       end
      end
   
      
        irow=l*M+1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F(i)=F(i)+gamma(alfa*n+1)/gamma(alfa*n+niu-1)*(i-k).^(niu+alfa*n-2)*B(irow);
 
       end
      end
      
      
           irow=l*M+l*M+1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F(i)=F(i)+gamma(alfa*n+1)/gamma(alfa*n+niu-2)*(i-k).^(niu+alfa*n-3)*B(irow);
 
       end
      end  
      
      
end

dat=0;

eps=0;
for i=1:t-1
   F(i)=P(i);
end
for i=t:m

eps=eps+abs((P(i)-F(i))/P(i));
end
eps=eps*100/(m-t+1);
epss=eps

FF=F;

BB=B;
size(B)
size(BB)
end

