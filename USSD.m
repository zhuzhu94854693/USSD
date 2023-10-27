clear
clc
e=imread('C:\Users\aaf\Desktop\USSD\DI\MSCDI.tif');
e2=imread('C:\Users\aaf\Desktop\USSD\DI\SCSDI.tif');
e=double(e);
e2=double(e2);

img=e;
img=mat2gray(img);
[m n]=size(img);
r=10;        
imgn=zeros(m+2*r+1,n+2*r+1);
imgn(r+1:m+r,r+1:n+r)=img;
imgn(1:r,r+1:n+r)=img(1:r,1:n);                 
imgn(1:m+r,n+r+1:n+2*r+1)=imgn(1:m+r,n:n+r);  
imgn(m+r+1:m+2*r+1,r+1:n+2*r+1)=imgn(m:m+r,r+1:n+2*r+1);   
imgn(1:m+2*r+1,1:r)=imgn(1:m+2*r+1,r+1:2*r);     


sigma_d=2;
sigma_r=0.1;
[x,y] = meshgrid(-r:r,-r:r);
w1=exp(-(x.^2+y.^2)/(2*sigma_d^2)); 

h=waitbar(0,'wait...');
for i=r+1:m+r
    for j=r+1:n+r        
        w2=exp(-(imgn(i-r:i+r,j-r:j+r)-imgn(i,j)).^2/(2*sigma_r^2));
        w=w1.*w2;
        
        g=imgn(i-r:i+r,j-r:j+r).*w;
        imgn(i,j)=sum(sum(g))/sum(sum(w));
    
    end
    waitbar(i/m);
end
close(h)

e=mat2gray(imgn(r+1:m+r,r+1:n+r));




img=e2;
img=mat2gray(img);
[m n]=size(img);
r=10;        
imgn=zeros(m+2*r+1,n+2*r+1);
imgn(r+1:m+r,r+1:n+r)=img;
imgn(1:r,r+1:n+r)=img(1:r,1:n);                 
imgn(1:m+r,n+r+1:n+2*r+1)=imgn(1:m+r,n:n+r);   
imgn(m+r+1:m+2*r+1,r+1:n+2*r+1)=imgn(m:m+r,r+1:n+2*r+1);    
imgn(1:m+2*r+1,1:r)=imgn(1:m+2*r+1,r+1:2*r);       


sigma_d=2;
sigma_r=0.1;
[x,y] = meshgrid(-r:r,-r:r);
w1=exp(-(x.^2+y.^2)/(2*sigma_d^2));   

h=waitbar(0,'wait...');
for i=r+1:m+r
    for j=r+1:n+r        
        w2=exp(-(imgn(i-r:i+r,j-r:j+r)-imgn(i,j)).^2/(2*sigma_r^2)); 
        w=w1.*w2;
        
        g=imgn(i-r:i+r,j-r:j+r).*w;
        imgn(i,j)=sum(sum(g))/sum(sum(w));
    
    end
    waitbar(i/m);
end
close(h)

e2=mat2gray(imgn(r+1:m+r,r+1:n+r));



k=double(e2).*double(e);
k = 255.*(k - min(k(:)))./(max(k(:))-min(k(:)));
imwrite(uint8(k),'C:\Users\aaf\Desktop\USSD\DI\combine.tif');

I=im2double(imread('C:\Users\aaf\Desktop\USSD\DI\combine.tif'));

[M,N]=size(I);                 
number_all=M*N;             
hui_all=0;                    
ICV_t=0;                         


for i=1:M
    for j=1:N
        hui_all=hui_all+I(i,j);
    end
end
all_ave=hui_all*255/number_all;   



for t=0:255                       

    hui_A=0;                      
    hui_B=0;        

    number_A=0;                  
    number_B=0;                 
 
    for i=1:M                   
        for j=1:N
            if (I(i,j)*255>=t)    
                number_A=number_A+1;  
                hui_A=hui_A+I(i,j);   
            elseif (I(i,j)*255<t) 
                number_B=number_B+1;  
                hui_B=hui_B+I(i,j);   
                
            end
        end
    end

    PA=number_A/number_all;           
    PB=number_B/number_all;          
  
    A_ave=hui_A*255/number_A;        
    B_ave=hui_B*255/number_B;         

    ICV=PA*((A_ave-all_ave)^2)+PB*((B_ave-all_ave)^2); 
   if (ICV>ICV_t)                     
        ICV_t=ICV;
        k=t;                          
    end

end

k                                     
I=imread('C:\Users\aaf\Desktop\USSD\DI\combine.tif');
Th=k;
for i=1:M
    for j=1:N
        if I(i,j)>=Th
            I(i,j)=255;
        else
            I(i,j)=0;
        end
    end
end
figure,imshow(I);
imwrite(I,'C:\Users\aaf\Desktop\USSD\result\USSD.tif');
