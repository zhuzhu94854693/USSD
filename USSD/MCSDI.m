clear
clc
img = imread('C:\Users\aaf\Desktop\USSD\data\Italy_1.bmp');
image = double (img);

image(abs(image)<=0) = min(image(abs(image)>0));
image = log(image+1);
img = (image - min(image(:)))./(max(image(:))-min(image(:)));


[row, col, high] = size(img);
row1=(row/10);row1=round(row1);
col1=(col/10);col1=round(col1);

new_row = ceil(row/row1) * row1;
new_col = ceil(col/col1) * col1;


r_img = img(:, :, 1); 
g_img = img(:, :, 2);
b_img = img(:, :, 3);



new_r_img = imresize(r_img, [new_row new_col]);
new_g_img = imresize(g_img, [new_row new_col]);
new_b_img = imresize(b_img, [new_row new_col]);
new_img(:, :, 1) = new_r_img;
new_img(:, :, 2) = new_g_img;
new_img(:, :, 3) = new_b_img;
[y_row,y_col,dim] = size(new_img);
row_blk_num = y_row/row1;  
col_blk_num = y_col/col1;  
blocks = 1;
for i = 1:row_blk_num
for j = 1:col_blk_num
disp(blocks);
block = new_img((i - 1) * row1 + 1 : i * row1, (j - 1) * col1 + 1 : j * col1, :);
blockall(:,blocks) =block(:);

blocks = blocks + 1;
end
end
X=blockall;
X1=blockall;
[N,D] = size(X);
for i = 1:D
for j = 1:D
local_sim(:,j)=abs(X(:,i)-X1(:,j));

local3 (:,j,i) =local_sim(:,j);
end
end
[x y z]=size(local3);
sum=zeros(x,y);
for i =1:z
    test1 = local3(:,:,i);
    sum = sum+test1;
end
k1 = sum./(z-1);

tt=k1;
for i =1:z
t1 = tt(:,i);
t3 = reshape(t1,[row1,col1,3]);
imwrite(t3, ['C:\Users\aaf\Desktop\USSD\block\' num2str(i) '.jpg']);
end
blocks = 1;
temp = zeros(row1,col_blk_num*col1);
for i = 1:row_blk_num
img = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str((i-1)*col_blk_num+1) '.jpg']);
for j = 1:col_blk_num
img2 = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str(blocks) '.jpg']);
img = [img,img2];
blocks = blocks + 1;
end
if(i==1)
temp = img;
end
if(i~=1)
temp=[temp;img];
end
end

temp=double(temp);
temp = (temp - min(temp(:)))./(max(temp(:))-min(temp(:)));
temp1=temp(:,:,1);
tend=temp1;
[m,n]=size(tend);
[p,q,w]=size(image);
l=(n-q)/2;
tend=tend(:,col1+1:end);
tend=imresize(tend,[p,q]);

imwrite(tend,'C:\Users\aaf\Desktop\USSD\DI\tend1.tif');








tt=blockall;
for i =1:z
t1 = tt(:,i);
t3 = reshape(t1,[row1,col1,3]);
imwrite(t3, ['C:\Users\aaf\Desktop\USSD\block\' num2str(i) '.jpg']);
end
blocks = 1;
temp = zeros(row1,col_blk_num*col1);
for i = 1:row_blk_num
img = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str((i-1)*col_blk_num+1) '.jpg']);
for j = 1:col_blk_num
img2 = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str(blocks) '.jpg']);
img = [img,img2];
blocks = blocks + 1;
end
if(i==1)
temp = img;
end
if(i~=1)
temp=[temp;img];
end
end

temp=double(temp);
temp = (temp - min(temp(:)))./(max(temp(:))-min(temp(:)));
temp1=temp(:,:,1);
tend=temp1;
[m,n]=size(tend);
[p,q,w]=size(image);
l=(n-q)/2;
tend=tend(:,col1+1:end);
tend=imresize(tend,[p,q]);


I=im2double(tend);  
           
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
I=(tend);
Th=k/255;
for i=1:M
for j=1:N
if I(i,j)>=Th
I(i,j)=255;
else
I(i,j)=0;
end
end
end

imwrite(I,'C:\Users\aaf\Desktop\USSD\DI\t3.tif');


clear
clc

img = imread('C:\Users\aaf\Desktop\USSD\data\Italy_2.bmp');
image = double (img);
img = (image - min(image(:)))./(max(image(:))-min(image(:)));




[row, col, high] = size(img);
row1=(row/10);row1=round(row1);
col1=(col/10);col1=round(col1);

new_row = ceil(row/row1) * row1;
new_col = ceil(col/col1) * col1;


r_img = img(:, :, 1); 
g_img = img(:, :, 2);
b_img = img(:, :, 3);




new_r_img = imresize(r_img, [new_row new_col]);
new_g_img = imresize(g_img, [new_row new_col]);
new_b_img = imresize(b_img, [new_row new_col]);
new_img(:, :, 1) = new_r_img;
new_img(:, :, 2) = new_g_img;
new_img(:, :, 3) = new_b_img;
[y_row,y_col,dim] = size(new_img);
row_blk_num = y_row/row1; 
col_blk_num = y_col/col1;  
blocks = 1;
for i = 1:row_blk_num
for j = 1:col_blk_num
disp(blocks);
block = new_img((i - 1) * row1 + 1 : i * row1, (j - 1) * col1 + 1 : j * col1, :);
blockall(:,blocks) =block(:);

blocks = blocks + 1;
end
end
X=blockall;
X1=blockall;
[N,D] = size(X);
for i = 1:D
for j = 1:D
local_sim(:,j)=abs(X(:,i)-X1(:,j));

local4 (:,j,i) =local_sim(:,j);
end
end
[x y z]=size(local4);
sum=zeros(x,y);
for i =1:z
    test1 = local4(:,:,i);
    sum = sum+test1;
end
k2 = sum./(z-1);


tt=k2;
for i =1:z
t1 = tt(:,i);
t3 = reshape(t1,[row1,col1,3]);
imwrite(t3, ['C:\Users\aaf\Desktop\USSD\block\' num2str(i) '.jpg']);
end
blocks = 1;
temp = zeros(row1,col_blk_num*col1);
for i = 1:row_blk_num
img = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str((i-1)*col_blk_num+1) '.jpg']);
for j = 1:col_blk_num
img2 = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str(blocks) '.jpg']);
img = [img,img2];
blocks = blocks + 1;
end
if(i==1)
temp = img;
end
if(i~=1)
temp=[temp;img];
end
end

temp=double(temp);
temp = (temp - min(temp(:)))./(max(temp(:))-min(temp(:)));
temp1=temp(:,:,1);
tend=temp1;
[m,n]=size(tend);
[p,q,w]=size(image);
l=(n-q)/2;
tend=tend(:,col1+1:end);
tend=imresize(tend,[p,q]);

imwrite(tend,'C:\Users\aaf\Desktop\USSD\DI\tend2.tif');




tt=blockall;
for i =1:z
t1 = tt(:,i);
t3 = reshape(t1,[row1,col1,3]);
imwrite(t3, ['C:\Users\aaf\Desktop\USSD\block\' num2str(i) '.jpg']);
end
blocks = 1;
temp = zeros(row1,col_blk_num*col1);
for i = 1:row_blk_num
img = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str((i-1)*col_blk_num+1) '.jpg']);
for j = 1:col_blk_num
img2 = imread(['C:\Users\aaf\Desktop\USSD\block\' num2str(blocks) '.jpg']);
img = [img,img2];
blocks = blocks + 1;
end
if(i==1)
temp = img;
end
if(i~=1)
temp=[temp;img];
end
end

temp=double(temp);
temp = (temp - min(temp(:)))./(max(temp(:))-min(temp(:)));
temp1=temp(:,:,1);
tend=temp1;
[m,n]=size(tend);
[p,q,w]=size(image);
l=(n-q)/2;
tend=tend(:,col1+1:end);
tend=imresize(tend,[p,q]);

I=im2double(tend);  
           
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
I=(tend);
Th=k/255;
for i=1:M
for j=1:N
if I(i,j)>=Th
I(i,j)=255;
else
I(i,j)=0;
end
end
end

imwrite(I,'C:\Users\aaf\Desktop\USSD\DI\t4.tif');



t1=double(imread('C:\Users\aaf\Desktop\USSD\DI\tend1.tif'));
t2=double(imread('C:\Users\aaf\Desktop\USSD\DI\tend2.tif'));
t3=double(imread('C:\Users\aaf\Desktop\USSD\DI\t3.tif'));
t4=double(imread('C:\Users\aaf\Desktop\USSD\DI\t4.tif'));

[M N]=size(t1);
for i=1:M
for j=1:N
if t3(i,j)>=255
t3(i,j)=1;
else
t3(i,j)=0;
end
end
end

[M N]=size(t1);
for i=1:M
for j=1:N
if t4(i,j)>=255
t4(i,j)=1;
else
t4(i,j)=0;
end
end
end

a=t1.*t3; 
b=t2.*t4; 


k1=abs(double(t1)-double(a));
k2=abs(double(t2)-double(b));
f=abs(k2-k1);


f=uint8(f); 
imshow(f);
imwrite(f,'C:\Users\aaf\Desktop\USSD\DI\MSCDI.tif');