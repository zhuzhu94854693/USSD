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
local1=local3;
local2=local4;
local1(isnan(local1)) = 0;
local2(isnan(local2)) = 0;
[x y z]=size(local1);
sum1=zeros(x,y);
sum2=zeros(x,y);
sum3=zeros(x,y);
sum4=zeros(x,y);
sum5=zeros(x,y);
sum6=zeros(x,y);
sum7=zeros(x,y);
for i =1:z
k1(:,:,i)=(local2(:,:,i).*local1(:,:,i));
sum1 = sum1+k1(:,:,i);
sum2=sum2+local2(:,:,i).*local2(:,:,i);
sum3=sum3+local1(:,:,i).*local1(:,:,i);
end
k1=acos(sum1./(sqrt(sum2).*sqrt(sum3)));
for i=1:z
sum4 = sum4+(local2(:,:,i)-local1(:,:,i)).*(local2(:,:,i)-local1(:,:,i));

end
k2= sqrt(sum4);
k2 = (k2 - min(k2(:)))./(max(k2(:))-min(k2(:)));

[x y z]=size(local1);
mean1=zeros(x,y);
mean2=zeros(x,y);
for i=1:x
for j=1:y
mean1(i,j)=mean(local1(i,j,:));
end
end
for i=1:x
for j=1:y
mean2(i,j)=mean(local2(i,j,:));
end
end
for i =1:z
sum5=sum5+(local2(:,:,i)-mean2).*(local1(:,:,i)-mean1);
sum6=sum6+(local2(:,:,i)-mean2).*(local2(:,:,i)-mean2);
sum7=sum7+(local1(:,:,i)-mean1).*(local1(:,:,i)-mean1);
end
k3=acos(sum5./(sqrt(sum6).*sqrt(sum7)));
k4=k3;
a=k4;
amax = max(max(a));   
amin = min(min(a));  
a=255*(a-amin)/(amax-amin);a=uint8(a);
k4=a;

tt=k4;
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
imshow(tend);
imwrite((tend),'C:\Users\aaf\Desktop\USSD\DI\SCSDI.tif');