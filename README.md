Linye Zhu, Wenbin Sun, Deqin Fan, Huaqiao Xing, Xiaoqi Liu, Unsupervised spatial self-similarity difference-based change detection method for multi-source heterogeneous images, Pattern Recognition, Volume 149, 2024, 110237, ISSN 0031-3203, https://doi.org/10.1016/j.patcog.2023.110237.
(https://www.sciencedirect.com/science/article/pii/S0031320323009342)


# USSD
Figure 1 shows the framework of the USSD multi-source heterogeneous change detection method. 
There are four main steps: (1) Constructing spatial self-difference images: The heterogeneous images are spatially blocked separately, and the spatially blocked images are subjected to self-differential calculations using the absolute distance method, which is used for comparing the differences between image blocks and image blocks. (2) Obtaining the heterogeneous spatial self-difference change magnitude map based on the magnitude difference: The heterogeneous spatial self-difference change magnitude map based on the magnitude difference is obtained by combining the binarised heterogeneous image results to constrain the spatial self-difference image and calculating the average difference between image blocks. (3) Obtaining the heterogeneous spatial self-difference change magnitude map based on the similarity difference: The correlation coefficient method is used to calculate the similarity of the two heterogeneous spatial self-difference images, thereby obtaining the heterogeneous spatial self-difference change magnitude map based on the similarity difference. (4) Acquiring the final heterogeneous spatial self-difference change magnitude map: The two types of heterogeneous spatial self-difference change magnitude maps are enhanced and denoised using bilateral filtering, and then the final heterogeneous spatial self-difference change magnitude map is obtained via dot product.
![Image text](https://github.com/zhuzhu94854693/USSD/blob/main/image.png)
Figure 1. Framework of the USSD method

