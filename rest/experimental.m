%Experiments with data
close all; clear all; clc;
nx=5;
ny=10;
cnt=0;
C=cell(nx,ny);
for i=1:nx
    for j=1:ny
        C{i,j}=[cnt cnt*2]
        cnt=cnt+1;
    end
end
