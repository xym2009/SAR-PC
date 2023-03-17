close all;clear;clc;

I=imread('C:\Users\xym\Desktop\XYMRecords\WLH\结果\1m机场\airport.tif');
I=double(I);
tic
% [cv] = cvTH(I,1.5,4);
[pcSum, U, M, Ms, or]=SARPC(I,1.5,4,6,cv);
toc
% or=or/180*pi;
% or(find(or<pi/2))=or(find(or<pi/2))+pi/2;
% or(find(or==pi/2))=or(find(or==pi/2))-pi/2;
% or(find(or>pi/2))=or(find(or>pi/2))-pi/2;
or=or*180/pi;
R1=nonmaxsup(M,or,1.5);
R1=R1./max(R1(:));
U1=hysthresh(R1,0.3,0.15);
figure;imshow(U1);

