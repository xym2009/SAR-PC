
function [ Res1] = cvTH(  M, sigma ,scale )

[rows,cols]=size(M);

Res1=cell(1,scale);
Res2=cell(1,scale);
l=[1,2,3,4];
for n=1:scale
    ax=sigma*1.4^(n-1);
    av=1.8/ax;
    ay=ax;
    h = fspecial('gaussian',l(n)*2+1);
    M2 = imfilter(M,h,'replicate');
    M2=M;
%     w=max(height,width);
    w=l(n);
    d=2*w+1;
    M1=[M2(w+1:-1:2,:);M2;M2(rows-1:-1:rows-w,:)];
    M1=[M1(:,w+1:-1:2),M1,M1(:,cols-1:-1:cols-w)];
    [r,c]=size(M1);
    cv=zeros(r,c);
    for i=(d+1)/2:r-(d-1)/2
        for j=(d+1)/2:c-(d-1)/2
            x1=max(i-(d-1)/2,1);
            x2=min(i+(d-1)/2,r);
            y1=max(j-(d-1)/2,1);
            y2=min(j+(d-1)/2,c);
            R=M1(x1:x2,y1:y2);
            cv(i,j) = sqrt(var(R(:)))/mean(R(:));
        end
    end
%     cv=cv./max(cv(:));

    Res1{n} = cv(w+1:rows+w,w+1:cols+w);
end
end