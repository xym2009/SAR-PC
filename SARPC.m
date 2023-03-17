function [ pcSum, U, Mome, Ms, or ] = SARPC( M, sigma ,scale, orient, cv1)
%---------%

epsilon   = .0001;            % Used to prevent division by zero.
M = M+epsilon;
[rows,cols] = size(M);
U=cell(1,orient);
X=cell(1,orient);
or = zeros(rows,cols); 
pcSum2= zeros(rows,cols); 
UR = cell(scale, orient);
UI = cell(scale, orient);
covx2 = zeros(rows,cols);                     % Matrices for covariance data
covy2 = zeros(rows,cols);
covxy = zeros(rows,cols);
sumAn_all = zeros(rows,cols);
EnergyV = zeros(rows,cols,3);  

pcSum=0;
for o=1:orient
    sumE_ThisOrient   = zeros(rows,cols);          % Initialize accumulator matrices.
    sumO_ThisOrient   = zeros(rows,cols);
    sumAn_ThisOrient  = zeros(rows,cols);
    Energy            = zeros(rows,cols); 
    theta=(o-1)*pi/orient;
    
    for n=1:scale
        UI{n,o}=zeros(rows,cols);
        UR{n,o}=zeros(rows,cols);
     
        ax=sigma*1.4^(n-1);
        sig=1.8;
        av=sig/ax;
        ay=ax;
        pm=exp(-sig*sig/2);
        width=ceil(1.5*ax)*2+1;
        height=ceil(1.5*ay)*2+1;
        d1=(pi/(2*av));
        w=max(height,width);
        w=(w-1)/2;
        M1=[M(w+1:-1:2,:);M;M(rows-1:-1:rows-w,:)];
        M1=[M1(:,w+1:-1:2),M1,M1(:,cols-1:-1:cols-w)];
        
        upperW = zeros(2*w+1, 2*w+1);
        lowerW = upperW;
        upperD = zeros(2*w+1, 2*w+1);
        ceterD = upperD;
        lowerD = upperD;
        for x = -w:w
            for y = -w:w
                rx = [cos(theta) sin(theta)]*[x; y];
                ry = [-sin(theta) cos(theta)]*[x; y];
                if rx >0
                    upperW(x+w+1, y+w+1) = exp(-0.5*(rx^2/ax^2+ry^2/ay^2))*abs(sin(av*(x*cos(theta)+y*sin(theta))));
                elseif rx<0
                    lowerW(x+w+1, y+w+1) = exp(-0.5*(rx^2/ax^2+ry^2/ay^2))*abs(sin(av*(x*cos(theta)+y*sin(theta))));
                end
                if rx >d1
                    upperD(x+w+1, y+w+1) = exp(-0.5*(rx^2/ax^2+ry^2/ay^2))*abs(cos(av*(x*cos(theta)+y*sin(theta)))-pm);
                elseif rx<-d1
                    lowerD(x+w+1, y+w+1) = exp(-0.5*(rx^2/ax^2+ry^2/ay^2))*abs(cos(av*(x*cos(theta)+y*sin(theta)))-pm);
                elseif rx<d1&&rx>-d1
                    ceterD(x+w+1, y+w+1) = exp(-0.5*(rx^2/ax^2+ry^2/ay^2))*abs(cos(av*(x*cos(theta)+y*sin(theta)))-pm);
                end
            end
        end
        upperW = upperW/sum(sum(upperW));
        lowerW = lowerW/sum(sum(lowerW));
        upperD = upperD/sum(sum(upperD));
        lowerD = lowerD/sum(sum(lowerD));
        ceterD = ceterD/sum(sum(ceterD));
        a=filter2(upperW,M1,'same');
        a=max(a,1);
        b=filter2(lowerW,M1,'same');
        b=max(b,1);
        c=filter2(upperD,M1,'same');
        c=max(c,1);
        d=filter2(lowerD,M1,'same');
        d=max(d,1);
        e=filter2(ceterD,M1,'same');
        e=max(e,1);
        
        R1=1-min(a./b,b./a);
        LI=R1(w+1:rows+w,w+1:cols+w);
        R2=1-min((c)./e,e./(c));
        R3=1-min((d)./e,e./(d));
        R2=min(R2,R3);
        LR=R2(w+1:rows+w,w+1:cols+w);
        mycv = cv1{n};
        UI{n,o}=LI;
        UI{n,o}=UI{n,o}./max(UI{n,o}(:));
        UR{n,o}=LR;
        UR{n,o}=1.0*UR{n,o}./max(UR{n,o}(:));
        An=sqrt(UR{n,o}.*UR{n,o}+UI{n,o}.*UI{n,o});
        sumAn_ThisOrient = sumAn_ThisOrient + An;  % Sum of amplitude responses.
        sumE_ThisOrient = sumE_ThisOrient + UR{n,o}; % Sum of even filter convolution results.
        sumO_ThisOrient = sumO_ThisOrient + UI{n,o}; % Sum of odd filter convolution results.
        EnergyV(:,:,1) = EnergyV(:,:,1) + sumE_ThisOrient;
        EnergyV(:,:,2) = EnergyV(:,:,2) + cos(theta)*sumO_ThisOrient;
        EnergyV(:,:,3) = EnergyV(:,:,3) + sin(theta)*sumO_ThisOrient;
        if n == 1
            maxAn = An;
            maxCv1 = mycv;
        else
            maxAn = max(maxAn,An);
            maxCv1 = min(maxCv1,mycv);
        end
    end
    width = (sumAn_ThisOrient./(maxAn + epsilon) - 1) / (scale-1);    
    weight = 1.0 ./ (1 + exp( (0.5 - width)*10)); 
    XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;  
    MeanE = sumE_ThisOrient ./ XEnergy; 
    MeanO = sumO_ThisOrient ./ XEnergy; 
    for n = 1:scale
        E = UR{n,o}; O = UI{n,o};    % Extract even and odd
        % convolution results.
        Energy = Energy + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
    end
    
    GT1=(0.5+0.25*log(1./(0.0001+maxCv1))); 
    T1=GT1;
%     T1=median(Energy(:));
    Energy = max(Energy-T1 , 0);  
    U{o}=weight.*Energy./sumAn_ThisOrient;
    X{o}=Energy./sumAn_ThisOrient;
    pcSum = pcSum+U{o};   
    pcSum2=pcSum2+X{o};
    sumAn_all = sumAn_all + sumAn_ThisOrient;
    covx = U{o}*cos(theta);
    covy = U{o}*sin(theta);
    covx2 = covx2 + covx.^2;
    covy2 = covy2 + covy.^2;
    covxy = covxy + covx.*covy;
end

covx2 = covx2/(orient/2);
covy2 = covy2/(orient/2);
covxy = 4*covxy/orient;   % This gives us 2*covxy/(norient/2)
denom = sqrt(covxy.^2 + (covx2-covy2).^2)+epsilon;
Mome = (covy2+covx2 + denom)/2;          % Maximum moment
Ms = (covy2+covx2 - denom)/2;          % ... and minimum moment

or = atan2(EnergyV(:,:,3), EnergyV(:,:,2));
or(or<0) = or(or<0)+pi;       % Wrap angles -pi..0 to 0..pi

end
          