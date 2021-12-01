%% Watershed with region merging by Denis Tsygankov

function G=WS_segmentation(F,Fcr)

    [SF, SI]=sort(F(:));     % arranging input pixels by the F-value
    z=find(SF>0,1,'first');  % index of the first non-zero value 
    MY=size(F,1);
    MX=size(F,2);
    L=MX*MY;
    G=zeros(size(F));    % output image
    PV=zeros(L,1);       % peak value for the initial pixel of each region

    cnt=1;               % to label regions as they appear
    
    x=ceil(SI(L)/MY);    % x-coordinate of the starting pixel (the pixel with the largest value)
    y=SI(L)-MY*(x-1);    % y-coordinate of the starting pixel (the pixel with the largest value)
    G(y,x)=cnt;
    PV(cnt)=F(y,x);
    
    % pixel by pixel addition to the regions 
    for i=(L-1):(-1):z
        x=ceil(SI(i)/MY);    % x-coordinate of the new pixel
        y=SI(i)-MY*(x-1);    % y-coordinate of the new pixel
    
        A=G((y-1):(y+1),(x-1):(x+1));   % 3x3 region around the new pixel
        B=unique(A(:));                 % unique values in the 3x3 region  
        LB=length(B);                   % the number of unique values in the 3x3 region
     
        if LB==1                        % if the new pixel is isolated from all regions, start the new region
            cnt=cnt+1;
            G(y,x)=cnt;
            PV(cnt)=F(y,x);                
        elseif LB==2                    % if the new pixel is next to one of the existing regions, add the pixel to this region
            G(y,x)=B(2);                
        elseif LB>=3                    % if the new pixel connects two or more existing regions, check the merging conditions 
            bV=F(y,x);
            tV=PV(B(2:LB));
            [Hv,Hi]=sort(tV-bV);        % sorted difference between the value of the new pixel and the peak values of each region
            if Hv(LB-2)<Fcr             % if all-but-one regions have the peak difference less than Fcr, merge these regions with the region that has the largest peak difference 
                for k=1:(LB-2)
                    G(G==B(Hi(k)+1))=B(Hi(LB-1)+1);                    
                end
                G(y,x)=B(Hi(LB-1)+1);  
            end
        end
    end

