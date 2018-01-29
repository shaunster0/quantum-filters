function Quantum_file(filtersfile,imfile,ptsfile,rotate,outfile)

filters = dlmread(filtersfile);
descdims = 16;
filters = filters(1:descdims,:);
regions = 2;
flength = size(filters,2);
fsize = sqrt(flength);
filtersSQ = reshape(filters,[size(filters,1) fsize fsize]);
%  filtersSQ = shiftdim(filtersSQ,1);
hsize = (fsize-1)/2;
filters1 = filtersSQ(:,1:hsize,:); filters2 = filtersSQ(:,hsize+1:end,:);
filters1 = reshape(filters1,[size(filters,1) size(filters1,2)*size(filters1,3)]);
filters2 = reshape(filters2,[size(filters,1) size(filters2,2)*size(filters2,3)]);
ptsfilerd = char(ptsfile);
pts = dlmread(ptsfilerd(2:end));
imfilerd = char(imfile);
im = im2single(imread(imfilerd(2:end)));
if (nargin < 5)
    if (~rotate)
        outfile = strcat(ptsfilerd(2:end),'Quantum_up.txt');
    elseif rotate
        outfile = strcat(ptsfilerd(2:end),'Quantum_rot.txt');
    end
end
flength = size(filters,2);
fsize = sqrt(flength);
        
    for k = 2:(size(pts,1))          
        centx = pts(k,1);
        centy = pts(k,2);
        psize = pts(k,3);
        if ~(rotate)
            y1 = round(centy-psize/2); y2 = round(centy+psize/2);
            x1 = round(centx-psize/2); x2 = round(centx+psize/2);
            impatch = im(y1:y2,x1:x2);
        else
            diagsize = psize*1.41;
            y1 = round(centy-diagsize/2); y2 = round(centy+diagsize/2);
            x1 = round(centx-diagsize/2); x2 = round(centx+diagsize/2);
            initpatch = im(y1:y2,x1:x2);
            rotpatch = imrotate(initpatch,-pts(k,4),'bilinear');
            starty = round((size(rotpatch,1) - psize)/2); startx = round((size(rotpatch,2) - psize)/2);
            endy = round((size(rotpatch,1) + psize)/2); endx = round((size(rotpatch,2) + psize)/2);
            impatch = rotpatch(starty:endy,startx:endx);
        end
        impatchR = imresize(impatch,[41 41]);
   %    impatchR = bsxfun(@minus,impatchR,mean(impatchR(:)));
   %    impatchR = bsxfun(@rdivide,impatchR,std(impatchR(:)));
        if (regions == 1)
            desc(k-1,:) = dot(repmat(impatchR(:)',descdims,1),filters,2);
        elseif (regions == 2)
            impatchR1 = impatchR(1:hsize,:); impatchR2 = impatchR(hsize+1:end,:);
            desc(k-1,1:descdims) = dot(repmat(impatchR1(:)',descdims,1),filters1,2);
            desc(k-1,descdims+1:2*descdims) = dot(repmat(impatchR2(:)',descdims,1),filters2,2);
        end
   %     desc(k,:) = desc(k,:)/norm(desc(k,:));   %    
    end

    dlmwrite(outfile, desc, ' ');
    
end

