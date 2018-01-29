function bestfilters = selectFilters(selectfac, patchnorm, savefile, savefile1, savefile2, parfile)

cd('C:\Dropbox\quantum');
load schrod_filtsBig3.mat  %fourquad.mat  %
descdims = 16; % 25;
flength = size(filters,2);
fsize = sqrt(flength);
filtersSQ = reshape(filters,[size(filters,1) fsize fsize]);
%  filtersSQ = shiftdim(filtersSQ,1);
hsize = (fsize-1)/2;
filters1 = filtersSQ(:,1:hsize,:); filters2 = filtersSQ(:,hsize+1:end,:);
filters1 = reshape(filters1,[size(filters,1) size(filters1,2)*size(filters1,3)]);
filters2 = reshape(filters2,[size(filters,1) size(filters2,2)*size(filters2,3)]);
filtersSQV = repmat(shiftdim(filtersSQ,1), [2 1]);
filtersSQD = shiftdim(imresize(filtersSQV, [fsize fsize]),2);
filters1D = filtersSQD(:,1:hsize,:); filters2D = filtersSQD(:,hsize+1:end,:);
filters1D = reshape(filters1D,[size(filters,1) size(filters1D,2)*size(filters1D,3)]);
filters2D = reshape(filters2D,[size(filters,1) size(filters2D,2)*size(filters2D,3)]);
direct = {'boat' 'wall' 'bark' 'graf' 'trees' 'leuven' 'bikes' 'ubc'};
scale = 1;
surfPts = true;
randPts = false;
rotate = true;
noseqs = 8; %8
noims = 4;  %6
imsch = [1 2 3 4]; invimsch = [0 1 2 3];  %imsch = [1 2 3 4]; invimsch = [0 1 2 3];
ent = @(j)-sum(j(j~=0).*log2(j(j~=0)));

i = 1;
for j = 1:noims
    im{i,j} = im2single(imread(char(strcat('./testing/',direct(i),'/img',num2str(imsch(j)),'.pgm'))));
end
for i = 2:noseqs
    for j = 1:noims 
        im{i,j} = im2single(rgb2gray(imread(char(strcat('./testing/',direct(i),'/img',num2str(imsch(j)),'.ppm')))));
    end
end

for i = 1:noseqs
    for j = 2:noims
        H(:,:,i,j-1) = load(char(strcat('./testing/',direct(i),'/H1to',num2str(imsch(j)),'pVL.txt')));
        the = atan2(H(2,1,i,j-1),H(1,1,i,j-1));
        sc = sqrt(H(1,1,i,j-1)*H(2,2,i,j-1) - H(1,2,i,j-1)*H(2,1,i,j-1));
        offx = H(1,3,i,j-1); offy = H(2,3,i,j-1); 
        hompars(:,i,j-1) = [the; sc; offx; offy];
    end
end

if (surfPts)
  i = 1;
  in = opencv_read(char(strcat('./testing/',direct(i),'/img1.pgm',direct(i),'.xml')))';
  input{i} = in(:,(in(3,:) > 16) & (in(3,:) < 80));
  surfsize(1,:,i) = size(input{i}');
  for i = 2:noseqs
    in = opencv_read(char(strcat('./testing/',direct(i),'/img1.ppm',direct(i),'.xml')))';
    input{i} = in(:,(in(3,:) > 16) & (in(3,:) < 80));
    surfsize(1,:,i) = size(input{i}');
  end
end

%im = imresize(im, scale);
%psize = 40*scale;
bestfilters = [];
matchrate = [];
flength = size(filters,2);
fsize = sqrt(flength);
regions = 2;
initpts = 2100;
notrainpts = 500; 
nofilters = round(size(filters,1)/3);
inord = 1:notrainpts;

w = ceil(765*scale);
h = ceil(512*scale);
iind(1) = 1; 

for i = 1:descdims
    
   randind = randperm(size(filters,1));
   if (regions == 1)
      useFilters = filters(randind(1:nofilters),:); % 1:nofilters
   elseif (regions == 2)
      useFilters1 = filters1(randind(1:nofilters),:); 
      useFilters2 = filters2(randind(1:nofilters),:);
   end
   
   useParams = params(randind(1:nofilters),:);
    
   % set sequence indices
   seqlen = round(initpts/noseqs);
   tryseqind = [];
   for n = 1:noseqs-1
       tryseqind = [tryseqind n*ones(1,seqlen)];
   end
   n = noseqs;
   tryseqind = [tryseqind n*ones(1,initpts-(n-1)*seqlen)];
   
   % set image indices
   tryimind = randi(imsch(noims)-1,1,initpts)+1;
    
   % get random points
   if(randPts)
       xpts = randi(w,1,initpts);
       ypts = randi(h,1,initpts);
       pts = [xpts/scale;ypts/scale;ones(1,initpts); ...
           40*ones(1,initpts);zeros(1,initpts)];
   end
   
   % get random surf points
   if(surfPts)
      inpts = [];
      for k = 1:noseqs-1
          order{k} = randperm(surfsize(1,1,k));
          list = order{k};
          indx = list(1:seqlen);
          ptss = input{k};
          inpts = [inpts ptss(:,indx)];
      end
      k = noseqs;
      order{k} = randperm(surfsize(1,1,k));
      list = order{k};
      indx = list(1:(initpts-(noseqs-1)*seqlen)); %   indx = 1:initpts;
      ptss = input{k};
      inpts = [inpts ptss(:,indx)];
      pts = [inpts(1:2,:); ones(1,initpts); inpts(3:4,:)];     
   end

   for k = 1:initpts
      corpts(1:3,k) = H(:,:,tryseqind(k),invimsch(tryimind(k)))*pts(1:3,k);
      corpts(4,k) = pts(4,k).*abs(hompars(2,tryseqind(k),invimsch(tryimind(k))));
      if (rotate)
         corpts(5,k) = pts(5,k)-hompars(1,tryseqind(k),invimsch(tryimind(k)))/pi()*180;
      end
   end
%   pts = [xpts*scale;ypts*scale;ones(1,initpts)];
   corpts(1,:) = corpts(1,:)./corpts(3,:)*scale;
   corpts(2,:) = corpts(2,:)./corpts(3,:)*scale;
   corpts(3,:) = 1;
   if ~(rotate)
       corpts(5,:) = pts(5,:);
   end
   noiscorpts(1,:) = corpts(1,:) + 1*scale*rand(1,initpts);
   noiscorpts(2,:) = corpts(2,:) + 1*scale*rand(1,initpts);
   noiscorpts(3:5,:) = corpts(3:5,:);
   edg = round(1.41*pts(4,:)/2+1);
   oneok = (pts(1,:) > edg) & (pts(2,:) > edg) & (pts(1,:) < (w-edg)) & (pts(2,:) < (h-edg));
   twook = (noiscorpts(1,:) > edg) & (noiscorpts(2,:) > edg) & (noiscorpts(1,:) < (w-edg)) & (noiscorpts(2,:) < (h-edg));
   bothok = (oneok & twook);
   permind = randperm(sum(bothok));
   okpts =  pts(:,bothok);
   spts = okpts(:,permind(1:notrainpts));
   noisokpts =  noiscorpts(:,bothok);
   spts(:,:,2) = noisokpts(:,permind(1:notrainpts));
   corokpts =  corpts(:,bothok);
   spts(:,:,3) = corokpts(:,permind(1:notrainpts));
   edgokpts =  edg(:,bothok);
   edg = edgokpts(:,permind(1:notrainpts));
   tryseqind = tryseqind(bothok);
   tryseqind = tryseqind(permind(1:notrainpts));
   tryimind = tryimind(bothok);
   tryimind = tryimind(permind(1:notrainpts));
   
   for m = 1:(i-1)
 
       if (regions == 1)
          curfilter = reshape(bestfilters(m,:),fsize,fsize);
          curfilter = imresize(curfilter, scale);
       elseif (regions == 2)
          curfilter1 = reshape(bestfilters1(m,:),hsize,fsize);
          curfilter2 = reshape(bestfilters2(m,:),hsize+1,fsize);
          curfilter1 = imresize(curfilter1, scale); 
          curfilter2 = imresize(curfilter2, scale); 
       end

       for k = 1:notrainpts %size(spts,2)

           iind(2) = invimsch(tryimind(k))+1;
           for imano = 1:2
               centx = spts(1,k,imano);
               centy = spts(2,k,imano);
               psize = spts(4,k,imano);
               if ~(rotate)
                   y1 = round(centy-psize/2); y2 = round(centy+psize/2);
                   x1 = round(centx-psize/2); x2 = round(centx+psize/2);
                   impatch = im{tryseqind(k), iind(imano)}(y1:y2,x1:x2);
               else
                   diagsize = psize*1.41;
                   y1 = round(centy-diagsize/2); y2 = round(centy+diagsize/2);
                   x1 = round(centx-diagsize/2); x2 = round(centx+diagsize/2);
                   initpatch = im{tryseqind(k), iind(imano)}(y1:y2,x1:x2);
                   rotpatch = imrotate(initpatch,-spts(5,k,imano),'bilinear');
                   starty = round((size(rotpatch,1) - psize)/2); startx = round((size(rotpatch,2) - psize)/2);
                   endy = round((size(rotpatch,1) + psize)/2); endx = round((size(rotpatch,2) + psize)/2);
                   impatch = rotpatch(starty:endy,startx:endx);
               end
               impatchR = imresize(impatch,[41 41]);
        %       impatchR = bsxfun(@minus,impatchR,mean(impatchR(:)));
               if patchnorm
                  impatchR = bsxfun(@rdivide,impatchR,std(impatchR(:)));
               end
               if (regions == 1)
                  desc(k,m,imano) = sum(sum(impatchR(:,:).*curfilter));
               elseif (regions == 2)
                  impatchR1 = impatchR(1:hsize,:); impatchR2 = impatchR(hsize+1:end,:);
                  desc(k,2*m-1,imano) = sum(sum(impatchR1.*curfilter1)); 
                  desc(k,2*m,imano) = sum(sum(impatchR2.*curfilter2));             
%                        desc4(k,4*m-3,imano) = sum(sum(impatch(1:psize/2,1:psize/2,imano).*curfilter(1:psize/2,1:psize/2)));
%                        desc4(k,4*m-2,imano) = sum(sum(impatch(1:psize/2,psize/2+1:end,imano).*curfilter(1:psize/2,psize/2+1:end)));
%                        desc4(k,4*m-1,imano) = sum(sum(impatch(psize/2+1:end,1:psize/2,imano).*curfilter(psize/2+1:end,1:psize/2)));
%                        desc4(k,4*m,imano) = sum(sum(impatch(psize/2+1:end,psize/2+1:end,imano).*curfilter(psize/2+1:end,psize/2+1:end)));
%                     end
               end
           end

       end
 
   end

       for k = 1:notrainpts %size(spts,2)
           
           iind(2) = invimsch(tryimind(k))+1;
           for imano = 1:2
               centx = spts(1,k,imano);
               centy = spts(2,k,imano);
               psize = spts(4,k,imano);
               if ~(rotate)
                   y1 = round(centy-psize/2); y2 = round(centy+psize/2);
                   x1 = round(centx-psize/2); x2 = round(centx+psize/2);
                   impatch = im{tryseqind(k), iind(imano)}(y1:y2,x1:x2);
               else
                   diagsize = psize*1.41;
                   y1 = round(centy-diagsize/2); y2 = round(centy+diagsize/2);
                   x1 = round(centx-diagsize/2); x2 = round(centx+diagsize/2);
                   initpatch = im{tryseqind(k), iind(imano)}(y1:y2,x1:x2);
                   rotpatch = imrotate(initpatch,-spts(5,k,imano),'bilinear');
                   starty = round((size(rotpatch,1) - psize)/2); startx = round((size(rotpatch,2) - psize)/2);
                   endy = round((size(rotpatch,1) + psize)/2); endx = round((size(rotpatch,2) + psize)/2);
                   impatch = rotpatch(starty:endy,startx:endx);
               end
               impatchR = imresize(impatch,[41 41]);
      %         impatchR = bsxfun(@minus,impatchR,mean(impatchR(:)));
               if patchnorm
                  impatchR = bsxfun(@rdivide,impatchR,std(impatchR(:)));
               end
               if (regions == 1)
                  convols(:,k,imano) = dot(repmat(impatchR(:)',nofilters,1),useFilters,2);
               elseif (regions == 2)
                  curimpat1 = impatchR(1:hsize,:); curimpat2 = impatchR(hsize+1:end,:);
                  convols1(:,k,imano) = dot(repmat(curimpat1(:)',nofilters,1),useFilters1,2);
                  convols2(:,k,imano) = dot(repmat(curimpat2(:)',nofilters,1),useFilters2,2);
%                        desc4(k,4*i-3,imano) = sum(sum(impatch(1:psize/2,1:psize/2).*curfilter(1:psize/2,1:psize/2)));
%                        desc4(k,4*i-2,imano) = sum(sum(impatch(1:psize/2,psize/2+1:end).*curfilter(1:psize/2,psize/2+1:end)));
%                        desc4(k,4*i-1,imano) = sum(sum(impatch(psize/2+1:end,1:psize/2).*curfilter(psize/2+1:end,1:psize/2)));
%                        desc4(k,4*i,imano) = sum(sum(impatch(psize/2+1:end,psize/2+1:end).*curfilter(psize/2+1:end,psize/2+1:end)));
%                     end
               end
           end
           
       end
       
       for p = 1:nofilters  % /5
            
          if (regions == 1)
              if (i == 1) 
                  tent1 = [convols(p,:,1)];
                  tent2 = [convols(p,:,2)];
              else
                  tent1  = [desc(:,:,1)'; convols(p,:,1)];
            %      tent1 = tent1./repmat(sqrt(sum(tent1(:,:).^2,1)),[i 1]);   %    
                  tent2  = [desc(:,:,2)'; convols(p,:,2)];
            %      tent2 = tent2./repmat(sqrt(sum(tent2(:,:).^2,1)),[i 1]);   % 
              end
          elseif (regions == 2)
              if (i == 1) 
                  tent1 = [convols1(p,:,1); convols2(p,:,1)];
                  tent2 = [convols1(p,:,2); convols2(p,:,2)];
              else
                  tent1  = [desc(:,:,1)'; convols1(p,:,1); convols2(p,:,1)];
                  tent2  = [desc(:,:,2)'; convols1(p,:,2); convols2(p,:,2)];
              end
          end
          
          D = VL_ALLDIST2(tent1,tent2);
          [~,indx] = min(D,[],2); 
          correct = (indx == inord');
          tabmatch = tabulate(tryseqind);
          seqmatch = tryseqind(correct);
          tabseqmat = tabulate(seqmatch);
          cormatch = sum(correct);
          matches(p) = cormatch;
          %AUC
          tmatch_nn = exp(round(linspace(0,log(size(D(:),1)),14)));
          [~,IndD] = sort(D(:));
          [~, index] = sort(IndD);
          Dind = reshape(index,size(D)); clear IndD;
          for g=1:length(tmatch_nn)
              dx=(Dind<tmatch_nn(g));
              wx=diag(dx); 
              cmatch(g) = sum(wx);
              match_dist(g) = max(wx.*diag(D));
              matx = (D <= match_dist(g));
              matchtot(g) = sum(sum(matx));
          end

          recall = cmatch/size(D,1);
          precision = cmatch./matchtot;
          AUPR = trapz(1-precision(~(isnan(precision))),recall(~(isnan(precision))));
          
          if ((size(tabseqmat,2)==3) & (size(tabseqmat,1)==8))
            seqrate = tabseqmat(:,2)./tabmatch(:,2);
            seqsucc = sort(seqrate);
            mratio = (seqsucc(1)+seqsucc(2))/(seqsucc(end-1)+seqsucc(end));
            switch(selectfac)
                case 1
                    attrac(p) = cormatch * mratio * mratio; % ent(tabseqmat(:,3)/100) * ent(tabseqmat(:,3)/100);   %AUPR * AUPR * AUPR 
                case 2
                    attrac(p) = cormatch * AUPR * mratio * mratio;
                case 3
                    attrac(p) = AUPR * mratio * mratio;
            end
          else
            attrac(p) = 0;
          end
       end
       
       matchrate = matches/notrainpts; 
       attracrate = attrac;  %/notrainpts; 
       
   i
   max(matchrate)
   [maxmatch, maxindx] = max(attracrate)
 %  maxindx = randind(maxindx);
   if (regions ==1)
      bestfilters(i,:) = useFilters(maxindx,:); %randind(maxindx)
   elseif (regions == 2)
      bestfilters1(i,:) = useFilters1(maxindx,:); 
      bestfilters2(i,:) = useFilters2(maxindx,:); 
   end
   bestparams(i,:) = useParams(maxindx,:);
   bestind(i) = maxindx;  %randind(maxindx);
   
   clear pts corpts  
end

dlmwrite(savefile1, bestfilters1, 'delimiter', ' ', 'precision', 4);
dlmwrite(savefile2, bestfilters2, 'delimiter', ' ', 'precision', 4);
bestfilters1 = reshape(bestfilters1,[size(bestfilters1,1) hsize fsize]);
bestfilters2 = reshape(bestfilters2,[size(bestfilters2,1) (hsize+1) fsize]);
bestfilters(:,1:hsize,:) = bestfilters1;
bestfilters(:,hsize+1:fsize,:) = bestfilters2;
dlmwrite(savefile, bestfilters, 'delimiter', ' ', 'precision', 4);
dlmwrite(parfile, bestparams, 'delimiter', ' ', 'precision', 4);
%readFiltersandTest();
%disp('done!')

end

% do same as line above for rotation 
% generate final descriptors
% generate final graphs
% write up results
% in paper, talk more about training and relationship to existing
% filters/families (with diagram eg. gradient)


