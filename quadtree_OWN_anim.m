clc; close all;

i1 = 'tanks.png';
i2 = 'street.jpg';
%i3 = 'girls.jpg';
i4 = 'toy4.png';
i5 = 'toy.png';
i6 = 'bench1.jpg';
i7 = 'bench2_1.jpg';
i8 = 'jellys2.tiff';
i9 = 'house.tiff';
i10 = 'peppers.tiff';
i11 = 'verytoy.png';
i12 = 'sign.png';

% parameters
predParam = [0.1,0,260]; %  input to predicate(), std,mean_low,mean_high
splitThresh = 0.1;
minMerge = 1;
mindim = 8;

file_path = i12;
f = imread(file_path);
orig = f;
f = rgb2gray(f);
q = 2^nextpow2(max(size(f)));
[m, n]=size(f);
% by default adds zeros, 'post' means add after
% 1 is larger size of image rounded to next pow of 2 i.e. 2,4,8,16,32
% forces f to be of size: pow(2,s) x pow(2,s), 
% where 's' is larger size rounded to next pow of 2
% e.g.: 3x5->8x8, 9x10->16x16
f = padarray(f,[q-m,q-n],'post');
orig = padarray(orig,[q-m,q-n],'post');

% nice example: 
% https://www.mathworks.com/help/images/ref/qtdecomp.html
s = qtdecomp(f, splitThresh, mindim);
fprintf('\nno. segments after only splitting = %d\n' , ...
        nnz(s)) % nnz - no. non-zero values in sparse mtrx
    

g = zeros(size(f));
marker = zeros(size(f));
labels = zeros(size(f));
labels_colors = zeros(size(f,1), size(f,2), 3);

p2 = @(x) 2^x;
prev_max = 0;
random_colors = [];
for k = 1:10;
    k = p2(k);
    [vals, r0, c0] = qtgetblk(f,s,k);
    for i=1:length(r0);
        xlow = r0(i);ylow=c0(i);regsize = s(xlow,ylow);
        labels(xlow:xlow+regsize-1, ylow:ylow+regsize-1) = i+prev_max;
        
        rand_color = rand(1,3)*255;
        random_colors(i+prev_max,1:3) = rand_color;
        labels_colors(xlow:xlow+regsize-1, ylow:ylow+regsize-1, 1) = rand_color(1);
        labels_colors(xlow:xlow+regsize-1, ylow:ylow+regsize-1, 2) = rand_color(2);
        labels_colors(xlow:xlow+regsize-1, ylow:ylow+regsize-1, 3) = rand_color(3);

    end
    prev_max = prev_max + length(r0);
end
figure, imshow(uint8(labels*255/max(max(labels))));

fprintf('length unique = %d\n length(s) = %d', [length(unique(labels)) length(s)])
% 
% x=[]; x(1,5);

fconst = figure;
fconst2 = figure;
pic_iter = 0;
for rep = 1:10;
    labels_begin = labels;
    
    fprintf('\n === rep=%d === \n', rep);
    len = length(unique(labels));
    %for label_1 = randsample(unique(labels), len)';
    for label_1 = unique(labels)';
      fprintf('L1: %d\n', label_1);
      if ~ismember(label_1, unique(labels)); continue; end
      len = length(unique(labels));
      %for label_2 = randsample(unique(labels), len)';
      for label_2 = unique(labels)';
        %fprintf('L1: %d, L2: %d \t', [label_1 label_2]);
        %fprintf('k=%d, L1=%d, L2=%d\n', [k, label_1, label_2])
        if label_1 == label_2; continue; end
        %fprintf('(xlow, ylow) = (%3.2f, %3.2f) \t (next_xlow, next_ylow)=(%3.2f, %3.2f)', [xlow,ylow,next_xlow,next_ylow])
        region = zeros(size(f));
        next_region = zeros(size(f));
        region(labels == label_1) = f(labels == label_1);
        next_region(labels == label_2) = f(labels == label_2);
        %fprintf('labels: %d \t %d \n', [label_1   label_2])

        % length(bwboundaries(..)) == no. separate areas
        % length(bwboundaries([1 1 1; 0 0 1; 0 1 1])) = 1
        % length(bwboundaries([1 1 1; 0 0 0; 0 1 1])) = 2
        are_neighbours = (length(bwboundaries(region + next_region,8, 'noholes')) < 2); % <==> are somehow connected
        if are_neighbours;
            flag = feval(@predicate_OWN,region,next_region);
            if flag
                labels(labels == label_2) = label_1;
                idx = labels == label_2;
                labels_colors(:,:,:) = zeros(size(labels_colors));
                for label = unique(labels)';
                   idx = labels == label;
                   for i_color = 1:3;
                        labels_colors(:,:,i_color) = labels_colors(:,:,i_color) + idx * random_colors(label,i_color);
                   end
                end

%                  if sum(sum(labels == label_1)) > sum(sum(labels == label_2)) % region_1 > region_2
%                      labels(labels == label_2) = label_1;
%                  else
%                      labels(labels == label_1) = label_2;
%                  end
                 %fprintf('merging')
                 %figure(), imshow(uint8(g*255)),title(['k=' num2str(k) ', i=' num2str(i)]);
                 %ff = figure(); 
                 %figure(fconst)
                 %imshow(uint8(labels*255/max(max(labels))));
                 set(fconst2,'visible','off');

                 %figure(fconst2)
                 imshow(uint8(labels_colors));
                 saveas(fconst2, ['pics/sign/anim2_' num2str(pic_iter) '.png']);
                 pic_iter = pic_iter + 1;
                 %fprintf('figure number: %d\n', ff);
            end

        end
      end
    end
    if (labels == labels_begin); break;  end
end

g = labels;
% g = bwlabel(imreconstruct(marker, g));
gc = zeros(size(g,1), size(g,2), 3);
gc_real = zeros(size(g,1), size(g,2), 3);

for label = unique(g)';
    label
    g == label;
    idx = (g == label);
    gray_color = uint8(mean(f(idx)));
    g(idx) = gray_color;

    random_color = rand(3,1)*255;
    
    real_R = mean(mean(orig(:,:,1).*uint8(idx)));
    real_G = mean(mean(orig(:,:,2).*uint8(idx)));
    real_B = mean(mean(orig(:,:,3).*uint8(idx)));
    real_color = [real_R, real_G, real_B] * size(idx,1)*size(idx,2)/sum(sum(idx));
    for ii = 1:3; 
        tmp = zeros(size(g,1));
        tmp(idx) = random_color(ii);
        gc(:,:,ii) = gc(:,:,ii) + tmp;
        
        tmp = zeros(size(g,1));
        tmp(idx) = real_color(ii);
        gc_real(:,:,ii) = gc_real(:,:,ii) + tmp;
    end
end
%g = g(1:m,1:n);
fprintf('no. segments after merging gray  = %d\n' , length(unique(g)) );
fprintf('no. segments after merging color = %d\n' , length(unique(labels)) );


g = uint8(g);
gc = uint8(gc);
gc_real = uint8(gc_real);
f = f(1:m,1:n);
figure, imshow(f),title('Original Image');
fig_gray = figure(); imshow(g);%,title(['Segmented Image; ' num2str(predParam(1))]);
%figure; imhist(g); set(gca, 'YScale', 'log');
fig_rand = figure(); imshow(gc),title(['Segmented Image in random colors; ' num2str(predParam(1))]);
fig_real = figure(); imshow(gc_real),title(['Segmented Image in original colors; ' num2str(predParam(1))]);
saveas(fig_gray, 'pics/sign/anim2_ending.png')
% if (length(unique(g)) > 2)
%     saveas(fig_gray, ['im_gray' strrep(strrep(file_path,'.jpg',''), '.png','') ...
%            '_std' num2str(predParam(1)) '_split' num2str(splitThresh) '_minMerge' num2str(minMerge)  '.png']);
%     saveas(fig_rand, ['im_rand' strrep(strrep(file_path,'.jpg',''), '.png','') ...
%            '_std' num2str(predParam(1)) '_split' num2str(splitThresh) '_minMerge' num2str(minMerge)  '.png']);
%     saveas(fig_real, ['im_real' strrep(strrep(file_path,'.jpg',''), '.png','') ...
%            '_std' num2str(predParam(1)) '_split' num2str(splitThresh) '_minMerge' num2str(minMerge)  '.png']);
% end