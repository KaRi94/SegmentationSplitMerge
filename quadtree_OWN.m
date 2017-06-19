clc; close all;

i1 = 'tanks.png';
i2 = 'street.jpg';
%i3 = 'girls.jpg';
i4 = 'toy4.png';
i5 = 'toy.png';
i6 = 'bench1.jpg';
i7 = 'baboon.jpg';
i8 = 'jellys.tiff';
i9 = 'house.tiff';
i10 = 'peppers.tiff';
i11 = 'verytoy.png';
i12 = 'sign.png';

% parameters
predParam = [0.1,0,260]; %  input to predicate(), std,mean_low,mean_high
splitThresh = 0.5;
minMerge = 1;
mindim = 4;

file_path = i9;
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

p2 = @(x) 2^x;
prev_max = 0;
random_colors = [];
for k = 1:10;
    k = p2(k);
    [vals, r0, c0] = qtgetblk(f,s,k);
    for i=1:length(r0);
        xlow = r0(i);ylow=c0(i);regsize = s(xlow,ylow);
        labels(xlow:xlow+regsize-1, ylow:ylow+regsize-1) = i+prev_max;
    end
    prev_max = prev_max + length(r0);
end
figure, imshow(uint8(labels*255/max(max(labels))));


fconst = figure;
fconst2 = figure;
pic_iter = 0;
for rep = 1:10;
    labels_begin = labels;
    
    fprintf('\n === rep=%d === \n', rep);
    len = length(unique(labels));
    for label_1 = unique(labels)';
      fprintf('L1: %d\n', label_1);
      if ~ismember(label_1, unique(labels)); continue; end
      len = length(unique(labels));
      for label_2 = unique(labels)';
        if label_1 == label_2; continue; end
        region = zeros(size(f));
        next_region = zeros(size(f));
        region(labels == label_1) = f(labels == label_1);
        next_region(labels == label_2) = f(labels == label_2);

        % length(bwboundaries(..)) == no. separate areas
        % length(bwboundaries([1 1 1; 0 0 1; 0 1 1])) = 1
        % length(bwboundaries([1 1 1; 0 0 0; 0 1 1])) = 2
        are_neighbours = (length(bwboundaries(region + next_region,8, 'noholes')) < 2); % <==> are somehow connected
        if are_neighbours;
            flag = feval(@predicate_OWN,region,next_region);
            if flag
                labels(labels == label_2) = label_1;
                
%                  figure(fconst2)
%                  imshow(uint8(labels*255/max(max(labels))));

            end
        end
      end
    end
    if (labels == labels_begin); break;  end
end

g = labels;
gc = zeros(size(g,1), size(g,2), 3);
gc_real = zeros(size(g,1), size(g,2), 3);

for label = unique(g)';
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
fprintf('no. segments after merging gray  = %d\n' , length(unique(g)) );
fprintf('no. segments after merging color = %d\n' , length(unique(labels)) );


g = uint8(g);
gc = uint8(gc);
gc_real = uint8(gc_real);
f = f(1:m,1:n);
figure, imshow(f),title('Original Image');
fig_gray = figure(); imshow(g),title(['Segmented Image; ' num2str(predParam(1))]);
fig_rand = figure(); imshow(gc),title(['Segmented Image in random colors; ' num2str(predParam(1))]);
fig_real = figure(); imshow(gc_real),title(['Segmented Image in original colors; ' num2str(predParam(1))]);
%saveas(fig_gray, 'pics/sign/gray1_result.png');
%saveas(fig_real, 'pics/sign/color1_result.png');
