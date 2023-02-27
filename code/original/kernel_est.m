function [dx dy c ] =  kernel_est(I_in) 

  I_in = rgb2gray(I_in); % change to gray scale image
  Laplacian=[0 -1 0; -1 4 -1; 0 -1 0]; % Laplacian operator
  resp = imfilter(I_in, Laplacian); % Laplacian I
  auto_corr = xcorr2(resp, resp); % cross correlation of matrices
  bdry = 370; %£¿
  auto_corr = auto_corr(bdry:end-bdry, bdry:end-bdry); % cut edge? 
  
  max_1 = ordfilt2(auto_corr, 25, true(5)); % max filter
  max_2 = ordfilt2(auto_corr, 24, true(5)); % second max filter

  auto_corr(end/2 - 4 : end/2 + 4, end/2 - 4 : end/2+4)=0; % remove local maxima within 4 pixels of the origin
  candidates = find((auto_corr == max_1) & ((max_1 - max_2)>70));
  candidates_val = auto_corr(candidates);

  cur_max = 0;
  dx = 0; 
  dy = 0;
  offset = size(auto_corr)/2 + 1; % starts from the center of autocorr matrix
  for i = 1 : length(candidates)
    if (candidates_val(i) > cur_max)  
      [dy dx] = ind2sub(size(auto_corr), candidates(i)); 
      dy = dy - offset(1);
      dx = dx - offset(2);
    end
  end
  c = est_attenuation(I_in, dx, dy);
