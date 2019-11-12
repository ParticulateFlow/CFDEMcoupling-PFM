fname='temp_histo_fine.txt';

% number of header lines
nheader=4;

inputarr = dlmread(fname,'' ,nheader,0);

% number of time steps
ntimes = 1000;

% number of bins
nbins = 224;
dt = 1.0;

lowerbound = 1;
upperbound = lowerbound + nbins - 1;

fid0 = fopen(strcat('stat_',fname), 'w');
fprintf(fid0,'# time index   ||   ave    ||    std    ||   maximum   ||   FWHF\n');

for i=0:ntimes

  maxval = 0.0;
  for j= lowerbound:upperbound
    if(inputarr(j,3) > maxval)
      maxval = inputarr(j,3);
    end
  end

  hm = 0.5*maxval
  lower = 0.0;
  upper = 0.0;

  mean = 0.0;
  meansqr = 0.0;
  counter = 0;

  for j = lowerbound:upperbound-1
    if(inputarr(j,3) < hm && inputarr(j+1,3) >= hm)
       lower = 0.5*(inputarr(j,2)+inputarr(j+1,2));
j
hm
inputarr(j,3)
inputarr(j+1,3)
lower
      % lower = inputarr(j,2) + (inputarr(j+1,2) - inputarr(j,2)) * (hm - inputarr(j,3)) / (inputarr(j+1,3) - inputarr(j,3));
    end
    if(inputarr(j,3) > hm && inputarr(j+1,3) <= hm)
      upper = 0.5*(inputarr(j,2)+inputarr(j+1,2));
      %upper = inputarr(j+1,2) - (inputarr(j+1,2) - inputarr(j,2)) * (hm - inputarr(j+1,3)) / (inputarr(j,3) - inputarr(j+1,3));
    end

    mean = mean + inputarr(j,2) * inputarr(j,3);
    meansqr = meansqr + inputarr(j,2)^2 * inputarr(j,3);
    counter = counter + inputarr(j,3);
  end

  mean = mean/double(counter);
  std = meansqr/double(counter) - mean^2;
  std = sqrt(std);
  fwhm = upper - lower;
upper
lower

  if(i>0)
    fprintf(fid0,'%f\t%f\t%f\t%f\t%f\n',i*dt,mean,std,maxval,fwhm);
  end


  lowerbound = upperbound + 2;
  upperbound = lowerbound + nbins - 1;

end

fclose(fid0);
