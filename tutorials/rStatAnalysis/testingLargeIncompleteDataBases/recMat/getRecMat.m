% open files
fid1 = fopen('../recurrenceMatrix');

% read dimensions
A1 = fscanf(fid1, '%g %g');

% skip the first two lines
tline = fgetl(fid1);
tline = fgetl(fid1);

% get dimensions
N1 = A1(1)
M1 = A1(2)

% allocate space
B0 = zeros(M1,N1);

% read data
for i=1:N1
    B0(:,i) = fscanf(fid1, '%g', inf);
    tline = fgetl(fid1);
    tline = fgetl(fid1);
end

% close files
fclose(fid1);


% skip this many leading entries
sle = 0;

B1 = zeros(M1-sle,N1-sle);

B1 = B0(1+sle:M1,1+sle:N1);

maxval=0.0;
%for i=1:M1
%  for j=1:N1
%    if(B1(i,j)>maxval)
%      maxval=B1(i,j);
%    endif
%  end
%end

maxval = 1.0;
for i=1:M1-sle
  for j=1:N1-sle
    B1(i,j)=1-B1(i,j)/maxval;
   % B1(i,j)=B1(i,j)/maxval;
  end
end



% write full matrix to simple text file
dlmwrite('myMatrix.txt',B1,'delimiter','\t','precision',3)



%plot(C1)
%saveas(gcf,'Plot','png')

%B1=B1*1;


%hold on


%colormap(jet(50))
%imagesc(B1)
%colorbar

%saveas(gcf,'Figure','png')


