#!/usr/bin/octave --silent

%% clear workspace
clear all
close all
clc

dirfile = '../../DEM/post';
filepattern = '*.dat';

% time column in the data matrix
col_t = 1;

% read all data

listfile = dir(fullfile(dirfile,filepattern));
data = readData(dirfile,listfile);
nFiles = length(listfile);

% init figures
hFig(1) = figure;

cmap = colormap(jet(16));
linS = {'-';'--';'-.';':';'-';'--';'-.';':';'-';'--';'-.';':'};
markers = {'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};

for ii=1:nFiles

    fname = data(ii).name;
    fbasename = strtrunc(fname,index(fname,".dat")-1);
    timesteps = data(ii).values(1:1:end,col_t);

    figure(hFig(1));
    clf reset;
    hold on;
    nColumns = columns(data(ii).values);

    for jj=2:nColumns
        xvalue = data(ii).values(1:1:end,jj);
        %if (strncmp(fbasename,"Aterm",5) || strncmp(fbasename,"Bterm",5))
        %    semilogy(timesteps,xvalue,'Color',cmap(jj,:),'LineStyle',linS{jj-1},'Marker',markers{jj-1},'MarkerSize',5);
        %else
            plot(timesteps,xvalue,'Color',cmap(jj,:),'LineStyle',linS{jj-1},'Marker',markers{jj-1},'MarkerSize',5);
        %endif
    end
    xlim auto;
    ylim auto;
    xlabel('time (s)');
    ylabel(fbasename);
    grid on;

    headerline = substr(data(ii).header{1},7);
    legend(strsplit(headerline, " "),'location','eastoutside');

    print(hFig(1),fullfile(dirfile,['figure_',fbasename,'.eps']),'-color','-deps');
end


