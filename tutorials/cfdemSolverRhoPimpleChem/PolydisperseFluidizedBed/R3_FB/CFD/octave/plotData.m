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
dataexp = importdata('R3_experiment.dat', ' ');

% init figures
hFig(1) = figure;

cmap = colormap(jet(16));
linS = {'-';'--';'-.';':';'-';'--';'-.';':';'-';'--';'-.';':'};
markers = {'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};

for ii=1:nFiles

    fname = data(ii).name;
    fbasename = strtrunc(fname,index(fname,".dat")-1);
    stepsize = 2;
    timesteps = data(ii).values(1:stepsize:end,col_t)*10;
    nColumns = columns(data(ii).values);

    figure(hFig(1));
    clf reset;
    hold on;

    for jj=2:nColumns
        xvalue = data(ii).values(1:stepsize:end,jj);
        %if (strncmp(fbasename,"Aterm",5) || strncmp(fbasename,"Bterm",5))
        %    semilogy(timesteps,xvalue,'Color',cmap(jj,:),'LineStyle',linS{jj-1},'Marker',markers{jj-1},'MarkerSize',5);
        %else
            plot(timesteps,xvalue,'Color',cmap(jj,:),'LineStyle',linS{jj-1},'Marker',markers{jj-1},'MarkerSize',5);
        %endif
    end
    if (strncmp(fbasename,"fr_OV",5))
        timesteps = dataexp.data(1:stepsize:end,col_t);
        xvalue = dataexp.data(1:stepsize:end,2)*0.01;
        plot(timesteps-65,xvalue,'Color','red','LineStyle','none','Marker','.','MarkerSize',10);%,'MarkerFaceColor','auto');
    endif
    xlim ([0.0, 1000.0]);
    ylim auto;
    xlabel('time (s)');
    ylabel(fbasename);
    grid on;

    headerline = substr(data(ii).header{1},7);

    if (strncmp(fbasename,"fr_OV",5))
        headerline = [headerline, " fr_exp_OV"];
    endif

    legend(strsplit(headerline, " "),'location','eastoutside');

    print(hFig(1),fullfile(dirfile,['figure_',fbasename,'.eps']),'-color','-deps');
end


