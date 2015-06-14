clear all
close all
fps = 30;
filename = ['bigmovie'];
aviobj = VideoWriter(filename,'MPEG-4');  
aviobj.FrameRate = fps;
aviobj.Quality = 100;
open(aviobj);
maxpop = 5; %Number of data-files in data/population
figure('Position',[1 0 800 600],'MenuBar','none','ToolBar','none','resize','off') % fullscreen
for k = 1:maxpop
    alldata = importdata(['data/population/data' num2str(k-1) '.dat']);
    % t, x, dx, th, dth, Fx, L, fps, arr->errorTh[0]);
    t  = alldata(:,1);
    L  = alldata(1,7);
    x1 = alldata(:,2);
    x2 = L*sin(alldata(:,4))+x1;
    y1 = zeros(length(x1));
    y2 = -L*cos(alldata(:,4));
    gen = alldata(1,10);
    set(gca,'FontSize',20)

    for i = 1:length(x1)
        h = plot([x1(i) x2(i)],[y1(i) y2(i)],'Linewidth',4); hold on
        plot(x1(i)-0.07, -0.03, '.k', 'MarkerSize',69),hold on
        plot(x1(i)+0.07, -0.03, '.k', 'MarkerSize',69),hold off
        rectangle('Position',[x1(i)-0.1 -0.0325 0.2 0.0325], 'FaceColor','r')
        axis([-0.8 0.8 -0.5 1])
        grid on;
        xlabel('Displacement [m]','FontSize', 20);
        ylabel('Height [m]','FontSize', 20)
        title({'Evolutionary Monte Carlo'; ['Generation ' num2str(gen)]},'FontSize', 24, 'fontweight','bold');
        movegui(h, 'onscreen');
        %hold all;
        drawnow;
        % Write each frame to the file.
        %currFrame = getframe;
        
        currFrame = getframe( gcf );
        writeVideo(aviobj,currFrame);
    end

end
close(aviobj);  
close all