function saveFigures(TITLE,savePath,inputMode,FIG_MODE)

if strcmp(FIG_MODE,'show') || strcmp(FIG_MODE,'hide')
    
    path = [savePath,'/figures'];
    if exist(path,'dir') ~= 7
        mkdir(path);
    end
    
    AA = load([savePath,'/results.mat']);
    summary = AA.summary;
    TYPE = AA.setting_alignment.data_type;
    if ~strcmp(TYPE,'d18O')
        CI_C14 = AA.CI_C14;
    end
    target = AA.target;
    
    clear AA;
    
    % returns n distinguishable colors for stack summary plot
    ncores = length(summary);
    colors = distinguishable_colors(ncores);
    cc(1,:) = [0 0 1]; % radiocarbon
    cc(2,:) = [0.6350 0.0780 0.1840]; % d18O
    cc(3,:) = [0.4660 0.6740 0.1880]; % additional ages
    cc(4,:) = [0.4940, 0.1840, 0.5560];% test
    
    fs = 14; % figure fontsize
    n = 2; % # error bars for d18O plot
    
    if isfield(target,'init_stack')
        istack = target.init_stack;
    end
    fstack = target.stack;
    
    if strcmp(inputMode,'stack_construction') || strcmp(inputMode,'stacking')
        figure(100); % stack summary figure
        hold on;
        shadebetweenlines(fstack(:,1),fstack(:,2)+1.96*fstack(:,3),fstack(:,2)-1.96*fstack(:,3),[0 0 0]);
        
        if isfield(target,'init_stack')
            plot(istack(:,1),istack(:,2),'k','handlevisibility','off');
            plot(istack(:,1),istack(:,2)+1.96*istack(:,3),'--k','handlevisibility','off');
            plot(istack(:,1),istack(:,2)-1.96*istack(:,3),'--k','handlevisibility','off');
        end
    end
    
    
    for i = 1:length(summary)
        as = summary(i).age_samples;
        d = summary(i).depth;
        med = summary(i).median;
        up = summary(i).upper_95;
        low = summary(i).lower_95;
        CI95 = up-low;
        if ~strcmp(TYPE,'d18O')
            cd = CI_C14(i).depth;
            cmed = CI_C14(i).median;
            l68 = CI_C14(i).lower_68;
            u68 = CI_C14(i).upper_68;
            l95 = CI_C14(i).lower_95;
            u95 = CI_C14(i).upper_95;
        end
        d18O = summary(i).d18O;
        shift = summary(i).d18O_shift;
        scale = summary(i).d18O_scale;
        aa = summary(i).additional_ages;
        
        fig = figure; % plot age vs. depth
        if sum(~isnan(d18O)) == 0
            h2 = subplot(1,3,[1 2]);
        else
            h2 = subplot(2,3,[1 2]);
        end
        hold on;
        title(summary(i).name);
        plot(as,d,'color',[0 0 0 0.015]);
        
        %radiocarbon ages
        if ~strcmp(TYPE,'d18O')
            for j = 1:length(cd)
                h = h2.YLim(2)*0.005;% scales the height of radiocarbon ages by figure size
                y = cd(j)-h/2; % centers rectangle
                x = l68(j);
                w = u68(j)-x; % centers rectangle
                rectangle('position',[x,y,w,h],'edgecolor','k','facecolor',cc(1,:));
                plot([l95(j) u95(j)],[cd(j) cd(j)],'color',cc(1,:),'linewidth',1);
            end
        end
        
        % additional ages
        
        indg = find(aa(:,3) == 1); % Gaussian
        indu = find(aa(:,3) == 0); % Uniform
        for j = 1:length(indg) % plot Gaussian
            h = h2.YLim(2)*0.005; % scale height by depth of core to keep consistent
            y = d(indg(j))-h/2;
            x = aa(indg(j),1)-aa(indg(j),2);
            w = 2*aa(indg(j),2);
            sig2l = aa(indg(j),1)-2*aa(indg(j),2);
            sig2u = aa(indg(j),1)+2*aa(indg(j),2);
            
            rectangle('position',[x y w h],'facecolor',cc(3,:),'edgecolor','k');
            plot([sig2l sig2u],[d(indg(j)) d(indg(j))], 'k');
        end
        
        for j = 1:length(indu) % plot uniform
            h = h2.YLim(2)*0.005;
            y = d(indu(j))-h/2;
            x = aa(indu(j),1)-aa(indu(j),2);
            w = 2*aa(indu(j),2);
            rectangle('position',[x y w h],'facecolor',cc(3,:),'edgecolor','k');
        end
        
        % plot 95upper, 95lower, and median age model
        plot(med,d,'r');
        plot(low,d,'--k');
        plot(up,d,'--k');
        
        xlabel('Age (kyr)');
        ylabel('Depth (m)');
        xlim([0 max(as,[],'all')]); % set xlim to max of age model samples
        set(gca,'fontsize',fs);
        grid on;
        
        % plot confidence interval width
        if sum(~isnan(d18O)) == 0
            subplot(1,3,3);
        else
            subplot(2,3,3);
        end
        hold on;
        if strcmp(TYPE,'C14')
            title('14C Only');
        elseif strcmp(TYPE,'d18O')
            title('d18O Only');
        elseif strcmp(TYPE,'both')
            title('dual proxy');
        end
        plot(CI95,d,'k','linewidth',2);
        if ~strcmp(TYPE,'d18O')
            scatter(zeros(length(cd),1),cd,'>','markerfacecolor',cc(1,:),'markeredgecolor','k','handlevisibility','off');
        end
        ind = find(~isnan(aa(:,1)));
        scatter(zeros(length(d(ind)),1),d(ind),'square','markerfacecolor',cc(3,:),'markeredgecolor','k','handlevisibility','off');
        xlabel('95% CI (kyr)');
        set(gca,'fontsize',fs);
        grid on;
        
        % plot shifted and scaled d18O
        if sum(~isnan(d18O)) == 0
            set(gcf,'position',[200 200 800 600],'color','w');
        else
            h = subplot(2,1,2);
            hold on;
            column = combine_columns(d, d18O);
            column(:,2) = (column(:,2)-shift)/scale;
            ages_median = interp1(d,med,column(:,1));
            ages_lower = interp1(d,low,column(:,1));
            ages_upper = interp1(d,up,column(:,1));
            shadebetweenlines(fstack(:,1),fstack(:,2)+1.96*fstack(:,3),fstack(:,2)-1.96*fstack(:,3),[0 0 0]);
            
            plot(fstack(:,1),fstack(:,2),'k');
            % plot error bars for every nth data point
            for j = 1:n:length(column(:,2))
                plot([ages_lower(j) ages_upper(j)],[column(j,2) column(j,2)],'r');%'color',cc(2,:))
            end
            scatter(ages_median,column(:,2),'*','markeredgecolor',cc(2,:));
            plot(ages_median,column(:,2),'color',cc(2,:));
            
            
            ylim(h.YLim);
            if ~strcmp(TYPE,'d18O')
                scatter(cmed,h.YLim(2)*ones(length(cmed),1),'^','markerfacecolor',cc(1,:),'markeredgecolor','k','handlevisibility','off');
            end
            scatter(aa(:,1),h.YLim(2)*ones(length(aa(:,1)),1),'square','markerfacecolor',cc(3,:),'markeredgecolor','k','handlevisibility','off');
            set(gca,'ydir','reverse','fontsize',fs);
            xlabel('Age (kyr)');
            ylabel(['\delta^1^8O (',char(8240),')']);
            title(['shift = ',num2str(round(shift,2)),' and scale = ',num2str(round(scale,2))]);
            xlim([0 max(as,[],'all')]);
            set(gcf,'position',[200 200 800 900],'color','w');
        end
        
        grid on;
        saveas(gcf,[savePath,'/figures/',summary(i).name,'.png']);
        
        movegui(fig,'center');
        
        
        
        % stack summary figure
        if strcmp(inputMode,'stack_construction') || strcmp(inputMode,'stacking')
            figure(100);
            h = gca;
            hold on;
            scatter(ages_median,column(:,2),'*','markeredgecolor',colors(i,:));
            if ~strcmp(TYPE,'d18O')
                scatter(cmed,h.YLim(2)*ones(length(cmed),1),'^','markerfacecolor',colors(i,:),'markeredgecolor','k','handlevisibility','off');
            end
            scatter(aa(:,1),h.YLim(2)*ones(length(aa(:,1)),1),'square','markerfacecolor',colors(i,:),'markeredgecolor','k','handlevisibility','off');
        end
    end
    %
    if strcmp(inputMode,'stack_construction') || strcmp(inputMode,'stacking')
        figure(100);
        hold on;
        title(TITLE);
        set(gca,'ydir','reverse','fontsize',fs);
        legend(summary.name,'location','northeastoutside');
        xlabel('Age (kyr)');
        ylabel(['\delta^1^8O (',char(8240),')']);
        grid on;
        set(gcf,'position',[200 200 1200 800],'color','w');
        saveas(gcf,[savePath,'/figures/stack_summary.png']);
        
        movegui(figure(100),'center');
    end
end


if strcmp(FIG_MODE,'hide')
    close all;
end


end