function AgeVsD18O(results_path)

results = load([results_path,'/results.mat']);
data = results.summary;
target = results.target;
setting = results.setting_alignment;

data_type = setting.data_type;

if ~strcmp(data_type,'C14')
    
    path = [results_path,'/figures/age_vs_d18O'];
    if exist(path,'dir') ~= 7
        mkdir(path);
    end
    
    clear results;
    
    L = length(data);
    
    QQQ = cell(L,3);
    
    for ll = 1:L
        if sum(isnan(data(ll).d18O)) < numel(data(ll).d18O)
            
            
            N = length(data(ll).depth);
            
            AMAX = size(data(ll).d18O,2);
            
            H = data(ll).isoutlier;
            H = (mean(H,2)>0.5);
            H = reshape(H,[N,AMAX]);
            
            m = 0;
            OLP = cell(N,1);
            for n = 1:N
                K = sum(~isnan(data(ll).d18O(n,:)));
                if K > 1
                    m = m + 1;
                    AA = quantile(data(ll).age_samples(n,:),0.5)*ones(1,K);
                    BB = (data(ll).d18O(n,1:K)-data(ll).d18O_shift)/data(ll).d18O_scale;
                    CC = H(n,1:K);
                    [~,order] = sort(BB,'ascend');
                    BB = BB(order);
                    CC = CC(order);
                    OLP{m} = {AA,BB,CC};
                end
            end
            OLP = OLP(1:m);
            
            M = 0;
            for n = 1:N
                if ~isnan(data(ll).d18O(n,1))
                    M = M + 1;
                end
            end
            
            LMU = zeros(M,3);
            d18O = zeros(M,AMAX);
            IsOutlier = zeros(M,AMAX);
            m = 0;
            for n = 1:N
                if ~isnan(data(ll).d18O(n,1))
                    m = m + 1;
                    LMU(m,1) = quantile(data(ll).age_samples(n,:),0.025);
                    LMU(m,2) = quantile(data(ll).age_samples(n,:),0.5);
                    LMU(m,3) = quantile(data(ll).age_samples(n,:),0.975);
                    d18O(m,:) = (data(ll).d18O(n,:)-data(ll).d18O_shift)/data(ll).d18O_scale;
                    % IsOutlier(m,:) = H(n,:);
                end
            end
            
            index = (~isnan(d18O));
            AAA = repmat(LMU(:,1),[1,AMAX]);
            BBB = repmat(LMU(:,2),[1,AMAX]);
            CCC = repmat(LMU(:,3),[1,AMAX]);
            
            AAA = AAA(index);
            BBB = BBB(index);
            CCC = CCC(index);
            LMU = [AAA,BBB,CCC];
            d18O = d18O(index);
            IsOutlier = IsOutlier(index);
            M = size(LMU,1);
            
            QQQ{ll,1} = LMU(:,2);
            QQQ{ll,2} = d18O;
            QQQ{ll,3} = IsOutlier;
            
            MAX = max(LMU(:,3));
            MAX = ceil(MAX/10)*10;
            if isempty(MAX)
                MAX = target.stack(end,1);
            end
            
            
            if isempty(OLP) == 1
                fig = figure;
                hold on;
                tt = [data(ll).name,' (',data_type,')'];
                tt(tt=='_') = '-';
                title(tt,'FontSize',16);
                xx = [target.stack(:,1);flipud(target.stack(:,1))];
                yy1 = [target.stack(:,2)-target.stack(:,3);flipud(target.stack(:,2)+target.stack(:,3))];
                yy2 = [target.stack(:,2)-2*target.stack(:,3);flipud(target.stack(:,2)+2*target.stack(:,3))];
                patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
                patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
                for m = 1:M
                    if IsOutlier(m) < 0.5
                        plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':b','LineWidth',2);
                        plot(LMU(m,2),d18O(m),'*r','LineWidth',2);
                    else
                        plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':r','LineWidth',2);
                        plot(LMU(m,2),d18O(m),'*b','LineWidth',2);
                    end
                end
                xlabel('age (ky)','FontSize',12);
                ylabel('{\delta}^{18}O (‰)','FontSize',12);
                xlim([target.stack(1,1) MAX]);
                
                set(gca,'YDir','rev');
                set(fig,'Position',[10 10 1200 500]);
                movegui(fig,'center');
                
                path = [results_path,'/figures/age_vs_d18O/',data(ll).name,'.fig'];
                savefig(fig,path);
            else
                fig = figure;
                
                subplot(2,1,1);
                hold on;
                tt = [data(ll).name,' (',data_type,')'];
                tt(tt=='_') = '-';
                title(tt,'FontSize',16);
                xx = [target.stack(:,1);flipud(target.stack(:,1))];
                yy1 = [target.stack(:,2)-target.stack(:,3);flipud(target.stack(:,2)+target.stack(:,3))];
                yy2 = [target.stack(:,2)-2*target.stack(:,3);flipud(target.stack(:,2)+2*target.stack(:,3))];
                patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
                patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
                for m = 1:M
                    if IsOutlier(m) < 0.5
                        plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':b','LineWidth',2);
                        plot(LMU(m,2),d18O(m),'*r','LineWidth',2);
                    else
                        plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':r','LineWidth',2);
                        plot(LMU(m,2),d18O(m),'*b','LineWidth',2);
                    end
                end
                xlabel('age (ky)','FontSize',12);
                ylabel('{\delta}^{18}O (‰)','FontSize',12);
                xlim([target.stack(1,1) MAX]);
                
                set(gca,'YDir','rev');
                
                subplot(2,1,2);
                hold on;
                title('Multiple Observations at the Same Depths','FontSize',16);
                xx = [target.stack(:,1);flipud(target.stack(:,1))];
                yy1 = [target.stack(:,2)-target.stack(:,3);flipud(target.stack(:,2)+target.stack(:,3))];
                yy2 = [target.stack(:,2)-2*target.stack(:,3);flipud(target.stack(:,2)+2*target.stack(:,3))];
                patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
                patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
                for m = 1:length(OLP)
                    plot(OLP{m}{1},OLP{m}{2},'--g','LineWidth',2);
                    index = (OLP{m}{3}==0);
                    plot(OLP{m}{1}(index),OLP{m}{2}(index),'*r','LineWidth',2);
                    plot(OLP{m}{1}(~index),OLP{m}{2}(~index),'*b','LineWidth',2);
                end
                xlabel('age (ky)','FontSize',12);
                ylabel('{\delta}^{18}O (‰)','FontSize',12);
                xlim([target.stack(1,1) MAX]);
                
                set(gca,'YDir','rev');
                
                set(fig,'Position',[10 10 1200 1100]);
                movegui(fig,'center');
                
                path = [results_path,'/figures/age_vs_d18O/',data(ll).name,'.fig'];
                savefig(fig,path);
            end
        end
    end
end

close all;


end