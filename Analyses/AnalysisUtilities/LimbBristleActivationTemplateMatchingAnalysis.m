%%

[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningPaperColormaps();

%% Compute the coherence and phase for each canonical gait mode at each timepoint


% Define list of gaits
gaitList = {'tripod','left tetrapod','right tetrapod','wave'};

% Define offsets for each canonical gait
psi = [
    0, 1/2, 0, 1/2, 0, 1/2;
    1/3, 2/3, 0, 0, 1/3, 2/3;
    2/3, 0, 1/3, 0, 1/3, 2/3;
    1/6, 1/3, 1/2, 2/3, 5/6, 0;
    ];

%% Iterate over canonical gaits
for gaitInd = 1:size(psi,1)
    
    % Compute the resultant
    phiBar = squeeze(nanmean(exp(1i* (Phisym + 2*pi*psi(gaitInd,:))),2));
    
    % Compute the modulus and global phase phase
    tripodCoherence = abs(phiBar);
    tripodPhase = angle(phiBar);
    
    %% Plot the distribution of tripod coherences at each timepoint
    
    xq = 0:0.01:1;
    fPDFr = nan(length(xq),length(t_ms),size(hitSym,2));
    for indH = 1:size(hitSym,2)
        idx = hitSym(:,indH) & vfIdx;
        for ind = 1:length(t_ms)
            fPDFr(:,ind,indH) = ksdensity(tripodCoherence(ind,idx), xq, 'kernel','epanechnikov');
        end
    end
    
    for indH = 1:size(hitSym,2)
        % Plot as contour
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        imagesc(t_ms,xq,fPDFr(:,:,indH));
        hold on;
        contour(t_ms,xq,fPDFr(:,:,indH),0:0.5:6, 'EdgeColor','k', 'linewidth',1);
        axis('xy','square','tight');
        ylim([0 1]);
        xlabel('time (ms)');
        ylabel(sprintf('coherence of %s global phase',gaitList{gaitInd}));
        ConfAxis('fontSize', 16);
        cbar = colorbar;
        ylabel(cbar, 'pdf(t)');
        caxis([0 4]);
        colormap(cmpRed);
        title(legendStr{indH});
    end
    
    %% Plot the distribution of tripod phases at each timepoint
    
    xq = (0:1/180:1)';
    
    fPDFpsi = nan(length(xq),length(t_ms),size(hitSym,2));
    for indH = 1:size(hitSym,2)
        idx = hitSym(:,indH) & vfIdx;
        for ind = 1:length(t_ms)
            fPDFpsi(:,ind,indH) = bqksdensity(mod(tripodPhase(ind,idx),2*pi), 2*pi*xq);
        end
    end
    
    for indH = 1:size(hitSym,2)
        
        % Plot as contour
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        imagesc(t_ms,xq,fPDFpsi(:,:,indH));
        hold on;
        contour(t_ms,xq,fPDFpsi(:,:,indH),0:0.05:0.5, 'EdgeColor','k');
        
        axis('xy','square','tight');
        ylim([0 1]);
        xlabel('time (ms)');
        ylabel(sprintf('global %s phase (cycles modulo 1)',gaitList{gaitInd}));
        ConfAxis('fontSize', 16);
        cbar = colorbar;
        ylabel(cbar, 'pdf(t) (1/cycles)');
        caxis([0 0.25]);
        cbar.Ticks = 0:0.05:0.5;
        colormap(cmpRed);
        title(legendStr{indH});
    end
    
    %% Locate the peak yaw rate
    
    vRpeakFrame = nan(size(V,3),1);
    vRpeakTime = nan(size(V,3),1);
    vRpeakValue = nan(size(V,3),1);
    
    phiPeakTurn = nan(size(V,3),6);
    tripodPhasePeakTurn = nan(size(V,3),1);
    idx = t>0;
    
    tPos = t_ms(idx);
    
    tic;
    for ind = 1:size(V,3)
        
        v = squeeze(Vsym(idx,1,ind));
        [~,l] = findpeaks(v, 'SortStr','descend');
        l = l(1);
        vRpeakFrame(ind) = l;
        vRpeakTime(ind) = tPos(l);
        vRpeakValue(ind) = v(l);
        
        p = Phisym(idx,:,ind);
        phiPeakTurn(ind,:) = p(l,:);
        p = tripodPhase(idx,ind);
        tripodPhasePeakTurn(ind) = p(l);
    end
    toc;
    
    %% Compute PDFs of tripod global phase
    
    % Set edges for discretization
    xq = (0:1/180:1)';
    
    % Iterate over limbs
    fPDF = zeros(length(xq),3);
    fPDFci = nan(length(xq),3,2);
    
    haP = nan(3,1);
    haM = nan(6,3);
    haN = nan(6,3);
    
    for indH = 1:3
        idx = hitSym(:,indH) & ~isnan(tripodPhasePeakTurn);
        fPDF(:, indH) = bqksdensity(mod(tripodPhasePeakTurn(idx),2*pi), 2*pi*xq, 1);
        
        ci = bootci(nboot, {@(x)  bqksdensity(mod(x,2*pi), 2*pi*xq, 1), tripodPhasePeakTurn(idx)});
        fPDFci(:,indH,:) = permute(ci, [2,3,1]);
    end
    
    MakeFigure;
    PlotAsymmetricErrorPatch(xq, fPDF,fPDFci(:,:,1),fPDFci(:,:,2), corder);
    ylim([0 0.3]);
    xlabel(sprintf('%s phase at yaw rate extremum (cycles)',gaitList{gaitInd}));
    ylabel('pdf (1/cycles)');
    legend(legendStr, 'location','eastoutside');
    yticks(0:0.1:0.3);
    axis('square');
    ConfAxis('fontSize', 16);
    
end