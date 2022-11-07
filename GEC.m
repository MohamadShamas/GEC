function gec = GEC(Sig,sWind,sStep,Fs,Fmin,Fmax)
%%%% testing %%%%

LWind = min(round(sWind*Fs),size(Sig,2));
LStep = round(sStep*Fs);
vWind = 0:LStep:(size(Sig,2)-LWind);
PeriWind = round(1/Fmin/2*Fs);
nWind = numel(vWind);
nChan = size(Sig,1);
nPairs = (nChan^2 - nChan)/2;
plot_fig= 0;

% Parameters
nOrder = 32;
if eq(Fs,1000)
    nbins = ceil(1000*1/Fmin)/(1000/(Fs/2))*2;
elseif eq(Fs,2000)
    nbins = (1000*1/Fmin)/(1000/(Fs/2))*2;
else nbins = ceil(Fs/Fmin)*2;
end
disp(['Sampling rate ' num2str(Fs)])
s_ShannonEMax = log(nbins);


% Filter gamma band (30-55Hz)
stt_Filt = fdesign.bandpass('n,f3db1,f3db2',nOrder,Fmin,Fmax,Fs);
h_Filt = design(stt_Filt,'butter');
nPad = min(2*nOrder,size(Sig,2));
SigPad = f_Padding(Sig,nPad,'symm',2);
FiltPad = f_IIRFilter(SigPad',h_Filt)';
SigHilb = hilbert(FiltPad);
SigFilt = FiltPad(:,nPad+1:end-nPad);
SigHilb = SigHilb(:,nPad+1:end-nPad);
SigPhase = angle(SigHilb);
SigAmpl = abs(SigHilb);
SigFilt = mapminmax(SigFilt,-1,1);

% Index to convert connectivity matrix into a connectivity vector
mTriu = triu(ones(nChan),1);
mTriu = find(mTriu);

% Initializing variables
[FConnPairs,FConnAmplPairs] = deal(zeros(nWind,nPairs));
% [FConnMatrix,FConnAmplMatrix,PhSyncMatrix,AmpCorrMatrix] = deal(zeros(nWind,nChan,nChan));
[FConnMatrix,FConnAmplMatrix] = deal(zeros(nWind,nChan,nChan));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iWind = 1:nWind
    
    % Cut variables in segments:
    if iWind < nWind
        WindFilt = SigFilt(:,vWind(iWind)+1:(vWind(iWind)+LWind));
        WindPhase = SigPhase(:,vWind(iWind)+1:(vWind(iWind)+LWind));
    else
        WindFilt = SigFilt(:,vWind(iWind)+1:end);
        WindPhase = SigPhase(:,vWind(iWind)+1:end);
    end
    
    strPairs = deal(cell(nChan));
    mFConn = zeros(nChan);
    mFConnAmpl = zeros(nChan);
    GammaEv = zeros(size(WindFilt));
    mPhaseSync = zeros(nChan);
    
    if plot_fig
    % choose 2 random channels to plot
    ch2plot = randperm(size(WindFilt,1),2);
    ch2plot = sort(ch2plot);
    figure
    ax1 = subplot(5,1,1);
    findpeaks(WindFilt(ch2plot(1),1:LWind),'Threshold',0.1);
    ax1.XGrid = 'on';
    ax1.YGrid = 'off';
    xticks(linspace(0,LWind,nbins));
    xticklabels({});
    
    ax3 = subplot(5,1,3);
    findpeaks(WindFilt(ch2plot(2),1:LWind),'Threshold',0.1);
    ax3.XGrid = 'on';
    ax3.YGrid = 'off';
    xticks(linspace(0,LWind,nbins));
    xticklabels({});
    end
    
    for iChan = 1:size(WindFilt,1)-1 %size is rows or no. of channels
        
        [v_GammaPks,ixPks] = findpeaks(WindFilt(iChan,1:LWind),'Threshold',0.1);
        ixDisc = find((ixPks < PeriWind)|(ixPks > (size(WindFilt,2)-PeriWind))); % Discard peaks with not enough time window before/after
        ixPks(ixDisc) = [];
        GammaEv(iChan,ixPks) = 1;
        
        if plot_fig
        % plot train of events
        if iChan == ch2plot(1)
            ax2 = subplot(5,1,2);
            train_events = zeros(length(WindFilt(iChan,1:LWind)),1);
            train_events(ixPks)= 1;
            plot(ax2,train_events);
            ax2.XGrid = 'on';
            ax2.YGrid = 'off';
            xticks(linspace(0,LWind,nbins));
            xticklabels({});
            xlim(ax2,[0 LWind]);
        end
        end
        for iChan1 = iChan+1:nChan % Loop pairs of channels across rows, i.e. 1-2, 1-3, 1-4,...1-nChan,2 2-3, 2-4,...2-nChan
            
            % Find peaks for all channels in the first iteration
            if iChan == 1
                [v_GammaPks1,ixPks1] = findpeaks(WindFilt(iChan1,1:LWind),'Threshold',0.1);
                ixDisc1 = find((ixPks1 < PeriWind)|(ixPks1 > (size(WindFilt,2)-PeriWind))); % Discard peaks with not enough time window before/after
                v_GammaPks1(ixDisc1) = [];
                ixPks1(ixDisc1) = [];
                GammaEv(iChan1,ixPks1) = 1;
                if plot_fig 
                % plot train of events for second channel
                if iChan1 == ch2plot(2)
                    ax4 = subplot(5,1,4);
                    train_events1 = zeros(length(WindFilt(iChan1,1:LWind)),1);
                    train_events1(ixPks1)= 1;
                    plot(ax4,train_events1);
                    ax4.XGrid = 'on';
                    ax4.YGrid = 'off';
                    xticks(linspace(0,LWind,nbins));
                    xticklabels({});
                    xlim(ax4,[0 LWind]);
                end
                end
            end
            
            
            % Peri-event histogram: matrix of gamma events x time
            [mPeriEv,mPeriEvAmpl] = deal(zeros(numel(ixPks),PeriWind*2));
            
            for iPks2 = 1:numel(ixPks)
                ixSegm = ixPks(iPks2) - PeriWind + 1 : ixPks(iPks2) + PeriWind;
                mPeriEv(iPks2,:)= GammaEv(iChan1,ixSegm);
                mPeriEvAmpl(iPks2,:)= round(GammaEv(iChan1,ixSegm).*abs(WindFilt(iChan1,ixSegm))); % Method 2, histogram with amplitude!
            end
            
            % Bin histogram - iChan x length of PeriWind*2
            v_Histbin = sum(mPeriEv,1); % zeros(1,nbins-1);
            v_HistbinAmpl = sum(mPeriEvAmpl,1);
            if plot_fig
            % plot histigram
            if (iChan1 == ch2plot(2)&& iChan == ch2plot(1))
            subplot(5,1,5);
            bar(v_Histbin);
            end
            end
            % Normalize histogram
            v_Histbin = v_Histbin + eps; % Avoid zeros
            v_Histbin = v_Histbin/sum(v_Histbin);
            v_HistbinAmpl = v_HistbinAmpl + eps;
            v_HistbinAmpl = v_HistbinAmpl/sum(v_HistbinAmpl);
            
            % Shannon entropy
            s_ShannonE = (s_ShannonEMax + sum(v_Histbin.*log(v_Histbin),2))/s_ShannonEMax;
            mFConn(iChan1,iChan) = s_ShannonE;
            mFConn(iChan,iChan1) = s_ShannonE;
            
        end
    end
    
    FConnPairs(iWind,:) = mFConn(mTriu);
    FConnMatrix(iWind,:,:) = mFConn;
    
    FConnAmplPairs(iWind,:) = mFConnAmpl(mTriu);
    FConnAmplMatrix(iWind,:,:) = mFConnAmpl;
    
end

gec = FConnMatrix;
