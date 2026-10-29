function ArrayFailureStudioN
%ARRAYFAILURESTUDION  Failure analysis and beam-preserving reconfiguration of a 1 x N patch array.
%
%   >> ArrayFailureStudioN
%
%   N is set in the Settings panel (default 20, range 2-32).
%
%   ALGORITHM
%     Failed set F, healthy set S (|S| >= 2), faulted weights w_F.
%     Active elements re-solve weights w_S by constrained weighted least squares
%        min_w || W^(1/2) (G_S w_S - b) ||^2
%        s.t.  peak value of ideal beam kept (hard), dF/dtheta = dF/dphi = 0 at peak (soft)
%     Solved through the KKT system.
%     Reconfiguration: exhaustive search over healthy subsets when there are few healthy
%     patches (Auto: <= 8), otherwise greedy backward elimination.
%     Failure analysis: all failure combinations when there are few, otherwise a
%     reproducible random sample.
%
%   FEATURES
%     - Polar comparison: beam plane and orthogonal plane, ideal / failed / compensated.
%     - theta-phi maps of the three states.
%     - S-parameters: parametric resonator model, or a Touchstone file (.sNp from CST)
%       that includes mutual coupling.
%
%   Requires MATLAB R2020a or later. No toolboxes.

% =========================== constants ====================================
C0         = 299792458;                 % speed of light [m/s]
DBF        = -40;                       % dB floor for pattern drawing
TOL        = 0.5;                       % allowed pointing drift [deg]
NFR        = 14;                        % animation frames per transition
EXH_MAX_H  = 8;                         % Auto: exhaustive search up to this many healthy patches
EXH_HARD_H = 12;                        % absolute limit when exhaustive is forced
NSP        = 201;                       % S-parameter frequency points (model mode)
NROWS      = 40;                        % rows in the settings grid

[PHg,THg] = meshgrid(deg2rad(0:2:358), deg2rad(0:1:90));
TH   = THg(:);
PH   = PHg(:);
AREA = sin(TH);                         % solid-angle weight
dTh  = deg2rad(1);
dPh  = deg2rad(2);
TCUT = deg2rad(0:0.5:90);
PCUT = deg2rad(0:15:345);
PHV  = 0:2:358;                         % map axes [deg]
THV  = (0:1:90).';

pal = struct('navy',[0.11 0.23 0.42],'ink',[0.08 0.08 0.10],'paper',[1 1 1], ...
    'panel',[0.96 0.97 0.98],'copper',[0.85 0.55 0.22],'fail',[0.85 0.16 0.12], ...
    'comp',[0.05 0.42 0.85],'off',[0.72 0.72 0.74],'ground',[0.80 0.80 0.82], ...
    'sub',[0.55 0.78 0.55],'amber',[0.98 0.66 0.08],'green',[0.10 0.60 0.30]);

% =============================== state ====================================
N        = 20;
mdl      = [];
EL       = struct('hasE',false,'hasH',false,'E',[],'H',[]);
cur      = [];                          % current displayed state
proof    = [];                          % beam-proof data
caseTbl  = {};
caseMask = false(0,N);
srcMode  = 'model';                     % 'model' or 'file'
Sfile    = [];                          % Touchstone S(f,i,j)
fFile    = [];                          % Touchstone frequency [GHz]
lastMethod = 'n/a';
busy     = false;
ro       = 0;                           % settings-panel row counter

% ============================== figure ====================================
scr = get(0,'ScreenSize');
fig = uifigure('Name','ArrayFailureStudioN | 1 x N Patch Array Failure Analysis', ...
    'Position',[20 40 min(1680,scr(3)-40) min(980,scr(4)-90)], ...
    'Color',pal.paper,'Resize','on');
try
    fig.Theme = 'light';
catch
end

root = uigridlayout(fig,[2 2],'RowHeight',{74,'1x'},'ColumnWidth',{390,'1x'}, ...
    'Padding',[10 10 10 10],'RowSpacing',8,'ColumnSpacing',10, ...
    'BackgroundColor',pal.paper);

% ---- header ----
hp = uipanel(root,'BackgroundColor',pal.navy,'BorderType','none');
place(hp,1,[1 2]);
hg = uigridlayout(hp,[2 2],'RowHeight',{'1x','1x'},'ColumnWidth',{'1x',560}, ...
    'BackgroundColor',pal.navy,'Padding',[18 4 18 4],'RowSpacing',0);
place(uilabel(hg,'Text','Array Element Reduction / Failure Analysis (1 x N)', ...
    'FontSize',22,'FontWeight','bold','FontColor',[1 1 1]),1,1);
place(uilabel(hg,'Text', ...
    'constrained least-squares reconfiguration  |  polar comparison  |  S-parameters', ...
    'FontSize',11,'FontColor',[0.82 0.90 1]),2,1);
lblStage = place(uilabel(hg,'Text','STAGE: idle','FontSize',14,'FontWeight','bold', ...
    'FontColor',pal.amber,'HorizontalAlignment','right'),1,2);
place(uilabel(hg,'Text','Rotate: drag in 3D views   Zoom: scroll wheel', ...
    'FontSize',11,'FontColor',[0.82 0.90 1],'HorizontalAlignment','right'),2,2);

% ---- left column: actions always on top, settings below (scrollable) ----
leftG = uigridlayout(root,[2 1],'RowHeight',{300,'1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[0 0 0 0],'RowSpacing',6);
place(leftG,2,1);

actP = uipanel(leftG,'Title','Actions','FontWeight','bold', ...
    'BackgroundColor',pal.paper,'ForegroundColor',pal.navy);
place(actP,1,1);

pnl = uipanel(leftG,'Title','Settings (scroll for more)','FontWeight','bold', ...
    'BackgroundColor',pal.paper,'ForegroundColor',pal.navy);
place(pnl,2,1);

cg = uigridlayout(pnl,[NROWS 2],'BackgroundColor',pal.paper,'Padding',[8 6 8 6], ...
    'RowSpacing',3,'ColumnWidth',{'1x','1x'});
rh = repmat({24},1,NROWS);
rh{NROWS} = '1x';
try
    cg.Scrollable = 'on';
catch
end

% ---- settings ----
sectionRow('ARRAY & SUBSTRATE');
edN    = numRow('Number of patches N',N,[2 32]);
edN.ValueChangedFcn = @(~,~)onNChanged();
edF    = numRow('Frequency f (GHz)',10,[1 100]);
edEr   = numRow('Substrate / surrounding eps_r',4.4,[1 20]);
edH    = numRow('Substrate height h (mm)',1.6,[0.1 10]);
edD    = numRow('Spacing d (lambda0)',0.5,[0.2 2]);
edTh0  = numRow('Steer theta0 (deg)',30,[0 90]);
edPh0  = numRow('Steer phi0 (deg)',0,[-360 360]);
ddTap  = ddRow('Amplitude taper',{'Uniform','Hamming'},'Uniform');
edMainW = numRow('Main-beam radius (deg)',15,[2 90]);

sectionRow('ELEMENT PATTERN');
ddEl   = ddRow('Element source',{'cos^q model','CST cuts'},'cos^q model');
edQ    = numRow('cos^q exponent q',1.5,[0 10]);
btnRow('Load CST E-plane cut (phi = 0)',@(~,~)loadCut('E'),pal.paper,pal.ink);
btnRow('Load CST H-plane cut (phi = 90)',@(~,~)loadCut('H'),pal.paper,pal.ink);
btnRow('Clear CST cuts',@(~,~)clearCuts,pal.paper,pal.ink);

sectionRow('S-PARAMETERS');
ddSrc  = ddRow('Source',{'Parametric model','Touchstone file'},'Parametric model');
ddSrc.ValueChangedFcn = @(~,~)onSrcChanged();
btnRow('Load Touchstone (.sNp) from CST',@(~,~)loadTouchstone,pal.paper,pal.ink);
edRin  = numRow('Feed input resistance (ohm)',50,[1 500]);
edQp   = numRow('Patch Q (bandwidth ~ 1/Q)',30,[2 200]);
edDet  = numRow('Resonance offset (%)',0,[-20 20]);
edSpan = numRow('Sweep span (GHz)',2,[0.1 20]);

sectionRow('FAILURE & ALGORITHM');
ro = ro + 1;
place(uilabel(cg,'Text','Failed patches (ctrl / shift click for several)','FontColor',pal.ink),ro,[1 2]);
ro = ro + 1;
lbFail = place(uilistbox(cg,'Items',patchItems(N),'Multiselect','on','Value',{}, ...
    'FontColor',pal.ink,'BackgroundColor',pal.paper),ro,[1 2]);
rh{ro} = 130;
sld    = sldRow('Failure severity (0 = healthy, 1 = dead)',1);
edMain = numRow('Main-beam emphasis (x)',10,[1 1000]);
edLock = numRow('Beam-lock weight',1000,[0 1e12]);
edEta  = numRow('Radiation efficiency (%)',90,[1 100]);
edSpd  = numRow('Animation speed (ms/frame)',30,[5 500]);
ddSearch = ddRow('Reconfiguration search',{'Auto','Exhaustive','Greedy'},'Auto');

sectionRow('FAILURE ANALYSIS');
edMaxF = numRow('Max failed patches per case',3,[1 10]);
edK    = numRow('Max cases (random sample if more)',30,[1 500]);
edSeed = numRow('Random seed',1,[0 1e6]);

cg.RowHeight = rh;

% ---- actions (always visible) ----
ag = uigridlayout(actP,[7 1],'RowHeight',{30,30,30,30,30,30,'1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[8 6 8 6],'RowSpacing',4);
actBtn(1,'RUN LIVE PROCESS (animated)',@(~,~)runProcess,pal.navy,[1 1 1]);
actBtn(2,'FULL FAILURE ANALYSIS (sampled cases)',@(~,~)runAnalysis,pal.navy,[1 1 1]);
actBtn(3,'Compute current case',@(~,~)computeCurrent,pal.panel,pal.ink);
actBtn(4,'Greedy vs exhaustive (current case)',@(~,~)compareGreedy,pal.panel,pal.ink);
actBtn(5,'Beam proof (all phi cuts)',@(~,~)beamProofButton,pal.panel,pal.ink);
actBtn(6,'Export results',@(~,~)exportCB,pal.panel,pal.ink);
lblStatus = uilabel(ag,'Text','Ready.','WordWrap','on','FontColor',pal.ink, ...
    'VerticalAlignment','top');
lblStatus.Layout.Row    = 7;
lblStatus.Layout.Column = 1;

% ---- tabs ----
tg = uitabgroup(root);
place(tg,2,2);
tabLive  = uitab(tg,'Title','Live 3D');
tabBA    = uitab(tg,'Title','Before / After');
tabProof = uitab(tg,'Title','Beam Proof');
tabIns   = uitab(tg,'Title','Angle Inspector');
tabPol   = uitab(tg,'Title','Polar & Maps');
tabS     = uitab(tg,'Title','S-Parameters');
tabAn    = uitab(tg,'Title','Failure Analysis');
tabRep   = uitab(tg,'Title','Report');

% ---- Live 3D tab ----
gL = uigridlayout(tabLive,[2 3],'RowHeight',{'1.05x','1x'},'ColumnWidth',{'1x','1x','1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',6);
axG = uiaxes(gL);  place(axG,1,[1 2]);
txtLog = uitextarea(gL,'Editable','off','FontName','Menlo','FontSize',11, ...
    'BackgroundColor',pal.panel,'FontColor',pal.ink,'Value',{'Ready. Press RUN LIVE PROCESS.'});
place(txtLog,1,3);
axLI = uiaxes(gL);  place(axLI,2,1);
axLF = uiaxes(gL);  place(axLF,2,2);
axLC = uiaxes(gL);  place(axLC,2,3);

% ---- Before / After tab ----
gB = uigridlayout(tabBA,[2 3],'RowHeight',{'1x','1x'},'ColumnWidth',{'1x','1x','1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',6);
axBI  = uiaxes(gB);  place(axBI,1,1);
axBIn = uiaxes(gB);  place(axBIn,1,2);
axBF  = uiaxes(gB);  place(axBF,1,3);
axBC  = uiaxes(gB);  place(axBC,2,1);
axBGD = uiaxes(gB);  place(axBGD,2,2);
tblKPI = uitable(gB,'ColumnName',{'Metric','Ideal','Failed (before)','Compensated (after)'}, ...
    'FontSize',11,'ForegroundColor',pal.ink);
place(tblKPI,2,3);

% ---- Beam proof tab ----
gP = uigridlayout(tabProof,[2 2],'RowHeight',{'1.1x','1x'},'ColumnWidth',{'1x','1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',6);
axPk = uiaxes(gP);  place(axPk,1,1);
axDv = uiaxes(gP);  place(axDv,1,2);
tblProof = uitable(gP,'ColumnName', ...
    {'phi (deg)','theta ideal','theta failed','theta compensated','dev failed (deg)','dev comp (deg)','status'}, ...
    'FontSize',10,'ForegroundColor',pal.ink);
place(tblProof,2,[1 2]);

% ---- Angle inspector tab ----
gI = uigridlayout(tabIns,[2 2],'RowHeight',{200,'1x'},'ColumnWidth',{'1x','1.5x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',8);
gIc = uigridlayout(gI,[5 1],'RowHeight',{22,30,22,30,30}, ...
    'BackgroundColor',pal.paper,'RowSpacing',4);
place(uilabel(gIc,'Text','Polar angle theta from vertical (deg)','FontColor',pal.ink),1,1);
sI = uislider(gIc,'Limits',[0 90],'Value',30);      place(sI,2,1);
place(uilabel(gIc,'Text','Azimuth phi (deg)','FontColor',pal.ink),3,1);
sP = uislider(gIc,'Limits',[0 360],'Value',0);      place(sP,4,1);
lblIns = place(uilabel(gIc,'Text','theta = 30.0 deg   phi = 0.0 deg', ...
    'FontWeight','bold','FontColor',pal.navy),5,1);
sI.ValueChangedFcn = @(s,~)updateInspector(s.Value,sP.Value);
sP.ValueChangedFcn = @(s,~)updateInspector(sI.Value,s.Value);
try
    sI.ValueChangingFcn = @(~,e)updateInspector(e.Value,sP.Value);
    sP.ValueChangingFcn = @(~,e)updateInspector(sI.Value,e.Value);
catch
end
tblIns = uitable(gI,'ColumnName',{'State','|E| dB re ideal','Local gain dBi', ...
    'Phase (deg)','Delta vs ideal (dB)'},'FontSize',11,'ForegroundColor',pal.ink);
place(tblIns,1,2);
axInsCut = uiaxes(gI);  place(axInsCut,2,[1 2]);

% ---- Polar & maps tab ----
gPo = uigridlayout(tabPol,[2 3],'RowHeight',{'1x','1x'},'ColumnWidth',{'1x','1x','0.8x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',6);
axPA = uiaxes(gPo);  place(axPA,1,1);
axPB = uiaxes(gPo);  place(axPB,1,2);
txtPo = uitextarea(gPo,'Editable','off','FontName','Menlo','FontSize',11, ...
    'BackgroundColor',pal.panel,'FontColor',pal.ink,'Value',{'Run a case to fill this panel.'});
place(txtPo,1,3);
axMI = uiaxes(gPo);  place(axMI,2,1);
axMF = uiaxes(gPo);  place(axMF,2,2);
axMC = uiaxes(gPo);  place(axMC,2,3);

% ---- S-parameter tab ----
gS = uigridlayout(tabS,[3 2],'RowHeight',{32,'1x','1x'},'ColumnWidth',{'1x','1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',6);
tb = uigridlayout(gS,[1 6],'ColumnWidth',{50,80,80,80,110,'1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[0 0 0 0]);
place(tb,1,[1 2]);
place(uilabel(tb,'Text','Port','FontColor',pal.ink),1,1);
ddSP = place(uidropdown(tb,'Items',patchItems(N),'Value','P1', ...
    'FontColor',pal.ink,'BackgroundColor',pal.paper),1,2);
place(uilabel(tb,'Text','Pair with','FontColor',pal.ink),1,3);
ddSQ = place(uidropdown(tb,'Items',patchItems(N),'Value','P2', ...
    'FontColor',pal.ink,'BackgroundColor',pal.paper),1,4);
place(uilabel(tb,'Text','Heatmap state','FontColor',pal.ink),1,5);
ddSS = place(uidropdown(tb,'Items',{'Compensated','Failed','Ideal'},'Value','Compensated', ...
    'FontColor',pal.ink,'BackgroundColor',pal.paper),1,6);
ddSP.ValueChangedFcn = @(~,~)updateSPlots();
ddSQ.ValueChangedFcn = @(~,~)updateSPlots();
ddSS.ValueChangedFcn = @(~,~)updateSPlots();
axS1 = uiaxes(gS);  place(axS1,2,1);
axS2 = uiaxes(gS);  place(axS2,2,2);
axS3 = uiaxes(gS);  place(axS3,3,1);
axS4 = uiaxes(gS);  place(axS4,3,2);

% ---- Failure analysis tab ----
gA = uigridlayout(tabAn,[2 2],'RowHeight',{'1.3x','1x'},'ColumnWidth',{'1.8x','1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6],'RowSpacing',6,'ColumnSpacing',8);
tblA = uitable(gA,'ColumnName',{'#','Failed patches','Healthy','Method','Active set', ...
    'Fidelity %','Gain dBi','Pointing err (deg)','Max sidelobe dB','Greedy fidelity %', ...
    'Verdict','No-fix fidelity %'},'FontSize',10,'ForegroundColor',pal.ink);
place(tblA,1,[1 2]);
tblA.CellSelectionCallback = @(~,evt)caseRowClicked(evt);
axBar = uiaxes(gA);  place(axBar,2,1);
txtBest = uitextarea(gA,'Editable','off','FontName','Menlo','FontSize',11, ...
    'BackgroundColor',pal.panel,'FontColor',pal.ink, ...
    'Value',{'Press FULL FAILURE ANALYSIS, or click a row after running it.'});
place(txtBest,2,2);

% ---- report tab ----
gR = uigridlayout(tabRep,[1 1],'RowHeight',{'1x'},'ColumnWidth',{'1x'}, ...
    'BackgroundColor',pal.paper,'Padding',[6 6 6 6]);
txtRep = uitextarea(gR,'Editable','off','FontName','Menlo','FontSize',11, ...
    'BackgroundColor',pal.panel,'FontColor',pal.ink,'Value',{'No report yet.'});

% ---- style axes ----
for ax = [axG axLI axLF axLC axBI axBF axBC]
    fix3d(ax);
end
styleAx(axBIn,'Excitation amplitude |w_n|');
styleAx(axBGD,'Gain and directivity');
styleAx(axPk,'Beam peak theta vs azimuth cut');
styleAx(axDv,'Beam-angle deviation vs ideal');
styleAx(axInsCut,'E-plane cut at selected phi (dB re ideal peak)');
styleAx(axBar,'Fidelity per failure case (higher = closer to ideal)');
styleAx(axMI,'Ideal');
styleAx(axMF,'Failed (before)');
styleAx(axMC,'Compensated (after)');

% ---- 3D pattern panels (each state appears in a live and a before/after view) ----
SI = {makeSurf(axLI,'Ideal (healthy)'), makeSurf(axBI,'Ideal (healthy)')};
SF = {makeSurf(axLF,'Failed, before algorithm'), makeSurf(axBF,'Failed, before algorithm')};
SC = {makeSurf(axLC,'Compensated, after algorithm'), makeSurf(axBC,'Compensated, after algorithm')};

% ---- startup ----
N = round(edN.Value);
prepareModel();
showAll(mdl.wI,mdl.wI,false(N,1),true(N,1),'ideal');
setStatus('Ready. Choose failed patches in the list, then press an action.');
setStage('idle');

% =========================================================================
%  NESTED: UI helpers
% =========================================================================

    function h = place(h,r,c)
        h.Layout.Row    = r;
        h.Layout.Column = c;
    end

    function sectionRow(txt)
        ro = ro + 1;
        place(uilabel(cg,'Text',txt,'FontWeight','bold','FontColor',pal.navy,'FontSize',11), ro, [1 2]);
    end

    function h = numRow(txt,val,lims)
        ro = ro + 1;
        place(uilabel(cg,'Text',txt,'FontColor',pal.ink),ro,1);
        h = place(uieditfield(cg,'numeric','Value',val,'Limits',lims, ...
            'FontColor',pal.ink,'BackgroundColor',pal.paper),ro,2);
    end

    function h = ddRow(txt,items,val)
        ro = ro + 1;
        place(uilabel(cg,'Text',txt,'FontColor',pal.ink),ro,1);
        h = place(uidropdown(cg,'Items',items,'Value',val, ...
            'FontColor',pal.ink,'BackgroundColor',pal.paper),ro,2);
    end

    function h = btnRow(txt,cb,bg,fc)
        ro = ro + 1;
        h = place(uibutton(cg,'Text',txt,'ButtonPushedFcn',cb, ...
            'BackgroundColor',bg,'FontColor',fc,'FontWeight','bold'),ro,[1 2]);
    end

    function h = sldRow(txt,val)
        ro = ro + 1;
        place(uilabel(cg,'Text',txt,'FontColor',pal.ink),ro,[1 2]);
        ro = ro + 1;
        h = place(uislider(cg,'Limits',[0 1],'Value',val),ro,[1 2]);
        rh{ro} = 30;
    end

    function h = actBtn(r,txt,cb,bg,fc)
        h = uibutton(ag,'Text',txt,'ButtonPushedFcn',cb, ...
            'BackgroundColor',bg,'FontColor',fc,'FontWeight','bold');
        h.Layout.Row    = r;
        h.Layout.Column = 1;
    end

    function setStatus(msg)
        lblStatus.Text = msg;
        drawnow limitrate;
    end

    function setStage(msg)
        lblStage.Text = ['STAGE: ' msg];
        drawnow limitrate;
    end

    function logMsg(varargin)
        msg = sprintf(varargin{:});
        txtLog.Value = [txtLog.Value(:); {msg}];
        drawnow limitrate;
    end

    function animPause(t)
        drawnow;
        pause(t);
    end

    function setBusy(v)
        busy = v;
    end

    function fm = getFail()
        fm = false(N,1);
        v = lbFail.Value;
        if ischar(v) || isstring(v)
            v = cellstr(v);
        end
        for k = 1:N
            fm(k) = any(strcmp(v,sprintf('P%d',k)));
        end
    end

    function setFailSelection(fm)
        items = patchItems(N);
        sel = items(logical(fm(:).'));
        lbFail.Value = sel;
    end

    function s = elemName()
        if strcmp(ddEl.Value,'CST cuts') && (EL.hasE || EL.hasH)
            parts = {};
            if EL.hasE, parts{end+1} = EL.E.file; end %#ok<AGROW>
            if EL.hasH, parts{end+1} = EL.H.file; end %#ok<AGROW>
            s = ['CST cuts: ' strjoin(parts,' + ')];
        else
            s = sprintf('cos^%.2g model', edQ.Value);
        end
    end

    function h = makeSurf(ax,ttl)
        Z0 = zeros(size(THg));
        h.s = surf(ax,Z0,Z0,Z0,Z0,'EdgeColor','none');
        hold(ax,'on');
        h.b = plot3(ax,[0 0],[0 0],[0 1.05],'-','Color',pal.amber,'LineWidth',2.5);
        h.m = plot3(ax,0,0,1,'p','MarkerSize',13,'MarkerFaceColor',pal.amber, ...
            'MarkerEdgeColor',pal.ink);
        h.q = plot3(ax,0,0,0,'o','MarkerSize',9,'MarkerFaceColor',pal.navy, ...
            'MarkerEdgeColor','w','Visible','off');
        hold(ax,'off');
        colormap(ax,jet(256));
        caxis(ax,[DBF 0]);
        axis(ax,'equal');
        xlim(ax,[-1.15 1.15]);  ylim(ax,[-1.15 1.15]);  zlim(ax,[0 1.2]);
        view(ax,[35 25]);
        styleAx(ax,ttl);
        xlabel(ax,'x','Color',pal.ink);
        ylabel(ax,'y','Color',pal.ink);
        zlabel(ax,'z','Color',pal.ink);
        cb = colorbar(ax);
        cb.Label.String = 'dB (0 at edge, -40 at centre)';
    end

    function setSurf(h,w)
        F  = mdl.G*w;
        A  = abs(reshape(F,size(THg)))/mdl.refV;
        Ec = max(20*log10(max(A,1e-6)),DBF);
        R  = max(0,(Ec-DBF)/(-DBF));
        set(h.s,'XData',R.*sin(THg).*cos(PHg),'YData',R.*sin(THg).*sin(PHg), ...
            'ZData',R.*cos(THg),'CData',Ec);
        [pT,pP] = peakOf(w);
        px = sin(pT)*cos(pP);  py = sin(pT)*sin(pP);  pz = cos(pT);
        set(h.b,'XData',[0 1.05*px],'YData',[0 1.05*py],'ZData',[0 1.05*pz]);
        set(h.m,'XData',px,'YData',py,'ZData',pz);
    end

    function setPanels(list,w)
        for k = 1:numel(list)
            setSurf(list{k},w);
        end
    end

    function styleAx(ax,ttl)
        ax.Color     = pal.paper;
        ax.XColor    = pal.ink;
        ax.YColor    = pal.ink;
        ax.ZColor    = pal.ink;
        ax.GridColor = [0.75 0.75 0.78];
        ax.FontSize  = 10;
        box(ax,'on');
        grid(ax,'on');
        title(ax,ttl,'Color',pal.ink,'FontWeight','bold');
    end

% =========================================================================
%  NESTED: element model
% =========================================================================

    function F = elemField(T,P)
        if strcmp(ddEl.Value,'CST cuts') && (EL.hasE || EL.hasH)
            if EL.hasE && EL.hasH
                wE = cos(P).^2;
                F  = wE.*cutVal(EL.E,T) + (1-wE).*cutVal(EL.H,T);
            elseif EL.hasE
                F  = cutVal(EL.E,T);
            else
                F  = cutVal(EL.H,T);
            end
        else
            F = max(cos(T),0).^edQ.Value;
        end
    end

    function F = cutVal(C,T)
        td = rad2deg(T);
        F  = C.fr(td) + 1j*C.fi(td);
    end

    function g = rowAt(T,P,x)
        % n x N matrix: element factor times each element's phase term
        T = T(:);  P = P(:);
        u = sin(T).*cos(P);
        g = exp(1j*2*pi*(u*x.')) .* elemField(T,P);
    end

    function loadCut(which)
        [f,p] = uigetfile({'*.txt;*.csv;*.dat','CST 1D export (*.txt, *.csv, *.dat)'}, ...
            sprintf('Select CST %s-plane 1D far-field file',which));
        if isequal(f,0), return; end
        try
            C = parseCSTCut(fullfile(p,f));
        catch err
            uialert(fig,err.message,'Could not read CST file','Icon','error');
            setStatus('CST import failed. See the message box.');
            return
        end
        C.file = f;
        if strcmp(which,'E')
            EL.E = C;  EL.hasE = true;
        else
            EL.H = C;  EL.hasH = true;
        end
        ddEl.Value = 'CST cuts';
        setStatus(sprintf('Loaded %s-plane cut %s (%d points). Press Compute current case.', ...
            which,f,C.npts));
    end

    function clearCuts()
        EL = struct('hasE',false,'hasH',false,'E',[],'H',[]);
        ddEl.Value = 'cos^q model';
        setStatus('CST cuts cleared. Using the cos^q model.');
    end

% =========================================================================
%  NESTED: physical model
% =========================================================================

    function prepareModel()
        f  = edF.Value*1e9;
        er = edEr.Value;
        h  = edH.Value*1e-3;
        Wm = C0/(2*f)*sqrt(2/(er+1));
        ee = (er+1)/2 + (er-1)/2*(1+12*h/Wm)^(-0.5);
        dL = 0.412*h*(ee+0.3)*(Wm/h+0.264)/((ee-0.258)*(Wm/h+0.8));
        Lm = C0/(2*f*sqrt(ee)) - 2*dL;
        lam0 = C0/f*1e3;

        N   = round(edN.Value);
        th0 = deg2rad(edTh0.Value);
        ph0 = deg2rad(edPh0.Value);
        d   = edD.Value;
        x   = ((1:N).' - (N+1)/2)*d;
        us  = sin(th0)*cos(ph0);
        if strcmp(ddTap.Value,'Hamming')
            tp = 0.54 - 0.46*cos(2*pi*(0:N-1).'/max(N-1,1));
        else
            tp = ones(N,1);
        end
        wI = tp.*exp(-1j*2*pi*x*us);

        M = struct('N',N,'x',x,'d',d,'th0',th0,'ph0',ph0,'wI',wI);
        M.lam0 = lam0;
        M.L    = Lm*1e3;
        M.W    = Wm*1e3;
        M.h    = h*1e3;
        M.ee   = ee;
        M.xm   = x*lam0;
        M.dmm  = d*lam0;
        M.mainWdeg = edMainW.Value;

        M.G  = rowAt(TH,PH,x);
        M.Fi = M.G*wI;
        [~,k] = max(abs(M.Fi));
        [pT,pP,pV] = refinePeak(x,wI,TH(k),PH(k));
        M.pT   = pT;
        M.pP   = pP;
        M.refV = pV;

        main  = angDist(TH,PH,pT,pP) <= M.mainWdeg;
        M.ws  = AREA.*(1+(edMain.Value-1)*main);

        hh = 1e-4;
        M.gP     = rowAt(pT,pP,x);
        M.gPt    = (rowAt(pT+hh,pP,x) - rowAt(pT-hh,pP,x))/(2*hh);
        M.gPp    = (rowAt(pT,pP+hh,x) - rowAt(pT,pP-hh,x))/(2*hh);
        M.target = M.gP*wI;
        mdl = M;
    end

    function [pT,pP,pV] = refinePeak(x,w,t0,p0)
        tt = linspace(max(t0-deg2rad(2),0), min(t0+deg2rad(2),pi/2), 41);
        pp = linspace(p0-deg2rad(2), p0+deg2rad(2), 41);
        [PP,TT] = meshgrid(pp,tt);
        V = abs(rowAt(TT(:),PP(:),x)*w);
        [pV,kk] = max(V);
        pT = TT(kk);
        pP = PP(kk);
    end

    function [pT,pP,pV] = peakOf(w)
        [~,k] = max(abs(mdl.G*w));
        [pT,pP,pV] = refinePeak(mdl.x,w,TH(k),PH(k));
    end

    function wF = failedWeights(fm,sev)
        wF = mdl.wI;
        wF(fm) = mdl.wI(fm)*(1-sev);
    end

% =========================================================================
%  NESTED: algorithm
% =========================================================================

    function [wC,info] = solveSubset(fm,act,sev)
        % Constrained weighted least squares, solved via the KKT system.
        wF = failedWeights(fm,sev);
        wF(~act & ~fm) = 0;
        F  = find(fm);
        S  = find(act);
        k  = numel(S);

        b  = mdl.Fi - mdl.G(:,F)*wF(F);
        sw = sqrt(mdl.ws);
        M  = sw.*mdl.G(:,S);
        bb = sw.*b;

        gFt = mdl.gPt(1,F)*wF(F);
        gFp = mdl.gPp(1,F)*wF(F);
        Cr  = mdl.gP(1,S);
        dC  = mdl.target - mdl.gP(1,F)*wF(F);

        lockW    = edLock.Value;
        bestErr  = inf;
        bestW    = wF;
        bestPass = 0;
        bestLock = lockW;
        bestRc   = NaN;
        for pass = 1:5
            lw  = sqrt(lockW);
            Dm  = lw*[mdl.gPt(1,S); mdl.gPp(1,S)];
            Db  = -lw*[gFt; gFp];
            A   = [M; Dm];
            y   = [bb; Db];
            R   = A'*A;
            reg = 1e-9*max(real(diag(R))) + eps;
            KK  = [R + reg*eye(k), Cr'; Cr, 0];
            sol = KK\[A'*y; dC];
            wS  = sol(1:k);

            wT = wF;
            wT(S) = wS;
            [pT,pP] = peakOf(wT);
            e = angDist(pT,pP,mdl.pT,mdl.pP);
            if e < bestErr
                bestErr  = e;
                bestW    = wT;
                bestPass = pass;
                bestLock = lockW;
                bestRc   = rcond(KK);
            end
            if e < 0.25, break; end
            lockW = lockW*10;
        end
        wC   = bestW;
        info = struct('err',bestErr,'pass',bestPass,'lock',bestLock,'rcond',bestRc,'nAct',k);
    end

    function m = metricsOf(w)
        F = mdl.G*w;
        [pT,pP,pV] = peakOf(w);
        m.gl    = 20*log10(max(pV,1e-12)/mdl.refV);
        m.pt    = rad2deg(pT);
        m.pp    = rad2deg(pP);
        m.ptErr = angDist(pT,pP,mdl.pT,mdl.pP);
        far = angDist(TH,PH,pT,pP) > mdl.mainWdeg;
        if any(far)
            m.sll = 20*log10(max(abs(F(far)))/mdl.refV);
        else
            m.sll = NaN;
        end
        m.prad = sum(AREA.*abs(F).^2)*dTh*dPh;
        m.dir  = 10*log10(4*pi*pV^2/m.prad);
        m.gain = m.dir + 10*log10(edEta.Value/100);
        m.corr = 100*abs(sum(AREA.*conj(mdl.Fi).*F))^2 / ...
                 (sum(AREA.*abs(mdl.Fi).^2)*sum(AREA.*abs(F).^2) + eps);
        m.rms  = 100*sqrt(sum(AREA.*abs(F-mdl.Fi).^2)/sum(AREA.*abs(mdl.Fi).^2));
        m.pwr  = 10*log10(sum(abs(w).^2)/sum(abs(mdl.wI).^2));
    end

    function s = scoreOf(m)
        s = m.corr + 1000*double(m.ptErr <= TOL);
    end

    function R = emptyR()
        R = struct('ok',false,'wC',[],'act',false(N,1),'m',[],'info',[], ...
            'n',0,'method','none','label','n/a','score',-inf);
    end

    function R = packR(wC,act,m,info,n,method)
        R = struct('ok',true,'wC',wC,'act',act,'m',m,'info',info, ...
            'n',n,'method',method,'label',patchList(act,', ','-'),'score',scoreOf(m));
    end

    function R = reconfigure(fm,sev,anim)
        % Chooses the search: exhaustive for few healthy patches, greedy otherwise.
        nH = sum(~fm);
        if nH < 2
            R = emptyR();
            return
        end
        switch ddSearch.Value
            case 'Exhaustive'
                useExh = nH <= EXH_HARD_H;
                if ~useExh
                    setStatus(sprintf('Exhaustive search limited to %d healthy patches; using greedy.', EXH_HARD_H));
                end
            case 'Greedy'
                useExh = false;
            otherwise
                useExh = nH <= EXH_MAX_H;
        end
        if useExh
            R = exhaustiveSearch(fm,sev,anim);
        else
            R = greedySearch(fm,sev,anim);
        end
    end

    function R = exhaustiveSearch(fm,sev,anim)
        % Every healthy subset with at least two elements.
        subs = subsetsOf(find(~fm).');
        R = emptyR();
        for s = 1:numel(subs)
            act = false(N,1);
            act(subs{s}) = true;
            [wC,info] = solveSubset(fm,act,sev);
            m  = metricsOf(wC);
            sc = scoreOf(m);
            if anim
                drawGeom(wC,fm,act,'search');
                logMsg('   subset {%s} fidelity %6.2f %%  pointing %5.3f deg  %s', ...
                    patchList(act,',','-'), m.corr, m.ptErr, ternStr(m.ptErr <= TOL,'holds beam','drifts'));
                animPause(edSpd.Value/1000*2);
            end
            if sc > R.score
                R = packR(wC,act,m,info,numel(subs),'Exhaustive');
            end
        end
        R.n = numel(subs);
    end

    function R = greedySearch(fm,sev,anim)
        % Greedy backward elimination: drop the element whose removal hurts least.
        S = find(~fm).';
        nEval = 0;
        R = emptyR();
        if numel(S) < 2, return; end
        act = false(N,1);  act(S) = true;
        [wC,info] = solveSubset(fm,act,sev);
        m = metricsOf(wC);  nEval = nEval + 1;
        R = packR(wC,act,m,info,nEval,'Greedy');
        while numel(S) > 2
            bestSc = -inf;  bestS = [];  bestA = [];  bestW = [];  bestM = [];  bestI = [];  drop = 0;
            for s = S
                Sn = S(S ~= s);
                actn = false(N,1);  actn(Sn) = true;
                [wn,in] = solveSubset(fm,actn,sev);
                mn = metricsOf(wn);  nEval = nEval + 1;
                if scoreOf(mn) > bestSc
                    bestSc = scoreOf(mn);  bestS = Sn;  bestA = actn;
                    bestW = wn;  bestM = mn;  bestI = in;  drop = s;
                end
            end
            S = bestS;
            if anim
                drawGeom(bestW,fm,bestA,'search');
                logMsg('   drop P%d -> {%s} fidelity %6.2f %%  pointing %5.3f deg  %s', ...
                    drop, patchList(bestA,',','-'), bestM.corr, bestM.ptErr, ...
                    ternStr(bestM.ptErr <= TOL,'holds beam','drifts'));
                animPause(edSpd.Value/1000*2);
            end
            if bestSc > R.score
                R = packR(bestW,bestA,bestM,bestI,nEval,'Greedy');
            end
        end
        R.n = nEval;
    end

    function [pT,pP] = cutPeakTheta(w,ph)
        % ph in radians; returns peak theta and phi in degrees
        tt = deg2rad(0:0.25:90);
        E  = abs(rowAt(tt(:),ph*ones(numel(tt),1),mdl.x)*w);
        [~,k] = max(E);
        t = tt(k);
        if k > 1 && k < numel(tt)
            a = E(k-1);  b = E(k);  c = E(k+1);
            den = a - 2*b + c;
            if den ~= 0
                t = tt(k) + 0.5*(a-c)/den*deg2rad(0.25);
            end
        end
        pT = rad2deg(t);
        pP = rad2deg(ph);
    end

    function [ang,dB] = cutPolar(w,phDeg)
        % Full azimuth plane: angle in (-90, 90) deg. Positive side is phi,
        % negative side is phi + 180 (theta mirrored). 0 deg = vertical.
        p1 = deg2rad(phDeg);
        n  = numel(TCUT);
        EA = abs(rowAt(TCUT(:),p1*ones(n,1),mdl.x)*w)/mdl.refV;
        EB = abs(rowAt(TCUT(:),(p1+pi)*ones(n,1),mdl.x)*w)/mdl.refV;
        dA  = max(20*log10(max(EA,1e-6)),DBF);
        dBv = max(20*log10(max(EB,1e-6)),DBF);
        th = rad2deg(TCUT(:));
        ang = [-flipud(th(2:end)); th];
        dB  = [flipud(dBv(2:end)); dA];
    end

% =========================================================================
%  NESTED: display
% =========================================================================

    function drawGeom(w,fm,act,mode)
        cla(axG);  hold(axG,'on');
        xm = mdl.xm;  L = mdl.L;  W = mdl.W;  h = mdl.h;
        xl = [min(xm)-L, max(xm)+L];
        yl = 1.8*W;
        labelAll = N <= 12;
        patch(axG,[xl(1) xl(2) xl(2) xl(1)],[-yl -yl yl yl],[0 0 0 0], ...
            'FaceColor',pal.ground,'EdgeColor','none','FaceAlpha',0.85);
        boxp(axG,xl(1),xl(2),-yl,yl,0,h,pal.sub,0.35);
        text(axG,xl(1),-yl,0,'  ground plane','FontSize',9,'Color',pal.ink);
        for n = 1:N
            if fm(n)
                c = pal.fail;
            elseif ~act(n)
                c = pal.off;
            elseif any(strcmp(mode,{'comp','search'}))
                c = pal.comp;
            else
                c = pal.copper;
            end
            X = xm(n) + [-1 1 1 -1]*L/2;
            Y = [-1 -1 1 1]*W/2;
            Z = h*[1 1 1 1];
            patch(axG,X,Y,Z,c,'EdgeColor',pal.ink,'LineWidth',1);
            plot3(axG,[xm(n) xm(n)],[-W/2-7 -W/2],[h h],'-','Color',pal.ink,'LineWidth',1.5);
            if labelAll
                text(axG,xm(n),W/2+2.5,h+0.8,sprintf('P%d',n),'HorizontalAlignment','center', ...
                    'FontWeight','bold','Color',pal.ink,'FontSize',9);
                text(axG,xm(n),-W/2-9,h+0.8,sprintf('|w| %.2f  %+.0f%s',abs(w(n)), ...
                    rad2deg(angle(w(n))),char(176)),'HorizontalAlignment','center', ...
                    'FontSize',8,'Color',pal.ink);
            end
            if fm(n)
                text(axG,xm(n),0,h+2.5,'FAILED','HorizontalAlignment','center', ...
                    'FontWeight','bold','Color',pal.fail,'FontSize',9);
            elseif ~act(n)
                text(axG,xm(n),0,h+2.5,'OFF','HorizontalAlignment','center', ...
                    'FontWeight','bold','Color',[0.4 0.4 0.4],'FontSize',9);
            end
        end
        hold(axG,'off');
        axis(axG,'equal');
        xlim(axG,xl+[-3 3]);
        ylim(axG,[-yl-12, yl+4]);
        zlim(axG,[0 max(h*7,5)]);
        view(axG,[-35 22]);
        title(axG,sprintf('Physical geometry (mm): N %d, L %.2f, W %.2f, h %.2f, spacing %.2f', ...
            N,L,W,h,mdl.dmm),'Color',pal.ink,'FontWeight','bold');
        xlabel(axG,'x (mm)','Color',pal.ink);
        ylabel(axG,'y (mm)','Color',pal.ink);
        zlabel(axG,'z (mm)','Color',pal.ink);
    end

    function showAll(wF,wC,fm,act,mode)
        cur = struct('wI',mdl.wI,'wF',wF,'wC',wC,'fm',fm,'act',act,'sev',sld.Value);
        setPanels(SI,mdl.wI);
        setPanels(SF,wF);
        setPanels(SC,wC);
        drawGeom(wC,fm,act,mode);

        mI = metricsOf(mdl.wI);
        mF = metricsOf(wF);
        mC = metricsOf(wC);
        tblKPI.Data = {
            'Gain (dBi)',            fx(mI.gain,'%.2f'), fx(mF.gain,'%.2f'), fx(mC.gain,'%.2f')
            'Directivity (dBi)',     fx(mI.dir,'%.2f'),  fx(mF.dir,'%.2f'),  fx(mC.dir,'%.2f')
            'Peak vs ideal (dB)',    fx(mI.gl,'%.2f'),   fx(mF.gl,'%.2f'),   fx(mC.gl,'%.2f')
            'Peak theta (deg)',      fx(mI.pt,'%.1f'),   fx(mF.pt,'%.1f'),   fx(mC.pt,'%.1f')
            'Peak phi (deg)',        fx(mI.pp,'%.1f'),   fx(mF.pp,'%.1f'),   fx(mC.pp,'%.1f')
            'Pointing error (deg)',  fx(mI.ptErr,'%.3f'),fx(mF.ptErr,'%.3f'),fx(mC.ptErr,'%.3f')
            'Max sidelobe (dB)',     fx(mI.sll,'%.1f'),  fx(mF.sll,'%.1f'),  fx(mC.sll,'%.1f')
            'Fidelity to ideal (%)', fx(mI.corr,'%.2f'), fx(mF.corr,'%.2f'), fx(mC.corr,'%.2f')
            'Pattern RMS error (%)', fx(mI.rms,'%.2f'),  fx(mF.rms,'%.2f'),  fx(mC.rms,'%.2f')
            'Excitation power (dB)', fx(mI.pwr,'%.2f'),  fx(mF.pwr,'%.2f'),  fx(mC.pwr,'%.2f')};

        cla(axBIn);
        bIn = bar(axBIn,1:N,[abs(mdl.wI) abs(wF) abs(wC)]);
        bIn(1).FaceColor = pal.ink;  bIn(2).FaceColor = pal.fail;  bIn(3).FaceColor = pal.comp;
        legend(axBIn,{'Ideal','Input (faulted)','After algorithm'},'Location','north');
        xlabel(axBIn,'Patch');  ylabel(axBIn,'|w_n|');
        styleAx(axBIn,'Excitation amplitude |w_n|');

        cla(axBGD);
        bG = bar(axBGD,1:3,[[mI.gain mF.gain mC.gain]' [mI.dir mF.dir mC.dir]']);
        bG(1).FaceColor = pal.comp;  bG(2).FaceColor = pal.copper;
        xticks(axBGD,1:3);
        xticklabels(axBGD,{'Ideal','Failed','Compensated'});
        legend(axBGD,{'Gain (dBi)','Directivity (dBi)'},'Location','northeast');
        ylabel(axBGD,'dBi');
        styleAx(axBGD,'Gain and directivity');

        updateInspector(sI.Value,sP.Value);
        updatePolar();
        updateMaps();
        updateSPlots();
    end

    function updateInspector(thv,phv)
        if isempty(cur), return; end
        th = deg2rad(thv);  ph = deg2rad(phv);
        lblIns.Text = sprintf('theta = %.1f deg    phi = %.1f deg',thv,phv);
        ws    = {cur.wI, cur.wF, cur.wC};
        names = {'Ideal','Failed (before)','Compensated (after)'};
        E0 = rowAt(th,ph,mdl.x)*ws{1};
        dB0 = 20*log10(max(abs(E0),1e-12)/mdl.refV);
        rows = cell(3,5);
        for k = 1:3
            E  = rowAt(th,ph,mdl.x)*ws{k};
            Pr = sum(AREA.*abs(mdl.G*ws{k}).^2)*dTh*dPh;
            dB = 20*log10(max(abs(E),1e-12)/mdl.refV);
            gi = 10*log10(4*pi*abs(E)^2/Pr * edEta.Value/100 + eps);
            rows(k,:) = {names{k}, fx(dB,'%.2f'), fx(gi,'%.2f'), ...
                fx(rad2deg(angle(E)),'%.1f'), fx(dB-dB0,'%.2f')};
        end
        tblIns.Data = rows;

        mk = [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)];
        set(SC{1}.q,'XData',mk(1),'YData',mk(2),'ZData',mk(3),'Visible','on');
        set(SC{2}.q,'XData',mk(1),'YData',mk(2),'ZData',mk(3),'Visible','on');

        cla(axInsCut);  hold(axInsCut,'on');
        cols = {pal.ink, pal.fail, pal.comp};
        styl = {'-','--','-'};
        for k = 1:3
            A = abs(rowAt(TCUT(:),ph*ones(numel(TCUT),1),mdl.x)*ws{k})/mdl.refV;
            plot(axInsCut,rad2deg(TCUT(:)),max(20*log10(max(A,1e-6)),DBF), ...
                styl{k},'Color',cols{k},'LineWidth',1.8);
        end
        plot(axInsCut,[thv thv],[DBF 3],':','Color',pal.amber,'LineWidth',2);
        hold(axInsCut,'off');
        xlim(axInsCut,[0 90]);  ylim(axInsCut,[DBF 3]);
        xlabel(axInsCut,'theta (deg)','Color',pal.ink);
        ylabel(axInsCut,'dB re ideal peak','Color',pal.ink);
        legend(axInsCut,{'Ideal','Failed','Compensated','Selected theta'}, ...
            'Location','southwest');
        styleAx(axInsCut,sprintf('E-plane cut at phi = %.1f deg (dB re ideal peak)',phv));
    end

    function updatePolar()
        if isempty(cur), return; end
        phB = edPh0.Value;
        thT = edTh0.Value;
        states = {cur.wI, cur.wF, cur.wC};
        names  = {'Ideal (target)','Failed (before)','Compensated (after)'};
        cols   = {pal.ink, pal.fail, pal.comp};
        styl   = {'-','--','-'};
        mks    = {'p','v','o'};

        % --- beam plane: contains the commanded direction ---
        drawPolarBase(axPA,sprintf('Beam plane, phi = %.0f deg (target inside)',phB));
        hold(axPA,'on');
        hA = gobjects(1,3);
        pk = zeros(3,2);
        for k = 1:3
            [ang,dB] = cutPolar(states{k},phB);
            [x,y] = polarXY(ang,dB);
            hA(k) = plot(axPA,x,y,styl{k},'Color',cols{k},'LineWidth',2);
            [pT,~] = cutPeakTheta(states{k},deg2rad(phB));
            Epk = abs(rowAt(deg2rad(pT),deg2rad(phB),mdl.x)*states{k})/mdl.refV;
            dBpk = 20*log10(max(Epk,1e-6));
            pk(k,:) = [pT dBpk];
            [xp,yp] = polarXY(pT,dBpk);
            plot(axPA,xp,yp,mks{k},'MarkerSize',9,'MarkerFaceColor',cols{k}, ...
                'MarkerEdgeColor',pal.ink);
        end
        th = deg2rad(thT);
        hT = plot(axPA,[0 1.08*sin(th)],[0 1.08*cos(th)],'-.','Color',pal.green,'LineWidth',1.6);
        hold(axPA,'off');
        legend(axPA,[hA hT],[names {'Commanded direction'}],'Location','south','FontSize',8);

        % --- orthogonal plane: beam is outside this plane ---
        drawPolarBase(axPB,sprintf('Orthogonal plane, phi = %.0f deg (beam is outside)',phB+90));
        hold(axPB,'on');
        hB = gobjects(1,3);
        for k = 1:3
            [ang,dB] = cutPolar(states{k},phB+90);
            [x,y] = polarXY(ang,dB);
            hB(k) = plot(axPB,x,y,styl{k},'Color',cols{k},'LineWidth',2);
        end
        hold(axPB,'off');
        legend(axPB,hB,names,'Location','south','FontSize',8);

        lines = {'POLAR COMPARISON (current case)';
            sprintf('Commanded : theta0 = %.1f deg, phi0 = %.1f deg', thT, phB);
            sprintf('Ideal peak   : theta %.2f deg  (%.2f dB)', pk(1,1), pk(1,2));
            sprintf('Failed peak  : theta %.2f deg  (%.2f dB)', pk(2,1), pk(2,2));
            sprintf('Fixed peak   : theta %.2f deg  (%.2f dB)', pk(3,1), pk(3,2));
            sprintf('Shift from target, failed : %+.2f deg', pk(2,1)-pk(1,1));
            sprintf('Shift from target, fixed  : %+.2f deg', pk(3,1)-pk(1,1));
            '';
            'How to read it:';
            ' - Radius = level in dB re ideal peak (outer ring = 0 dB).';
            ' - Black solid = ideal, red dashed = failed, blue = after';
            '   reconfiguration. Markers show each peak.';
            ' - Green dash-dot = commanded direction.';
            ' - Right plot: same state in the perpendicular plane.'};
        txtPo.Value = lines;
    end

    function updateMaps()
        if isempty(cur), return; end
        maps   = {axMI, axMF, axMC};
        ttl    = {'Ideal','Failed (before)','Compensated (after)'};
        states = {cur.wI, cur.wF, cur.wC};
        for k = 1:3
            ax = maps{k};
            F  = mdl.G*states{k};
            A  = abs(reshape(F,size(THg)))/mdl.refV;
            Ec = max(20*log10(max(A,1e-6)),DBF);
            cla(ax);
            imagesc(ax,PHV,THV,Ec);
            set(ax,'YDir','normal');
            colormap(ax,jet(256));
            caxis(ax,[DBF 0]);
            hold(ax,'on');
            plot(ax,mod(edPh0.Value,360),edTh0.Value,'p','MarkerSize',12, ...
                'MarkerFaceColor',pal.amber,'MarkerEdgeColor',pal.ink);
            hold(ax,'off');
            colorbar(ax);
            xlabel(ax,'phi (deg)');
            ylabel(ax,'theta (deg, from vertical)');
            styleAx(ax,ttl{k});
        end
    end

    function updateSPlots()
        if isempty(cur), return; end
        [~,fG] = sBase();
        Sid = sMatrix(false(N,1),true(N,1),0);
        Sfa = sMatrix(cur.fm,~cur.fm,cur.sev);
        Sco = sMatrix(cur.fm,cur.act,cur.sev);
        p = find(strcmp(patchItems(N),ddSP.Value),1);
        q = find(strcmp(patchItems(N),ddSQ.Value),1);
        if isempty(p), p = 1; end
        if isempty(q), q = 2; end

        % 1. reflection of the selected port
        cla(axS1); hold(axS1,'on');
        plot(axS1,fG,dbOf(Sid(:,p,p)),'-','Color',pal.ink,'LineWidth',2);
        plot(axS1,fG,dbOf(Sfa(:,p,p)),'--','Color',pal.fail,'LineWidth',1.8);
        plot(axS1,fG,dbOf(Sco(:,p,p)),'-','Color',pal.comp,'LineWidth',1.8);
        hold(axS1,'off');
        ylim(axS1,[-40 1]);
        xlabel(axS1,'Frequency (GHz)');  ylabel(axS1,'|S| (dB)');
        legend(axS1,{'Ideal','Failed (before)','Compensated (after)'},'Location','south');
        styleAx(axS1,sprintf('Reflection S%d%d of selected port',p,p));

        % 2. average reflected power of the active ports
        cla(axS2); hold(axS2,'on');
        plot(axS2,fG,meanRefl(Sid,true(N,1)),'-','Color',pal.ink,'LineWidth',2);
        plot(axS2,fG,meanRefl(Sfa,~cur.fm),'--','Color',pal.fail,'LineWidth',1.8);
        plot(axS2,fG,meanRefl(Sco,cur.act),'-','Color',pal.comp,'LineWidth',1.8);
        hold(axS2,'off');
        xlabel(axS2,'Frequency (GHz)');  ylabel(axS2,'dB');
        legend(axS2,{'Ideal','Failed (before)','Compensated (after)'},'Location','south');
        styleAx(axS2,'Mean reflected power of active ports (dB)');

        % 3. all ports vs frequency
        switch ddSS.Value
            case 'Failed',  Sh = Sfa;
            case 'Ideal',   Sh = Sid;
            otherwise,      Sh = Sco;
        end
        Mx = zeros(N,numel(fG));
        for k = 1:N
            Mx(k,:) = dbOf(Sh(:,k,k)).';
        end
        cla(axS3);
        imagesc(axS3,fG,1:N,Mx);
        set(axS3,'YDir','normal');
        colormap(axS3,jet(256));
        caxis(axS3,[-40 0]);
        colorbar(axS3);
        xlabel(axS3,'Frequency (GHz)');  ylabel(axS3,'Port');
        styleAx(axS3,sprintf('Reflection |S_kk| (dB), all ports, %s state',ddSS.Value));

        % 4. mutual coupling (Touchstone only)
        cla(axS4);
        if strcmp(srcMode,'file') && p ~= q
            hold(axS4,'on');
            plot(axS4,fG,dbOf(Sid(:,p,q)),'-','Color',pal.ink,'LineWidth',2);
            plot(axS4,fG,dbOf(Sfa(:,p,q)),'--','Color',pal.fail,'LineWidth',1.8);
            plot(axS4,fG,dbOf(Sco(:,p,q)),'-','Color',pal.comp,'LineWidth',1.8);
            hold(axS4,'off');
            ylim(axS4,[-80 0]);
            xlabel(axS4,'Frequency (GHz)');  ylabel(axS4,'|S| (dB)');
            legend(axS4,{'Ideal','Failed','Compensated'},'Location','south');
            styleAx(axS4,sprintf('Mutual coupling |S%d%d| (Touchstone)',p,q));
        else
            xlim(axS4,[0 1]);  ylim(axS4,[0 1]);
            if p == q
                msg = 'Pick two different ports in "Pair with" to see coupling.';
            else
                msg = {'Mutual coupling S_ij needs a Touchstone file.', ...
                       'Load the CST .sNp export (Settings > Load Touchstone).', ...
                       'The parametric model has S11 only.'};
            end
            text(axS4,0.5,0.5,msg,'HorizontalAlignment','center','Color',pal.ink);
            styleAx(axS4,'Mutual coupling S_ij');
        end
    end

    function [S0,fG] = sBase()
        % Base S-matrix before failures: Touchstone data, or a resonator model.
        if strcmp(srcMode,'file') && ~isempty(Sfile)
            S0 = Sfile;
            fG = fFile;
        else
            fG  = linspace(edF.Value - edSpan.Value/2, edF.Value + edSpan.Value/2, NSP).';
            f0  = edF.Value*(1 + edDet.Value/100);
            dlt = (fG - f0)/f0;
            Z   = edRin.Value./(1 + 1j*2*edQp.Value*dlt);       % parallel resonator
            Gm  = (Z - 50)./(Z + 50);                           % reflection vs 50 ohm
            S0  = zeros(numel(fG),N,N);
            for n = 1:N
                S0(:,n,n) = Gm;
            end
        end
    end

    function S = sMatrix(fm,act,sev)
        % Failed port: reflects (1-sev) of its match plus sev of a short.
        % Switched-off healthy port: reflects fully and does not couple.
        [S0,~] = sBase();
        S = S0;
        for n = 1:N
            if fm(n)
                sn = sev;
            elseif ~act(n)
                sn = 1;
            else
                continue
            end
            dg = (1-sn)*S0(:,n,n) - sn;
            S(:,n,:) = (1-sn)*S0(:,n,:);
            S(:,:,n) = (1-sn)*S0(:,:,n);
            S(:,n,n) = dg;
        end
    end

    function onSrcChanged()
        if strcmp(ddSrc.Value,'Touchstone file')
            if isempty(Sfile)
                ddSrc.Value = 'Parametric model';
                srcMode = 'model';
                uialert(fig,'Load a Touchstone file first (button in the S-PARAMETERS section).', ...
                    'No S-parameter file');
            else
                srcMode = 'file';
            end
        else
            srcMode = 'model';
        end
        updateSPlots();
    end

    function loadTouchstone()
        [f,p] = uigetfile({'*.s*p','Touchstone S-parameter (*.sNp)'}, ...
            'Select Touchstone file exported from CST');
        if isequal(f,0), return; end
        try
            [fG,S,Np] = parseTouchstone(fullfile(p,f));
        catch err
            uialert(fig,err.message,'Could not read Touchstone file','Icon','error');
            setStatus('Touchstone import failed. See the message box.');
            return
        end
        Sfile = S;  fFile = fG;  srcMode = 'file';
        ddSrc.Value = 'Touchstone file';
        edN.Value = Np;
        onNChanged();
        tg.SelectedTab = tabS;
        setStatus(sprintf('Loaded %s: %d ports, %d frequency points (%.2f to %.2f GHz). Mutual coupling included.', ...
            f, Np, numel(fG), fG(1), fG(end)));
    end

    function onNChanged()
        N = round(edN.Value);
        edN.Value = N;
        if strcmp(srcMode,'file') && size(Sfile,2) ~= N
            srcMode = 'model';
            ddSrc.Value = 'Parametric model';
        end
        items = patchItems(N);
        lbFail.Items = items;
        lbFail.Value = {};
        ddSP.Items = items;  ddSP.Value = items{1};
        ddSQ.Items = items;  ddSQ.Value = items{min(2,N)};
        caseMask = false(0,N);
        cur = [];
        prepareModel();
        showAll(mdl.wI,mdl.wI,false(N,1),true(N,1),'ideal');
        setStatus(sprintf('Array set to N = %d patches (%s).',N, ...
            ternStr(strcmp(srcMode,'file'),'Touchstone S-matrix','parametric S11 model')));
    end

    function out = fx(v,fmt)
        if isempty(v) || ~isfinite(v)
            out = '--';
        else
            out = sprintf(fmt,v);
        end
    end

    function updateReport()
        if isempty(cur), return; end
        mI = metricsOf(mdl.wI);
        mF = metricsOf(cur.wF);
        mC = metricsOf(cur.wC);
        L = {
            'ARRAY FAILURE & RECONFIGURATION REPORT'
            repmat('=',1,84)
            sprintf('N = %d patches | f = %.2f GHz | eps_r = %.2f | h = %.2f mm | patch L = %.2f mm, W = %.2f mm, eps_eff = %.3f', ...
                N,edF.Value,edEr.Value,edH.Value,mdl.L,mdl.W,mdl.ee)
            sprintf('spacing = %.2f lambda0 (%.2f mm) | steer (theta0, phi0) = (%.1f, %.1f) deg | taper %s', ...
                edD.Value,mdl.dmm,edTh0.Value,edPh0.Value,ddTap.Value)
            sprintf('element model: %s | radiation efficiency %.0f %% | S-source: %s', ...
                elemName(), edEta.Value, ternStr(strcmp(srcMode,'file'),'Touchstone','parametric model'))
            sprintf('failed: %s | severity %.2f | active set (%d patches): %s', ...
                patchList(cur.fm,', ','none'), cur.sev, sum(cur.act), patchList(cur.act,', ','none'))
            sprintf('reconfiguration method: %s', lastMethod)
            repmat('-',1,84)
            sprintf('%-24s %12s %12s %16s','Metric','Ideal','Failed','Compensated')
            sprintf('%-24s %12.2f %12.2f %16.2f','Gain (dBi)',mI.gain,mF.gain,mC.gain)
            sprintf('%-24s %12.2f %12.2f %16.2f','Directivity (dBi)',mI.dir,mF.dir,mC.dir)
            sprintf('%-24s %12.2f %12.2f %16.2f','Peak vs ideal (dB)',mI.gl,mF.gl,mC.gl)
            sprintf('%-24s %12.3f %12.3f %16.3f','Pointing error (deg)',mI.ptErr,mF.ptErr,mC.ptErr)
            sprintf('%-24s %12.1f %12.1f %16.1f','Max sidelobe (dB)',mI.sll,mF.sll,mC.sll)
            sprintf('%-24s %12.2f %12.2f %16.2f','Fidelity to ideal (%%)',mI.corr,mF.corr,mC.corr)
            sprintf('%-24s %12.2f %12.2f %16.2f','Pattern RMS error (%%)',mI.rms,mF.rms,mC.rms)
            repmat('-',1,84)
            'Excitations: |w| / phase (deg)'
            };
        for k = 1:N
            L{end+1,1} = sprintf('P%-3d ideal %.3f/%+7.1f   input %.3f/%+7.1f   after %.3f/%+7.1f', ...
                k,abs(mdl.wI(k)),rad2deg(angle(mdl.wI(k))), ...
                abs(cur.wF(k)),rad2deg(angle(cur.wF(k))), ...
                abs(cur.wC(k)),rad2deg(angle(cur.wC(k)))); %#ok<AGROW>
        end
        txtRep.Value = L(:);
    end

% =========================================================================
%  NESTED: actions
% =========================================================================

    function computeCurrent()
        fm = getFail();
        if all(fm)
            uialert(fig,'Keep at least one healthy patch.','Check failures');
            return
        end
        prepareModel();
        sev = sld.Value;
        wF  = failedWeights(fm,sev);
        R   = reconfigure(fm,sev,false);
        if ~R.ok
            showAll(wF,wF,fm,~fm,'fail');
            setStatus('Not recoverable: fewer than two healthy patches remain.');
            return
        end
        lastMethod = R.method;
        showAll(wF,R.wC,fm,R.act,'comp');
        updateReport();
        setStage('Computed');
        setStatus(sprintf(['Best active set (%d of %d patches, %s): fidelity %.2f %%, ' ...
            'gain %.2f dBi, pointing %.3f deg'],sum(R.act),N,R.method,R.m.corr,R.m.gain,R.m.ptErr));
    end

    function runProcess()
        if busy, return; end
        setBusy(true);
        cu = onCleanup(@()setBusy(false)); %#ok<NASGU>
        try
            fm = getFail();
            if all(fm)
                uialert(fig,'Keep at least one healthy patch.','Check failures');
                return
            end
            prepareModel();
            sev = sld.Value;
            spd = edSpd.Value/1000;
            tg.SelectedTab = tabLive;
            txtLog.Value = {'ARRAY FAILURE & RECONFIGURATION: LIVE PROCESS';''};

            % ---- STEP 1 ----
            setStage('1/6  Ideal array');
            logMsg('STEP 1  Ideal 1 x %d patch array', N);
            logMsg('        f = %.2f GHz, eps_r = %.2f, h = %.2f mm', edF.Value, edEr.Value, edH.Value);
            logMsg('        patch L = %.2f mm, W = %.2f mm, eps_eff = %.3f', mdl.L, mdl.W, mdl.ee);
            logMsg('        spacing = %.2f mm (%.2f lambda0)', mdl.dmm, edD.Value);
            logMsg('        steer theta0 = %.1f deg, phi0 = %.1f deg', edTh0.Value, edPh0.Value);
            if mdl.dmm < mdl.L + 1
                logMsg('        WARNING: spacing close to patch length; coupling not modelled.');
            end
            showAll(mdl.wI,mdl.wI,false(N,1),true(N,1),'ideal');
            animPause(NFR*spd);

            % ---- STEP 2 ----
            setStage('2/6  Failure injection');
            logMsg('STEP 2  Inject failure: %s, severity %.2f', patchList(fm,', ','none'), sev);
            for t = linspace(0,1,NFR)
                wT = failedWeights(fm,sev*t);
                setPanels(SF,wT);
                drawGeom(wT,fm,~fm,'fail');
                animPause(spd);
            end
            wF = failedWeights(fm,sev);
            mU = metricsOf(wF);
            logMsg('        before algorithm: peak %.2f dB, pointing %.3f deg, fidelity %.2f %%', ...
                mU.gl, mU.ptErr, mU.corr);
            showAll(mdl.wI,wF,fm,~fm,'fail');
            animPause(NFR*spd/2);

            % ---- STEP 3 ----
            setStage('3/6  Reconfiguration search');
            R = reconfigure(fm,sev,true);
            if ~R.ok
                setStage('Not recoverable');
                logMsg('        NOT RECOVERABLE: fewer than two healthy patches remain.');
                setStatus('Not recoverable: fewer than two healthy patches remain.');
                return
            end
            lastMethod = R.method;
            logMsg('STEP 3  %s search: %d candidate evaluations', R.method, R.n);
            logMsg('        selected active set (%d patches): fidelity %.2f %%, pointing %.3f deg', ...
                sum(R.act), R.m.corr, R.m.ptErr);
            animPause(NFR*spd/2);

            % ---- STEP 4 ----
            setStage('4/6  Constrained least squares (KKT)');
            logMsg('STEP 4  Solve  [A''A  C''; C  0] [w; lambda] = [A''y; d]');
            logMsg('        active unknowns %d, peak constraint row 1, slope rows 2', R.info.nAct);
            logMsg('        lock weight %g (used %g) in %d pass(es), residual %.3f deg', ...
                edLock.Value, R.info.lock, R.info.pass, R.info.err);
            for t = linspace(0,1,NFR)
                wT = (1-t)*wF + t*R.wC;
                setPanels(SC,wT);
                drawGeom(wT,fm,R.act,'comp');
                animPause(spd);
            end

            % ---- STEP 5 ----
            setStage('5/6  Verification');
            mC = R.m;
            logMsg('STEP 5  Verification against ideal');
            logMsg('        gain     : %.2f dBi  (before algorithm %.2f dBi)', mC.gain, mU.gain);
            logMsg('        peak     : %.2f dB re ideal  (before %.2f dB)', mC.gl, mU.gl);
            logMsg('        pointing : %.3f deg  (before %.3f deg)', mC.ptErr, mU.ptErr);
            logMsg('        fidelity : %.2f %%  (before %.2f %%)', mC.corr, mU.corr);
            logMsg('        sidelobe : %.1f dB  (before %.1f dB)', mC.sll, mU.sll);
            beamProof_compute(mdl.wI, wF, R.wC);
            if mC.ptErr <= TOL
                logMsg('        RESULT: PASS. Beam direction held within %.2f deg.', TOL);
            else
                logMsg('        RESULT: CHECK. Pointing drift above %.2f deg.', TOL);
            end
            animPause(NFR*spd/2);

            % ---- STEP 6 ----
            setStage('6/6  Complete');
            showAll(mdl.wI,wF,fm,R.act,'comp');
            updateReport();
            tg.SelectedTab = tabLive;
            setStatus(sprintf(['Compensated with %d patches (%s): fidelity %.2f %%, gain %.2f dBi, ' ...
                'pointing %.3f deg.'],sum(R.act),R.method,mC.corr,mC.gain,mC.ptErr));
        catch err
            setStatus(['Error: ' err.message]);
            uialert(fig,err.message,'Process error','Icon','error');
            disp(getReport(err,'extended'));
        end
    end

    function runAnalysis()
        if busy, return; end
        setBusy(true);
        cu = onCleanup(@()setBusy(false)); %#ok<NASGU>
        try
            setStage('Failure analysis');
            prepareModel();
            sev = sld.Value;
            [cases,isSample] = makeCases(N,round(edMaxF.Value),round(edK.Value),round(edSeed.Value));
            nC = size(cases,1);
            caseMask = cases;

            rows = cell(nC,12);
            fidB = nan(nC,1);  fidG = nan(nC,1);  fidU = nan(nC,1);
            isExh = false(nC,1);
            lbl = cell(nC,1);
            for c = 1:nC
                fm = cases(c,:).';
                nH = sum(~fm);
                lbl{c} = sprintf('#%d',c);
                wF = failedWeights(fm,sev);
                mU = metricsOf(wF);
                fidU(c) = mU.corr;
                R = reconfigure(fm,sev,false);
                if R.ok
                    if strcmp(R.method,'Greedy')
                        G = R;
                    else
                        G = greedySearch(fm,sev,false);
                    end
                    fidB(c) = R.m.corr;
                    fidG(c) = G.m.corr;
                    isExh(c) = strcmp(R.method,'Exhaustive');
                    v = verdictOf(R.m);
                    rows(c,:) = {c, patchList(fm,', ','none'), nH, R.method, ...
                        sprintf('%d: %s', sum(R.act), patchList(R.act,', ','-')), ...
                        fx(R.m.corr,'%.2f'), fx(R.m.gain,'%.2f'), fx(R.m.ptErr,'%.3f'), ...
                        fx(R.m.sll,'%.1f'), fx(G.m.corr,'%.2f'), v, fx(mU.corr,'%.2f')};
                else
                    rows(c,:) = {c, patchList(fm,', ','none'), nH, 'n/a', 'n/a', '--', '--', '--', ...
                        '--', '--', 'Not recoverable (<2 healthy)', fx(mU.corr,'%.2f')};
                end
                setStatus(sprintf('Analysing case %d of %d ...',c,nC));
            end
            caseTbl = rows;
            tblA.Data = rows;
            tg.SelectedTab = tabAn;

            cla(axBar);
            bb = bar(axBar,1:nC,[fidB fidG fidU]);
            bb(1).FaceColor = pal.comp;
            bb(2).FaceColor = pal.green;
            bb(3).FaceColor = pal.fail;
            xticks(axBar,1:nC);
            if nC <= 40
                xticklabels(axBar,lbl);
            else
                xticklabels(axBar,{});
            end
            ylim(axBar,[0 110]);
            legend(axBar,{'Reconfiguration (best)','Greedy','No fix'},'Location','north');
            ylabel(axBar,'Fidelity (%)');
            xlabel(axBar,'Failure case (see table)');
            styleAx(axBar,'Fidelity per failure case (higher = closer to ideal)');

            ok = ~isnan(fidB);
            L = {'FAILURE & RECONFIGURATION ANALYSIS';
                 sprintf('N = %d | steer (%.0f, %.0f) deg | d = %.2f lambda0 | severity %.2f | eta %.0f %%', ...
                    N,edTh0.Value,edPh0.Value,edD.Value,sev,edEta.Value);
                 sprintf('Cases: %d (%s), up to %d failed per case', nC, ...
                    ternStr(isSample,'random sample','all combinations'), round(edMaxF.Value));
                 repmat('-',1,78)};
            sel1 = find(sum(cases,2) == 1 & ok);
            if ~isempty(sel1)
                [wv,iw] = min(fidB(sel1));
                [bv,ib] = max(fidB(sel1));
                L{end+1,1} = sprintf('Best single failure : %s failed -> fidelity %.2f %%', rows{sel1(ib),2}, bv);
                L{end+1,1} = sprintf('Worst single failure: %s failed -> fidelity %.2f %%', rows{sel1(iw),2}, wv);
            else
                L{end+1,1} = 'No single-failure cases in this set.';
            end
            L{end+1,1} = sprintf('Recoverable cases: %d of %d', sum(ok), nC);
            if any(ok)
                L{end+1,1} = sprintf('Mean fidelity over recoverable cases: %.2f %%', mean(fidB(ok)));
            end
            ex = ok & isExh;
            if any(ex)
                gap = fidB(ex) - fidG(ex);
                L{end+1,1} = sprintf('Exhaustive search used in %d cases; greedy matched it in %d of them', ...
                    sum(ex), sum(abs(gap) < 0.05));
                L{end+1,1} = sprintf('Mean fidelity given up by greedy (where compared): %.2f %%', mean(gap));
            else
                L{end+1,1} = 'Exhaustive search not used here (too many healthy patches). Greedy elimination was used.';
            end
            L{end+1,1} = 'Click any row in the table to load that failure case into the live views.';
            txtBest.Value = L(:);

            setStage('Analysis complete');
            setStatus(sprintf('Analysed %d failure cases. Click a row to load it.',nC));
        catch err
            setStatus(['Error: ' err.message]);
            uialert(fig,err.message,'Analysis error','Icon','error');
            disp(getReport(err,'extended'));
        end
    end

    function caseRowClicked(evt)
        if isempty(evt.Selection) || isempty(caseMask), return; end
        r = evt.Selection(1,1);
        if r < 1 || r > size(caseMask,1), return; end
        fm = caseMask(r,:).';
        setFailSelection(fm);
        computeCurrent();
        tg.SelectedTab = tabLive;
    end

    function compareGreedy()
        fm = getFail();
        if all(fm)
            uialert(fig,'Keep at least one healthy patch.','Check failures');
            return
        end
        prepareModel();
        sev = sld.Value;
        nH = sum(~fm);
        if nH > EXH_HARD_H
            txtBest.Value = {'GREEDY vs EXHAUSTIVE', repmat('-',1,78), ...
                sprintf('Exhaustive search over %d healthy patches is not practical.', nH), ...
                'Use more failures, or reduce N, to compare the two methods.'};
            tg.SelectedTab = tabAn;
            return
        end
        tA = tic;
        Rx = exhaustiveSearch(fm,sev,false);
        tEx = toc(tA);
        tB = tic;
        Rg = greedySearch(fm,sev,false);
        tGr = toc(tB);
        if ~Rx.ok
            setStatus('Not recoverable: fewer than two healthy patches remain.');
            return
        end
        if abs(Rx.m.corr - Rg.m.corr) < 0.05
            verdict = 'matched the exhaustive optimum';
        else
            verdict = sprintf('fell short of the exhaustive optimum by %.2f %%', Rx.m.corr - Rg.m.corr);
        end
        L = {'GREEDY vs EXHAUSTIVE (current failure case)'; repmat('-',1,78);
            sprintf('Exhaustive: %2d patches {%s} fidelity %6.2f %%  pointing %6.3f deg  evals %4d  time %.3f s', ...
                sum(Rx.act),Rx.label,Rx.m.corr,Rx.m.ptErr,Rx.n,tEx);
            sprintf('Greedy    : %2d patches {%s} fidelity %6.2f %%  pointing %6.3f deg  evals %4d  time %.3f s', ...
                sum(Rg.act),Rg.label,Rg.m.corr,Rg.m.ptErr,Rg.n,tGr);
            sprintf('Greedy %s.',verdict);
            'Exhaustive is the optimum over every healthy subset.'};
        txtBest.Value = L(:);
        tg.SelectedTab = tabAn;
        setStatus('Greedy vs exhaustive comparison written to the Failure Analysis tab.');
    end

    function beamProofButton()
        if isempty(cur)
            setStatus('Run a case first (Compute current case).');
            return
        end
        beamProof_compute(cur.wI,cur.wF,cur.wC);
        tg.SelectedTab = tabProof;
        setStatus(sprintf('Beam proof: max compensated deviation %.2f deg over %d cuts.', ...
            max(proof.devC),numel(PCUT)));
    end

    function beamProof_compute(wI_,wF_,wC_)
        nC = numel(PCUT);
        tI = zeros(nC,1);  tF = tI;  tC = tI;
        for k = 1:nC
            [tI(k),~] = cutPeakTheta(wI_,PCUT(k));
            [tF(k),~] = cutPeakTheta(wF_,PCUT(k));
            [tC(k),~] = cutPeakTheta(wC_,PCUT(k));
        end
        proof = struct('phi',rad2deg(PCUT(:)),'tI',tI,'tF',tF,'tC',tC, ...
            'devF',abs(tF-tI),'devC',abs(tC-tI));
        cla(axPk);  hold(axPk,'on');
        plot(axPk,proof.phi,proof.tI,'-','Color',pal.ink,'LineWidth',2);
        plot(axPk,proof.phi,proof.tF,'--','Color',pal.fail,'LineWidth',1.8);
        plot(axPk,proof.phi,proof.tC,'-o','Color',pal.comp,'LineWidth',1.8,'MarkerSize',4);
        hold(axPk,'off');
        xlim(axPk,[0 345]);  ylim(axPk,[0 90]);
        xticks(axPk,0:45:345);
        xlabel(axPk,'phi cut (deg)');  ylabel(axPk,'peak theta (deg)');
        legend(axPk,{'Ideal','Failed','Compensated'},'Location','best');
        styleAx(axPk,'Beam peak theta vs azimuth cut');

        cla(axDv);  hold(axDv,'on');
        bar(axDv,proof.phi,[proof.devF proof.devC],'grouped');
        yline(axDv,TOL,'--','Color',pal.green,'LineWidth',1.5);
        hold(axDv,'off');
        xlim(axDv,[-10 355]);
        xticks(axDv,0:45:345);
        xlabel(axDv,'phi cut (deg)');  ylabel(axDv,'|delta theta| vs ideal (deg)');
        legend(axDv,{'Failed','Compensated','Tolerance'},'Location','northeast');
        styleAx(axDv,'Beam-angle deviation vs ideal (every cut)');

        rows = cell(nC,7);
        for k = 1:nC
            stat = 'OK';
            if proof.devC(k) > TOL, stat = 'CHECK'; end
            rows(k,:) = {sprintf('%.0f',proof.phi(k)), fx(tI(k),'%.2f'), fx(tF(k),'%.2f'), ...
                fx(tC(k),'%.2f'), fx(proof.devF(k),'%.2f'), fx(proof.devC(k),'%.2f'), stat};
        end
        tblProof.Data = rows;
    end

    function exportCB()
        if isempty(cur)
            setStatus('Run a case first (Compute current case or Live process).');
            return
        end
        [f,p] = uiputfile({'*.mat','MATLAB data (*.mat)'},'Save results', ...
            'array_failure_studio_N.mat');
        if isequal(f,0), return; end
        base = fullfile(p,regexprep(f,'\.mat$',''));

        dBI = 20*log10(max(abs(mdl.G*cur.wI),1e-12)/mdl.refV);
        dBF = 20*log10(max(abs(mdl.G*cur.wF),1e-12)/mdl.refV);
        dBC = 20*log10(max(abs(mdl.G*cur.wC),1e-12)/mdl.refV);

        results = struct('cur',cur,'model',mdl,'table',{caseTbl}, ...
            'proof',proof,'settings',struct('N',N,'f_GHz',edF.Value,'eps_r',edEr.Value, ...
            'h_mm',edH.Value,'d_lambda',edD.Value,'theta0',edTh0.Value,'phi0',edPh0.Value, ...
            'efficiency',edEta.Value,'element',elemName(),'S_source',srcMode, ...
            'search',ddSearch.Value,'method',lastMethod));
        save([base '.mat'],'results');

        C = [{'theta_deg','phi_deg','ideal_dB','failed_dB','compensated_dB'}; ...
             num2cell([rad2deg(TH) rad2deg(PH) dBI dBF dBC])];
        writecell(C,[base '_pattern_3D.csv']);

        if ~isempty(caseTbl)
            hdr = tblA.ColumnName;
            writecell([hdr(:).'; caseTbl],[base '_failure_analysis.csv']);
        end
        if ~isempty(proof)
            P = [{'phi_deg','theta_ideal','theta_failed','theta_comp','dev_failed','dev_comp'}; ...
                 num2cell([proof.phi proof.tI proof.tF proof.tC proof.devF proof.devC])];
            writecell(P,[base '_beam_proof.csv']);
        end

        [~,fG] = sBase();
        Sst = sMatrix(cur.fm,cur.act,cur.sev);
        Sdb = zeros(numel(fG),N);
        for k = 1:N
            Sdb(:,k) = dbOf(Sst(:,k,k));
        end
        hdrS = [{'freq_GHz'}, arrayfun(@(k)sprintf('S%d%d_dB',k,k),1:N,'UniformOutput',false)];
        writecell([hdrS; num2cell([fG Sdb])],[base '_S11_current.csv']);

        exportgraphics(fig,[base '_dashboard.png'],'Resolution',200);

        fid = fopen([base '_report.txt'],'w');
        fprintf(fid,'%s\n',txtRep.Value{:});
        fclose(fid);
        setStatus(['Exported ' base ' (.mat, pattern CSV, analysis CSV, beam proof CSV, S11 CSV, PNG, report).']);
    end

end  % ====================== end of main function ======================

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function items = patchItems(N)
items = arrayfun(@(k)sprintf('P%d',k),1:N,'UniformOutput',false);
end

function [cases,isSample] = makeCases(N,maxF,K,seed)
% Every failure combination (1..maxF failed) if there are at most K of them,
% otherwise a reproducible random sample of K distinct combinations.
maxF = min(maxF,N-1);
total = 0;
for k = 1:maxF
    total = total + nchoosek(N,k);
end
cases = false(0,N);
isSample = total > K;
if ~isSample
    for k = 1:maxF
        C = nchoosek(1:N,k);
        for r = 1:size(C,1)
            m = false(1,N);
            m(C(r,:)) = true;
            cases(end+1,:) = m; %#ok<AGROW>
        end
    end
else
    st = rng;
    rng(seed);
    tries = 0;
    while size(cases,1) < K && tries < 100*K
        tries = tries + 1;
        k = randi(maxF);
        idx = randperm(N,k);
        m = false(1,N);
        m(idx) = true;
        if ~any(all(cases == m,2))
            cases(end+1,:) = m; %#ok<AGROW>
        end
    end
    rng(st);
end
[~,ord] = sort(sum(cases,2),'ascend');
cases = cases(ord,:);
end

function subs = subsetsOf(H)
% all subsets of H with at least two elements
subs = {};
nh = numel(H);
for mask = 1:(2^nh-1)
    bits = logical(bitget(mask,1:nh));
    if sum(bits) >= 2
        subs{end+1} = H(bits); %#ok<AGROW>
    end
end
end

function s = patchList(mask,sep,emptyTxt)
idx = find(mask).';
if isempty(idx)
    s = emptyTxt;
else
    parts = arrayfun(@(k)sprintf('P%d',k),idx,'UniformOutput',false);
    s = strjoin(parts,sep);
end
end

function v = verdictOf(m)
if m.ptErr > 0.5
    v = 'Pointing drift';
elseif m.corr >= 95
    v = 'Excellent';
elseif m.corr >= 85
    v = 'Good';
else
    v = 'Limited';
end
end

function s = ternStr(c,a,b)
if c
    s = a;
else
    s = b;
end
end

function fix3d(ax)
try
    ax.Interactions = [rotateInteraction zoomInteraction];
catch
end
end

function boxp(ax,x0,x1,y0,y1,z0,z1,c,al)
V = [x0 y0 z0; x1 y0 z0; x1 y1 z0; x0 y1 z0; x0 y0 z1; x1 y0 z1; x1 y1 z1; x0 y1 z1];
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch(ax,'Faces',F,'Vertices',V,'FaceColor',c,'FaceAlpha',al,'EdgeColor',[0.35 0.35 0.35]);
end

function a = angDist(t1,p1,t2,p2)
c = sin(t1).*sin(t2).*cos(p1-p2) + cos(t1).*cos(t2);
a = rad2deg(acos(min(max(c,-1),1)));
end

function d = dbOf(z)
d = 20*log10(max(abs(z),1e-6));
end

function pm = meanRefl(S,act)
% mean reflected power of the active ports, in dB, vs frequency
idx = find(act);
if isempty(idx)
    pm = nan(size(S,1),1);
    return
end
p = zeros(size(S,1),1);
for k = idx(:).'
    p = p + abs(S(:,k,k)).^2;
end
pm = 10*log10(p/numel(idx) + eps);
end

function [x,y] = polarXY(angDeg,dB)
% dB -> radius (outer ring = 0 dB, centre = -40 dB); angle from vertical
r = min(max((dB+40)/40,0),1);
a = deg2rad(angDeg);
x = r.*sin(a);
y = r.*cos(a);
end

function drawPolarBase(ax,ttl)
% Polar-style template drawn with plain axes so it works in a uifigure.
cla(ax);  hold(ax,'on');
tt = linspace(-pi/2,pi/2,181);
for L = [0 -10 -20 -30]
    r = (L+40)/40;
    plot(ax,r*sin(tt),r*cos(tt),'-','Color',[0.84 0.84 0.87],'LineWidth',0.8);
    text(ax,0.02,r+0.015,sprintf('%d dB',L),'FontSize',8,'Color',[0.45 0.45 0.45]);
end
for a = -90:30:90
    plot(ax,[0 sin(deg2rad(a))],[0 cos(deg2rad(a))],'-','Color',[0.86 0.86 0.88],'LineWidth',0.8);
    text(ax,1.13*sin(deg2rad(a)),1.13*cos(deg2rad(a)),sprintf('%d%s',a,char(176)), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.3 0.3 0.3]);
end
hold(ax,'off');
xlim(ax,[-1.3 1.3]);  ylim(ax,[-0.05 1.25]);
axis(ax,'equal');
axis(ax,'off');
title(ax,ttl,'Color',[0.08 0.08 0.10],'FontWeight','bold');
end

function C = parseCSTCut(fname)
% Reads a CST 1D farfield ASCII export: skips the text header, finds the
% numeric block, and extracts theta, magnitude and optional phase.
txt   = fileread(fname);
lines = regexp(txt,'\r\n|\n|\r','split');
hdr = '';
blk = {};
started = false;
for k = 1:numel(lines)
    s = strtrim(lines{k});
    if isempty(s)
        if started, break; end
        continue
    end
    tok = regexp(s,'[\s,;\t]+','split');
    tok = tok(~cellfun(@isempty,tok));
    v   = str2double(tok);
    if numel(v) >= 2 && all(isfinite(v))
        blk{end+1} = v(:).'; %#ok<AGROW>
        started = true;
    else
        if started, break; end
        hdr = [hdr ' ' s]; %#ok<AGROW>
    end
end
if isempty(blk)
    error('parseCSTCut:noData', ...
        'No numeric data found in %s. Export a CST ASCII 1D cut (theta, magnitude, phase).', fname);
end
nc  = mode(cellfun(@numel,blk));
blk = blk(cellfun(@numel,blk) == nc);
D   = vertcat(blk{:});
rawTh = D(:,1);
mag   = D(:,2);
if nc >= 3
    ph = D(:,3);
    hasPhase = true;
else
    ph = zeros(size(mag));
    hasPhase = false;
end
if sum(rawTh >= 0) >= 3
    sel = rawTh >= 0 & rawTh <= 90;
else
    sel = abs(rawTh) <= 90;
end
theta = abs(rawTh(sel));
mag   = mag(sel);
ph    = ph(sel);
h = lower(hdr);
isdB = contains(h,'db') || any(mag < 0);
if contains(h,'rad') && ~contains(h,'deg')
    ph = rad2deg(ph);
end
if isdB
    amp = 10.^(mag/20);
else
    amp = mag;
end
Fc = amp.*exp(1j*deg2rad(ph));
[tu,~,ic] = unique(round(theta,4));
cnt = accumarray(ic,1);
Fu  = (accumarray(ic,real(Fc))./cnt) + 1j*(accumarray(ic,imag(Fc))./cnt);
if numel(tu) < 3
    error('parseCSTCut:few','Only %d usable angle points found.',numel(tu));
end
Fu = Fu/max(abs(Fu));
C.theta      = tu(:);
C.fr         = griddedInterpolant(tu(:),real(Fu(:)),'linear','nearest');
C.fi         = griddedInterpolant(tu(:),imag(Fu(:)),'linear','nearest');
C.isdB       = isdB;
C.hasPhase   = hasPhase;
C.npts       = numel(theta);
C.thetaRange = [tu(1) tu(end)];
end

function [fG,S,Np] = parseTouchstone(fname)
% Reads a Touchstone 1.x file (.sNp) from CST. Frequencies are returned in GHz.
tok = regexp(fname,'\.s(\d+)p$','tokens','ignorecase');
if isempty(tok)
    error('parseTouchstone:ext','File must have a .sNp extension (for example .s20p).');
end
Np  = str2double(tok{1}{1});
txt = fileread(fname);
lines = regexp(txt,'\r\n|\n|\r','split');
unitS = 1;          % GHz
fmt   = 'MA';
vals  = [];
for k = 1:numel(lines)
    s = strtrim(regexprep(lines{k},'!.*$',''));
    if isempty(s), continue; end
    if s(1) == '#'
        up = upper(strsplit(strtrim(s(2:end))));
        for j = 1:numel(up)
            switch up{j}
                case 'HZ',  unitS = 1e-9;
                case 'KHZ', unitS = 1e-6;
                case 'MHZ', unitS = 1e-3;
                case 'GHZ', unitS = 1;
                case {'MA','DB','RI'}, fmt = up{j};
            end
        end
        continue
    end
    t = regexp(s,'[\s,]+','split');
    v = str2double(t);
    if any(isnan(v))
        error('parseTouchstone:token','Non-numeric entry in data line: %s', s);
    end
    vals = [vals v(:).']; %#ok<AGROW>
end
per = 1 + 2*Np*Np;
if isempty(vals) || mod(numel(vals),per) ~= 0
    error('parseTouchstone:format', ...
        'Found %d numbers, which does not fit %d ports (%d per frequency).', numel(vals), Np, per);
end
D  = reshape(vals,per,[]).';
fG = D(:,1)*unitS;
A  = D(:,2:end);
a  = A(:,1:2:end);
b  = A(:,2:2:end);
nF = size(D,1);
S  = zeros(nF,Np,Np);
for k = 1:nF
    switch fmt
        case 'RI'
            Z = a(k,:) + 1j*b(k,:);
        case 'DB'
            Z = 10.^(a(k,:)/20).*exp(1j*deg2rad(b(k,:)));
        otherwise
            Z = a(k,:).*exp(1j*deg2rad(b(k,:)));
    end
    if Np == 2
        M = reshape(Z,2,2);          % 2-port order: S11 S21 S12 S22 (column-major)
    else
        M = reshape(Z,Np,Np).';     % N-port order: row by row
    end
    S(k,:,:) = M;
end
end