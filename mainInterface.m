function varargout = mainInterface(varargin)
% MAININTERFACE MATLAB code for mainInterface.fig
%      MAININTERFACE, by itself, creates a new MAININTERFACE or raises the existing
%      singleton*.
%
%      H = MAININTERFACE returns the handle to a new MAININTERFACE or the handle to
%      the existing singleton*.
%
%      MAININTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAININTERFACE.M with the given input arguments.
%
%      MAININTERFACE('Property','Value',...) creates a new MAININTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mainInterface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mainInterface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mainInterface

% Last Modified by GUIDE v2.5 29-Mar-2016 10:16:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mainInterface_OpeningFcn, ...
    'gui_OutputFcn',  @mainInterface_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mainInterface is made visible.
function mainInterface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mainInterface (see VARARGIN)

% Choose default command line output for mainInterface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mainInterface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mainInterface_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function varedit_Callback(hObject, ~, handles)
str=get(hObject,'string');
num=str2double(str);
vartag=regexp(get(hObject,'tag'),'[^_]+$','match');
if isnan(num)
    set(hObject,'string',[]);
    setappdata(handles.figure_main,char(vartag),[]);
else
    setappdata(handles.figure_main,char(vartag),num);
end

function varcheckbox_Callback(hObject,~,handles)
vartag=regexp(get(hObject,'tag'),'[^_]+$','match');
hedit=handles.(['edit_',vartag{1}]);
if get(hObject,'value')==0
    set(hedit,'enable','off','string',[]);
else
    set(hedit,'enable','on');
end
varedit_Callback(hedit, [], handles);

function checkbox_hW_Callback(hObject,~,handles)
if get(hObject,'value')==0
    set(handles.edit_hW,'enable','off','string',[]);
    set(handles.edit_hL,'enable','on');
else
    set(handles.edit_hW,'enable','on');
    set(handles.edit_hL,'enable','off','string',[]);
end
varedit_Callback(handles.edit_hW,[],handles);
varedit_Callback(handles.edit_hL,[],handles);

function checkbox_lW_Callback(hObject,~,handles)
if get(hObject,'value')==0
    set(handles.edit_lW,'enable','off','string',[]);
    set(handles.edit_lWD,'enable','on');
else
    set(handles.edit_lW,'enable','on');
    set(handles.edit_lWD,'enable','off','string',[]);
end
varedit_Callback(handles.edit_lW,[],handles);
varedit_Callback(handles.edit_lWD,[],handles);

function checkbox_YanXing_Callback(hObject,~,handles)
if get(hObject,'value')==0
    set(handles.popupmenu_YanXing,'enable','off','value',1);
    popupmenu_YanXing_Callback(handles.popupmenu_YanXing,[],handles);
else
    set(handles.popupmenu_YanXing,'enable','on');
end

function popupmenu_YanXing_Callback(hObject,~,handles)
if get(hObject,'value')==1
    set(handles.text_hn,'enable','off');
    set(handles.edit_hn,'enable','off');
else
    set(handles.text_hn,'enable','on');
    set(handles.edit_hn,'enable','on');
end

function checkbox_ho_hwso_Callback(hObject,~,handles)
tg=regexp(get(hObject,'tag'),'[^_]+$','match');
tg=tg{1};
tgc={'ho','hwso'};
tgc=tgc{~strcmp(tg,tgc)};
if get(hObject,'value')==0
    set(handles.(['edit_',tg]),'enable','off','string',[]);
else
    set(handles.(['edit_',tg]),'enable','on');
    set(handles.(['checkbox_',tgc]),'value',0);
    set(handles.(['edit_',tgc]),'enable','off','string',[]);
end
setappdata(handles.figure_main,tg,[]);
setappdata(handles.figure_main,tgc,[]);

function popupmenu_BanXing_Callback(hObject,~,handles)
if get(hObject,'value')==2
    set(handles.popupmenu_FuFaJiHe,'value',1,'enable','off');
    set(handles.text_delta,'enable','on');
    set(handles.edit_delta,'enable','on');
    set(handles.edit_Fo,'enable','off','string',[]);
    set(handles.text_Fo,'enable','off');
    set(handles.checkbox_fkt,'enable','inactive','value',1);
    set(handles.edit_fkt,'enable','on');
    % 流体力学部分
    set(handles.checkbox_fhFo,'enable','on','value',0);
    set(handles.edit_fhFo,'enable','off','string',[]);
    set(handles.checkbox_co,'enable','on');
    set(handles.checkbox_fhfloodL,'enable','off','value',0);
    set(handles.edit_fhfloodL,'enable','off','string',[]);
    set(handles.text_hyK,'enable','off');
    set(handles.edit_hyK,'enable','off','string',[]);
    set(handles.checkbox_CF,'enable','off','value',0);
    set(handles.edit_CF,'enable','off','string',[]);
    set(handles.checkbox_flh,'enable','off','value',0);
    set(handles.edit_flh,'enable','off','string',[]);
    set(handles.checkbox_fheV,'enable','on');
    set(handles.checkbox_hfrL,'enable','on');
else
    set(handles.popupmenu_FuFaJiHe,'enable','on');
    set(handles.text_delta,'enable','off');
    set(handles.edit_delta,'enable','off','string',[]);
    set(handles.edit_Fo,'enable','on');
    set(handles.text_Fo,'enable','on');
    set(handles.checkbox_fkt,'enable','on');
    % 流体力学部分
    set(handles.checkbox_fhFo,'enable','inactive','value',1);
    set(handles.edit_fhFo,'enable','on');
    set(handles.checkbox_co,'enable','off','value',0);
    set(handles.edit_co,'enable','off','string',[]);
    set(handles.checkbox_fhfloodL,'enable','on');
    set(handles.text_hyK,'enable','on');
    set(handles.edit_hyK,'enable','on');
    set(handles.checkbox_CF,'enable','on');
    set(handles.edit_CF,'enable','on');
    set(handles.checkbox_flh,'enable','on');
    set(handles.checkbox_fheV,'enable','off');
    set(handles.edit_fheV,'enable','off','string',[]);
    set(handles.checkbox_hfrL,'enable','off');
    set(handles.edit_hfrL,'enable','off','string',[]);
end

function dlg_errvar(emptyvar)
he=errordlg({'以下参数未设置:';emptyvar},'参数错误!','modal');
hm=findobj(he,'tag','MessageBox');
set(hm,'fontsize',14);

function dlg_warnvar(varargin)
if nargin==1
    he=warndlg(varargin{1},'','modal');
else
    he=warndlg(varargin{1},varargin{2},'modal');
end
hm=findobj(he,'tag','MessageBox');
set(hm,'fontsize',14);

function [varargout]=geteditVar(cstr,handles)
% 检查参数函数
varargout=cell(1,length(cstr));
for ii=1:length(cstr)
    tmv=getappdata(handles.figure_main,cstr{ii});
    if isempty(tmv)
        dlg_errvar(cstr{ii});
        return
    else
        varargout{ii}=tmv;
    end
end

function [C20H,scLV] = SmithCorrelation( HT,hL,L,V,rhoL,rhoV )
%史密斯关联图其回归式
% HT 塔板间距 m
% hL 板上液层高度 m
% L 塔内液体流量 {L/V}
% V 塔内气体流量 {V/V}
% rhoL 液相密度 {rhoL/rhoV}
% rhoV 气相密度 {rhoV/rhoV}
%%
scH=HT-hL;
scLV=L./V.*sqrt(rhoL./rhoV);
C20H=exp((-4.531+1.6562.*scH+5.5496.*scH.^2-6.4695.*scH.^3)* ...
    ones(size(scLV)) + ...
    (-0.474675+0.079.*scH-1.39.*scH.^2+1.3212.*scH.^3)*log(scLV) + ...
    (-0.07291+0.088307.*scH-0.49123.*scH.^2+0.43196.*scH.^3)* ...
    (log(scLV)).^2);
if (nargout==2)
    C20H=scH;
end

function [ tesVs ] = Sawtoothweir( tesLs,lW,hn,yfVsfx)
%齿形堰分段计算函数
%  Ls   数组，计算点液体流量 m3/s
% lW 堰长 m
% hn 齿深 m
% yfVsfx 函数句柄，参数形式 @(cxhOWf)@(Ls)
%         其中 参数 cxhOWf 为 @(Ls) 函数句柄
%%
tesVs=tesLs;
crisLs=zeros(1,2);
crisLs(1)=lW*hn^(3/2)/0.0442^(5/2)/3600;
crisLs(2)=2646*lW*hn^(3/2)/3600;
%%
cxhOW2=@(hn,lW)@(Ls) ...
    0.0442*(Ls*3600*hn/lW)^(2/5);
cxhOWf2=cxhOW2(hn,lW); % 齿形堰 未过顶
% cxhOWf2=@(Ls) 0;
cxhOW25=@(hn)@(Ls) hn+0;
cxhOWf25=cxhOW25(hn); % 齿形堰 齐顶
cxyLs=@(lW,hn)@(Ls)@(hOW) 2464*lW/hn* ...
    (hOW.^(5/2)-(hOW-hn).^(5/2)) -3600*Ls;
cxyLsp=@(lW,hn)@(hOW) 2464*lW/hn*(5/2)*(hOW.^(3/2)-(hOW-hn).^(3/2)) ;
cxyLsL=cxyLs(lW,hn);
cxyLspf=cxyLsp(lW,hn);
cxhOW3=@(cxyLsL,cxyLspf,ihOW)@(Ls) ...
    NewtonIteration(ihOW,1e-8,cxyLsL(Ls),cxyLspf);
cxhOWf3=cxhOW3(cxyLsL,cxyLspf,2*hn);% 齿形堰 过顶
% cxhOWf3=@(Ls)0;
%%
roco=size(tesLs);
for ii=1:(roco(1)*roco(2))
    if tesLs(ii) < crisLs(1)
        yfVsf=yfVsfx(cxhOWf2);
    elseif tesLs(ii) > crisLs(2)
        yfVsf=yfVsfx(cxhOWf3);
    else
        yfVsf=yfVsfx(cxhOWf25);
    end
    tesVs(ii)=yfVsf(tesLs(ii));
end

function tespzyE=weircorrectionfactor(tesLs,D,lW)
tespzyE=tesLs;
bet=(D/lW)^2;
xLsfx=@(bet,lW)@(Ls)@(x) x.*(sqrt(bet-x.^3)-sqrt(bet-1))- ...
    nthroot((3600/2147*Ls).^2/lW^5,3);
xLsf=xLsfx(bet,lW);
xEpfx=@(bet)@(x) sqrt(bet-x^3)-sqrt(bet-1)-3/2*x^3/sqrt(bet-x^3);
xEpf=xEpfx(bet);
roco=size(tesLs);
for ii=1:(roco(1)*roco(2))
    xEf=xLsf(tesLs(ii));
    tespzyE(ii)=1/NewtonIteration(1,1e-4,xEf,xEpf);
end

function CAF0 = vapourcapacityfactor(HT,rhoV)
AD0=[-1.4982473955,-2.4192704887e-1,-1.0356208239e-1,2.5273175766e-2;
    3.5944717249e-2,-1.1336244695e-1,5.1983585356e-2,-7.5266626935e-3;
    -1.2669642026e-3,3.2363453585e-3,-1.4380132494e-3,2.0385337432e-4;
    1.0833195587e-5,-2.6578038392e-5,1.1448458759e-5,-1.5881537563e-6];
AD1=[-2.0061624669,-4.4536753835e-2,-1.5939275543e-2,-4.1818163355e-3;
    3.6938451245,-2.7915912829,1.0875796258,-1.392892479e-1;
    -9.830397717,9.3192984056,-4.1749660055,5.8801828854e-1;
    6.2752931905,-6.31166623,2.8579088752,-4.0378014777e-1];
tsfx=@(AD)@(Ts) AD(1)+AD(2)/Ts+AD(3)./Ts.^2+AD(4)/Ts.^3;
if rhoV < 1
    AD=AD0;
else
    AD=AD1;
end
af=tsfx(AD(1,:));
bf=tsfx(AD(2,:));
cf=tsfx(AD(3,:));
df=tsfx(AD(4,:));
CAF0=exp(af(HT) + bf(HT)/rhoV + cf(HT)/rhoV.^2 + ...
    df(HT)/rhoV.^3);

function pushbutton_SizeOfPlate(~, ~, handles)
% %工艺尺寸计算按钮回调函数

% 固定与选定参数检查
[Vs,Ls,rhoV,rhoL,sigma, ...
    saf,HT,uop,do,Wc,Ws]=geteditVar(...
    {'Vs','Ls','rhoV','rhoL','sigma', ...
    'saf','HT','uop','do','Wc','Ws'}, ...
    handles);
if isempty(Ws)
    return % 参数缺失终止回调
end
svvarx=@(hh)@(nvar,val)setappdata(hh,nvar,val);
svvar=svvarx(handles.figure_main); % 保存变量函数
%% 主要工艺尺寸
hL=getappdata(handles.figure_main,'hL');
if get(handles.checkbox_hW,'value')==1 && isempty(hL)
    hL=0.05; % 取一初始值用于递归初值
elseif isempty(hL)
    geteditVar({'hL'},handles); % 显示参数缺失对话框
    return
end
% 塔径
[scH,scLV]=SmithCorrelation(HT,hL,Ls,Vs,rhoL,rhoV );
if get(handles.checkbox_C20,'value')==0
    C20=SmithCorrelation(HT,hL,Ls,Vs,rhoL,rhoV ); %负荷参数
else
    C20=geteditVar({'C20'},handles);
    if isempty(C20)
        return
    end
end
C=C20.*(sigma/20)^0.2;
umax=C.*sqrt((rhoL-rhoV)./rhoV);
u=saf.*umax; % 空塔气速
if get(handles.checkbox_D,'value')==0
    D=sqrt(4*Vs/pi./u); %塔径
    set(handles.edit_D,'string',sprintf('%1.3f',D));
else
    D=geteditVar({'D'},handles);
    if isempty(D)
        return
    end
end
AT=pi/4*D.^2; %塔截面积
u=Vs./AT; % 实际空塔气速
% 溢流装置
if get(handles.checkbox_lW,'value')==0
    lWD=geteditVar({'lWD'},handles);
    if isempty(lWD)
        return
    end
    lW=lWD.*D;
else
    lW=geteditVar({'lW'},handles);
    if isempty(lW)
        return
    end
    lWD=lW*D;
end
tf.YanXing=get(handles.popupmenu_YanXing,'value');
% 单溢流平直堰
bet=(D/lW)^2;
kap=nthroot(Ls^2/lW^5*(3600/2147)^2,3);
zetfx=@(bet)@(zet) 2/3*(bet-zet-sqrt(bet-zet)*sqrt(bet-1));
zetf=zetfx(bet);
zet=iterationSim(zetf,1,1e-6);
kapmax=nthroot(zet,3)*(sqrt(bet-zet)-sqrt(bet-1));
ELsmax=sqrt(kap^3*lW^5/(3600/2147)^2); % 液流收缩系数推荐线最大流量
if kap > kapmax
    warnLs=sprintf('Ls = %.3e 过大',Ls);
    dlg_warnvar(warnLs);
    return
end
if get(handles.checkbox_pzyE,'value')==0
    pzyE=weircorrectionfactor(Ls,D,lW);
else
    pzyE=geteditVar({'pzyE'},handles);
    if isempty(pzyE)
        return
    end
end
hOW=2.84/1000*pzyE*(Ls*3600/lW)^(2/3);
if tf.YanXing==2
    hn=geteditVar({'hn'},handles);
    if isempty(hn)
        return
    end
    cxyhOWcal=@(cxhOWf)@(Ls) cxhOWf(Ls);
    hOW=Sawtoothweir(Ls,lW,hn,cxyhOWcal);
elseif hOW < 0.006
    warnhOW=sprintf(['hOW=%0.4g < 0.006 (m)\n', ...
        '不适用于平直堰，宜改用齿形堰'],hOW);
end
if get(handles.checkbox_hW,'value')==0
    hW=hL-hOW; % 堰高
else
    hW=geteditVar({'hW'},handles);
    if isempty(hW)
        return
    end
    if abs(hL-hOW-hW)  >1e-5
        hL=hW+hOW;
        setappdata(handles.figure_main,'hL',hL);
        pushbutton_SizeOfPlate([], [], handles);
        return
    else
        set(handles.edit_hL,'string',sprintf('%0.2e',hL));
    end
end
if exist('warnhOW','var')
    dlg_warnvar(warnhOW,'Warning')
end
if get(handles.checkbox_Wd,'value')==0
    Wd=D/2-sqrt((D/2).^2-(lW/2).^2); % 弓形降液管宽度
else
    Wd=geteditVar({'Wd'},handles);
    if isempty(Wd)
        return
    end
end
if get(handles.checkbox_Af,'value')==0
    % 降液管截面积
    Af=asin(lW/D)*(D/2).^2- sqrt((D/2).^2-(lW/2).^2).*lW/2;
else
    Af=geteditVar({'Af'},handles);
    if isempty(Af)
        return
    end
end
theta=Af*HT/Ls; % 降液管内液体停留时间 s
if theta<3
    dlg_warnvar(sprintf('theta=%.2f s <3s ,\n 应调整堰长lW',theta), ...
        'Waring');
end
if get(handles.checkbox_hwso,'value')==1
    hwso=geteditVar({'hwso'},handles);
    if isempty(hwso)
        return
    end
    ho=hW-hwso; % 液封算法
elseif get(handles.checkbox_ho,'value')==0
    ho=Ls/lW/uop; % 降液管底隙高度
end
hwso=hW-ho; % 液封深度
%% 塔板排布
tbx=D/2-(Wd+Ws);
tbR=D/2-Wc;
Aa=2*(tbx*sqrt(tbR^2-tbx^2)+tbR^2*asin(tbx/tbR)); %鼓泡区面积
if (get(handles.checkbox_fkt,'value')==1)
    fkt=geteditVar({'fkt'},handles);
    if isempty(fkt)
        return
    end
end
if get(handles.popupmenu_BanXing,'value')==1
    Fo=geteditVar({'Fo'},handles);
    if isempty(Fo)
        return
    end
    uo=Fo/sqrt(rhoV); % 气体通过阀孔时的速度 m/s
    N=floor(Vs/(pi/4*do^2*uo)); % 每层板上的阀孔数
    Ao=Vs/uo; %阀孔总面积 m2
else
    N=floor(2/sqrt(3)*Aa/fkt^2); % 筛孔数
    uo=Vs/(pi/4*do^2*N); %孔气速
    Ao=pi/4*do^2*N; %筛孔总面积 m2
    Fo=uo*sqrt(rhoV); %孔动能因数
end
if (get(handles.popupmenu_FuFaJiHe,'value')==1) &&  ...
        (get(handles.popupmenu_BanXing,'value')==1) && ...
        (get(handles.checkbox_fkt,'value')==0)
    fkt=do*sqrt(pi/4/sin(pi/3)*Aa/Ao); % 叉排 液流向垂直阀孔中心间距 m
end
if (get(handles.popupmenu_BanXing,'value')==1)
    fktp=Aa/N/fkt; % 叉排 液流向平行阀孔中心间距 m
    setappdata(handles.figure_main,'fktp',fktp);
end
%重新排布时
if get(handles.checkbox_fkN,'value')~=0
    fkN=geteditVar({'fkN'},handles);
    if isempty(fkN)
        return
    end
    N=fkN; %阀孔数
    uo=Vs/(pi/4*do^2*N); %孔气速
    Fo=uo*sqrt(rhoV); %阀孔动能因数
    Ao=pi/4*do^2*N; %阀孔总面积 m2
end
portrate=Ao/Aa; % 开孔率
%% 保存变量
cellfun(svvar,{'scH','scLV','umax','u','theta','Aa','N','uo', ...
    'Ao','portrate','AT','ELsmax'}, ...
    {scH,scLV,umax,u,theta,Aa,N,uo,Ao,portrate,AT,ELsmax});
cellfun(svvar,{'C20','lW','lWD','hOW','hW','Wd','Af','hwso', ...
    'Fo','ho','fkt','pzyE','hL'}, ...
    {C20,lW,lWD,hOW,hW,Wd,Af,hwso,Fo,ho,fkt,pzyE,hL});
%% 显示结果
opstr={'工艺尺寸计算成功！'};
set(handles.edit_output,'string',opstr);

function uipanel_tab_SelectionChangeFcn(hObject,~,handles)
slobj=get(hObject,'SelectedObject');
slcl=[83 82 78]/100;
uscl=[1 1 1]*0.94;
if strcmp(get(slobj,'tag'),'radiobutton_size')
    set(handles.uipanel_hydromech,'visible','off');
    set(handles.uipanel_size,'visible','on');
    set(slobj,'BackgroundColor',slcl);
    set(handles.radiobutton_hydro,'BackgroundColor',uscl);
else
    set(handles.uipanel_size,'visible','off');
    set(handles.uipanel_hydromech,'visible','on');
    set(slobj,'BackgroundColor',slcl);
    set(handles.radiobutton_size,'BackgroundColor',uscl);
end

function [stat]=HydrodynamicsOfpF(~, ~, handles)
% 单溢流浮阀塔流体力学验算函数

% 选定参数检查
[upvarepsilon,upphi,hyK]=geteditVar(...
    {'upvarepsilon','upphi','hyK'}, ...
    handles);
if isempty(hyK)
    stat=-1;
    return
end
% 获取工艺尺寸结果
[rhoV,rhoL,uo,sigma,Ls,lW,ho,hL,HT,hW,D,Wd,AT,Af,Vs]=geteditVar( ...
    {'rhoV','rhoL','uo','sigma','Ls','lW','ho','hL','HT','hW','D', ...
    'Wd','AT','Af','Vs'},handles);
if isempty(Vs)
    stat=-1;
    return
end
%% 单溢流浮阀塔
if get(handles.checkbox_g,'value')==0
    g=9.81;   % 重力加速度 m/s
else
    g=geteditVar({'g'},handles);
    if isempty(g)
        stat=-1;
        return
    end
end
if get(handles.checkbox_CF,'value')==0
    CF=vapourcapacityfactor(HT,rhoV);
else
    CF=geteditVar({'CF'},handles);
    if isempty(CF)
        stat=-1;
        return
    end
end
% F1型重阀
uoc=nthroot(19.9*2/5.34*g/rhoV,1.825); % 临界孔速 m/s
% 干板阻力
if uo<=uoc
    hc=19.9*uo^0.175/rhoL; % 阀全开前干板压头 m液柱
else
    hc=5.34*rhoV*uo^2/2/rhoL/g; % 阀全开后干板压头 m液柱
end
hl=upvarepsilon*hL; % 板上充气液层阻力 m液柱
% 液体表面张力所造成的阻力 m液柱
if get(handles.checkbox_flh,'value')==1
    flh=geteditVar({'flh'},handles);
    if isempty(flh)
        stat=-1;
        return
    end
    hsigma=2*sigma*1e-3/flh/rhoL/g;
else
    hsigma=0; % 忽略该阻力
end
hp=hc+hl+hsigma; % 单板压头 m液柱
dltpp=hp*rhoL*g; % 单板压降
dltpc=hc*rhoL*g;
dltpl=hl*rhoL*g;
dltpsigma=hsigma*rhoL*g;
%% 淹塔
% 液体流过降液管压强降液柱当量 m液柱
if get(handles.popupmenu_JinKouYan,'value')==0
    hd=0.153*(Ls/lW/ho)^2;
else
    hd=0.2*(Ls/lW/ho)^2;
end
Hd=hp+hL+hd; % 降液管内清液层高度 m
if Hd > upphi*(HT+hW)
    tem=sprintf(['不符合液泛安全：\n     Hd > upphi*(HT+hW) \n',...
        '     5.4f > 5.4f'],Hd,upphi*(HT+hW));
    dlg_warnvar(tem,'警告!');
end
%%　雾沫夹带
ZL=D-2*Wd; % 板上液体流径长度 m
Ab=AT-2*Af; % 板上液流面积 m2
% 泛点率 1
floodL=(Vs*sqrt(rhoV/(rhoL-rhoV))+1.36*Ls*ZL)/(hyK*CF*Ab);
% 泛点率 2
floodS=(Vs*sqrt(rhoV/(rhoL-rhoV)))/(0.78*hyK*CF*AT);
%% 保存变量
stat=0;
svvarx=@(hh)@(nvar,val)setappdata(hh,nvar,val);
svvar=svvarx(handles.figure_main);
cellfun(svvar,{'uoc','hc','hl','hsigma','hp','dltpp','dltpc', ...
    'dltpl','dltpsigma','hd','Hd','ZL','Ab','floodL','floodS'}, ...
    {uoc,hc,hl,hsigma,hp,dltpp,dltpc,dltpl,dltpsigma,hd,Hd,ZL,Ab, ...
    floodL,floodS});
cellfun(svvar,{'g','CF'},{g,CF});

function [stat]=HydrodynamicsOfpS(~, ~, handles)
% 单溢流筛板塔流体力学验算函数

% 选定参数检查
[upvarepsilon,upphi]=geteditVar(...
    {'upvarepsilon','upphi'}, ...
    handles);
if isempty(upphi)
    stat=-1;
    return
end
% 获取工艺尺寸结果
[rhoV,rhoL,uo,sigma,do,delta,Ls,lW,ho,HT,hW,AT,Af,Ao,Aa,hL,Vs]= ...
    geteditVar({'rhoV','rhoL','uo','sigma','do','delta','Ls','lW', ...
    'ho','HT','hW','AT','Af','Ao','Aa','hL','Vs'}, ...
    handles);
if isempty(Vs)
    stat=-1;
    return
end
%% 单溢流筛板塔
if get(handles.checkbox_g,'value')==0
    g=9.81;   % 重力加速度 m/s
else
    g=geteditVar({'g'},handles);
    if isempty(g)
        stat=-1;
        return
    end
end
% 流量系数
if get(handles.checkbox_co,'value')==0
    co=0.87054163-0.0646663*(do/delta)+0.6756572e-2*(do/delta)^2- ...
        0.2989456e-3*(do/delta)^3;
    if do >= 0.010
        co=co*1.15;
    end
else
    co=geteditVar({'co'},handles);
    if isempty(co)
        stat=-1;
        return
    end
end
hc=0.051*(uo/co)^2*(rhoV/rhoL)*(1-(Ao/Aa)^2); % 干板阻力
hl=upvarepsilon*hL; % 板上充气液层阻力 m液柱
ua=Vs/(AT-Af);
Fo=ua*sqrt(rhoV); % 气相动能因子
hsigma=4*sigma*1e-3/rhoL/g/do; % 表面张力阻力 m液柱
hp=hc+hl+hsigma; % 单板压头 m液柱
dltpp=hp*rhoL*g; % 单板压降
dltpc=hc*rhoL*g;
dltpl=hl*rhoL*g;
dltpsigma=hsigma*rhoL*g;
% % 降液管液泛
% 液体流过降液管压降液柱当量 m液柱
if get(handles.popupmenu_JinKouYan,'value')==1
    hd=0.153*(Ls/lW/ho)^2;
else
    hd=0.2*(Ls/lW/ho)^2;
end
Hd=hp+hL+hd; % 降液管内清液层高度 m
if Hd > upphi*(HT+hW)
    warnHd=sprintf(['不符合液泛安全：Hd > upphi*(HT+hW)',...
        '            %5.4f > %5.4f'],Hd,upphi*(HT+hW));
    dlg_warnvar(warnHd,'警告');
    stat=-1;
    return
end
% %　雾沫夹带
if get(handles.checkbox_hfrL,'value')==0
    hfrL=2.5;
else
    hfrL=geteditVar({'hfrL'},handles);
    if isempty(hfrL)
        stat=-1;
        return
    end
end
hf=hfrL*hL; % 塔板上鼓泡层高度 m
eV=5.7e-3/sigma*(ua/(HT-hf))^3.2; % 液沫夹带量 kg液/kg气
if eV > 0.1
    warneV=sprintf(['液沫夹带量\n\t\t eV=%.5f  ', ...
        '> 0.1  kg液/kg气 ,\n 不在允许范围内'],eV);
    dlg_warnvar(warneV,'警告');
end
if get(handles.checkbox_fhFo,'value')==1
    uomin=fhFo/sqrt(rhoV);%   漏液点气速 m/s
elseif hL < 0.030  || do < 0.003
    uomin=4.4*co*sqrt((0.01+0.13*hL-hsigma)*rhoL/rhoV);
else
    uomin=4.4*co*sqrt((0.0056+0.13*hL-hsigma)*rhoL/rhoV);
end
lyK=uo/uomin; % 稳定系数
%% 保存变量
stat=0;
svvarx=@(hh)@(nvar,val)setappdata(hh,nvar,val);
svvar=svvarx(handles.figure_main);
cellfun(svvar,{'co','hc','hl','hsigma','hp','dltpp','dltpc', ...
    'dltpl','dltpsigma','hd','Hd','hf','eV','uomin','lyK'}, ...
    {co,hc,hl,hsigma,hp,dltpp,dltpc,dltpl,dltpsigma,hd,Hd,hf,eV, ...
    uomin,lyK});
cellfun(svvar,{'g','Fo','hfrL'},{g,Fo,hfrL});

function LoadperformPlot(handles,Ls,Vs,minLs,maxLs,opminLs,opmaxLs, ...
    opminVs,opmaxVs,ymjdVsLf,yfVsf,lyVsf)
% 负荷性能图绘制函数
rgLs=linspace(minLs,maxLs,100);
% rgLs=linspace( 8.200166110939185e-04,8.575051217048211e-04,80);
%     [minLs maxLs],[lyVsf(minLs),lyVsf(maxLs)],'-k', ...
plot(handles.axeslp, ...
    rgLs,ymjdVsLf(rgLs),'-r',rgLs,yfVsf(rgLs),'-b', ...
    [maxLs maxLs],[lyVsf(maxLs),ymjdVsLf(maxLs)],'-g', ...
    rgLs,lyVsf(rgLs),'-k', ...
    [minLs minLs],[lyVsf(minLs) ymjdVsLf(minLs)],'-m', ...
    [0 maxLs],[0 maxLs*Vs/Ls],'-c');
box(handles.axeslp,'on');
legend(handles.axeslp,{' 液沫夹带线','液泛线','液相负荷上限', ...
    '气相负荷下限','液相负荷下限','操作线'},'Location','northwest');
xlabel(handles.axeslp,'{\it L_s}/m^3\cdot s^{-1}');
ylabel(handles.axeslp,'{\it V_s}/m^3\cdot s^{-1}');
hold(handles.axeslp,'on');
plot(handles.axeslp,Ls,Vs,'kp','Markersize',3);
plot(handles.axeslp,[opminLs opmaxLs],[opminVs opmaxVs], ...
    'ko','Markersize',3);
hold(handles.axeslp,'off');

function [stat]=LoadperformOfpF(~, ~, handles)
% 单溢流 浮阀塔 塔板负荷性能 计算绘图函数

% 选定参数检查
[fhhOW,fhtheta,fhFo,hyK,CF,upphi,upvarepsilon]=geteditVar(...
    {'fhhOW','fhtheta','fhFo','hyK','CF','upphi', ...
    'upvarepsilon'},handles);
if isempty(upvarepsilon)
    stat=-1;
    return
end
% 获取工艺尺寸结果
[pzyE,lW,Af,HT,do,N,rhoV,rhoL,ho,Ls,hW,D,Vs]=geteditVar( ...
    {'pzyE','lW','Af','HT','do','N','rhoV', ...
    'rhoL','ho','Ls','hW', 'D','Vs'}, ...
    handles);
if isempty(Vs)
    stat=-1;
    return
end
%获取水力学结果
[ZL,Ab,g,hsigma]=geteditVar(...
    {'ZL','Ab','g','hsigma'}, ...
    handles);
if isempty(Ab)
    stat=-1;
    return
end
%% 负荷性能图
% 可默认参数
if get(handles.checkbox_fhfloodL,'value')==0
    fhfloodL=0.8; % 按 80% 的泛点率算
else
    fhfloodL=geteditVar({'fhfloodL'},handles);
    if isempty(fhfloodL)
        return
    end
end
% 液相负荷下限线
if get(handles.popupmenu_YanXing,'value')==1
    minLs=(fhhOW/pzyE*1000/2.84)^(3/2)*lW/3600;% 单溢流平直堰
    minLsfx=@(fhhOW,D,lW)@(Ls) (fhhOW/weircorrectionfactor(Ls,D,lW)* ...
        1000/2.84)^(3/2)*lW/3600;
    minLsf=minLsfx(fhhOW,D,lW);
    minLs=iterationSim(minLsf,minLs,1e-6);
else
    hn=geteditVar({'hn'},handles);
    if isempty(hn)
        stat=-1;
        return
    end
    if fhhOW <= hn
        minLs=lW/hn*(fhhOW/0.0442)^(5/2)/3600;
    else
        minLs=2464*lW/hn*(fhhOW.^(5/2)-(fhhOW-hn).^(5/2))/3600;
    end
end
% 液相负荷上限线
maxLs=Af*HT/fhtheta;
% 漏液线
minVs=pi/4*do^2*N*fhFo/sqrt(rhoV);% 浮阀
lyVsfx=@(minVs)@(Ls)minVs+Ls.*0; % 转化为匿名函数
lyVsf=lyVsfx(minVs);
% 雾沫夹带线
ymjdVs=@(floodL,hyK,CF,Ab,ZL,rhoV,rhoL)@(Ls) ...
    (floodL*hyK*CF*Ab-1.36*Ls*ZL)/sqrt(rhoV/(rhoL-rhoV));
ymjdVsf=ymjdVs(fhfloodL,hyK,CF,Ab,ZL,rhoV,rhoL);
% 液泛线
% 单溢流 平直堰 齿形堰
yfVs=@(rhoL,g,do,N,rhoV,upphi,HT,hW,hsigma, ...
    upvarepsilon,lW,ho,coehd)@(cxhOWf)@(Ls) ...
    sqrt(2*rhoL*g*(pi/4*do^2*N)^2/5.34/rhoV*( ...
    upphi*(HT+hW)-hsigma-(upvarepsilon+1)* ...
    (hW+cxhOWf(Ls))- coehd*(Ls/lW/ho).^2));
if get(handles.popupmenu_JinKouYan,'value')==1
    yfVsfx=yfVs(rhoL,g,do,N,rhoV,upphi,HT,hW,hsigma, ...
        upvarepsilon,lW,ho,0.153); % 无进口堰
else
    yfVsfx=yfVs(rhoL,g,do,N,rhoV,upphi,HT,hW,hsigma, ...
        upvarepsilon,lW,ho,0.2); % 有进口堰
end
if get(handles.popupmenu_YanXing,'value')==1
    cxhOW1=@(D,lW)@(Ls) ...
        2.84e-3*weircorrectionfactor(Ls,D,lW).*(3600*Ls/lW).^(2/3);
    cxhOWf1=cxhOW1(D,lW);
    yfVsf=yfVsfx(cxhOWf1); % 平直堰
else
    yfVscx=@(lW,hn,yfVsfx)@(Ls)Sawtoothweir(Ls,lW,hn,yfVsfx);
    yfVsf=yfVscx(lW,hn,yfVsfx); % 齿形堰
end
%% 操作线气相负荷 上限 下限
tanVL=Vs/Ls;
ymjdmLmax=fhfloodL*hyK*CF*Ab/ ...
    (tanVL*sqrt(rhoV/(rhoL-rhoV))+1.36*ZL);
yfLmaxfx=@(tanVL,yfVsf)@(Ls) tanVL*Ls-yfVsf(Ls);
yfLmaxf=yfLmaxfx(tanVL,yfVsf);
yfLmax=fzero(yfLmaxf,maxLs/2,optimset('TolX',1e-6)); %操作线与液泛线相交
Ind=zeros(1,2);
[opmaxLs,Ind(2)]=min([ymjdmLmax,yfLmax,maxLs]); %操作线液相负荷 上限
[opminLs,Ind(1)]=max([minVs/tanVL minLs]);%操作线液相负荷 下限
opmaxVs=opmaxLs*tanVL;%操作线气相负荷 上限
opminVs=opminLs*tanVL;%操作线气相负荷 下限
flexib=opmaxLs/opminLs; % 操作弹性
%% 作图
LoadperformPlot(handles,Ls,Vs,minLs,maxLs,opminLs,opmaxLs, ...
    opminVs,opmaxVs,ymjdVsf,yfVsf,lyVsf);
%% 保存结果
stat=0;
svvarx=@(hh)@(nvar,val)setappdata(hh,nvar,val);
svvar=svvarx(handles.figure_main);
cellfun(svvar,{'Ind1','Ind2','minLs','maxLs','opminLs','opmaxLs', ...
    'ymjdVsf','yfVsf','lyVsf'}, ...
    {Ind(1),Ind(2),minLs,maxLs,opminLs,opmaxLs, ...
    ymjdVsf,yfVsf,lyVsf});
cellfun(svvar,{'opmaxVs','opminVs','flexib'}, ...
    {opmaxVs,opminVs,flexib});
%% 显示结果
set(handles.text_opmaxVs,'string',sprintf('%.3f',opmaxVs));
set(handles.text_opminVs,'string',sprintf('%.3f',opminVs));
set(handles.text_flexib,'string',sprintf('%.1f',flexib));

function [stat]=LoadperformOfpS(~, ~,handles)
% 单溢流 浮阀塔 塔板负荷性能 计算及绘图函数

% 选定参数检查
[fhhOW,fhtheta,upphi,upvarepsilon]=geteditVar(...
    {'fhhOW','fhtheta','upphi','upvarepsilon'},handles);
if isempty(upvarepsilon)
    stat=-1;
    return
end
% 获取工艺尺寸结果
[pzyE,lW,Af,HT,do,N,rhoV,hL,Ao,hW,rhoL,ho,Ls,sigma,AT,Aa,D,Vs]= ...
    geteditVar( {'pzyE','lW','Af','HT','do','N','rhoV', ...
    'hL','Ao','hW','rhoL','ho','Ls','sigma','AT','Aa','D','Vs'}, ...
    handles);
if isempty(Vs)
    stat=-1;
    return
end
%获取水力学结果
[co,hsigma,hfrL]=geteditVar({'co','hsigma','hfrL'},handles);
if isempty(hfrL)
    stat=-1;
    return
end
%% 负荷性能图
if get(handles.popupmenu_YanXing,'value')==1
    cxhOW1=@(D,lW)@(Ls) ...
        2.84e-3*weircorrectionfactor(Ls,D,lW).*(3600*Ls/lW).^(2/3);
    cxhOWf1=cxhOW1(D,lW);
else
    hn=geteditVar({'hn'},handles);
    if isempty(hn)
        stat=-1;
        return
    end
    Vsfcx=@(lW,hn,Vsfx)@(Ls)Sawtoothweir(Ls,lW,hn,Vsfx);% 齿形堰
end
% 液相负荷下限线
if get(handles.popupmenu_YanXing,'value')==1
    minLs=(fhhOW/pzyE*1000/2.84)^(3/2)*lW/3600;% 单溢流 平直堰
    minLsfx=@(fhhOW,D,lW)@(Ls) (fhhOW/weircorrectionfactor(Ls,D,lW)* ...
        1000/2.84)^(3/2)*lW/3600;
    minLsf=minLsfx(fhhOW,D,lW);
    minLs=iterationSim(minLsf,minLs,1e-6);
elseif fhhOW <= hn
    minLs=lW/hn*(fhhOW/0.0442)^(5/2)/3600;
else
    minLs=2464*lW/hn*(fhhOW.^(5/2)-(fhhOW-hn).^(5/2))/3600;
end
% 液相负荷上限线
maxLs=Af*HT/fhtheta;
% 漏液线
lyVs=@(co,Ao,hW,hsigma,rhoL,rhoV,coeuo)@(cxhOWf)@(Ls) ...
    4.4*co*Ao*sqrt((coeuo+0.13*(hW+cxhOWf(Ls))-hsigma)*rhoL/rhoV);
if hL < 0.030  || do < 0.003
    lyVsfx=lyVs(co,Ao,hW,hsigma,rhoL,rhoV,0.01);
else
    lyVsfx=lyVs(co,Ao,hW,hsigma,rhoL,rhoV,0.0056);
end
if get(handles.popupmenu_YanXing,'value')==1
    lyVsf=lyVsfx(cxhOWf1);% 平直堰
else
    lyVsf=Vsfcx(lW,hn,lyVsfx); % 齿形堰
end
if  get(handles.checkbox_fhFo,'value')==1
    fhFo=geteditVar({'fhFo'},handles);
    if isempty(fhFo)
        stat=-1;
        return
    end
    minVs=pi/4*do^2*N*fhFo/sqrt(rhoV);
    lyVsfx=@(minVs)@(Ls)minVs;
    lyVsf=lyVsfx(minVs);
end
% 液沫夹带线
if get(handles.checkbox_fheV,'value')==0
    fheV=0.1;
else
    fheV=geteditVar({'fheV'},handles);
    if isempty(fheV)
        stat=-1;
        return
    end
end
ymjdVs=@(fheV,sigma,AT,Af,HT,hfrL,hW)@(cxhOWf)@(Ls) ...
    (fheV*sigma*1e-3/5.7e-6)^(1/3.2)*(AT-Af)* ...
    (HT-hfrL*(hW+cxhOWf(Ls)));
ymjdVsfx=ymjdVs(fheV,sigma,AT,Af,HT,hfrL,hW);
if  get(handles.popupmenu_YanXing,'value')==1
    ymjdVsf=ymjdVsfx(cxhOWf1);% 平直堰
else
    ymjdVsf=Vsfcx(lW,hn,ymjdVsfx); % 齿形堰
end
% 液泛线
% 单溢流 平直堰 齿形堰 匿名函数
yfVs=@(co,rhoL,rhoV,Ao,Aa,upphi,HT,hW,hsigma, ...
    upvarepsilon,lW,ho,coehd)@(cxhOWf)@(Ls) ...
    sqrt(Ao^2*co^2*rhoL/0.051/rhoV/(1-(Ao/Aa)^2)* ...
    (upphi*(HT+hW)-hsigma-(upvarepsilon+1)* ...
    (hW+cxhOWf(Ls))- coehd*(Ls/lW/ho).^2));
if get(handles.popupmenu_JinKouYan,'value')==1
    yfVsfx=yfVs(co,rhoL,rhoV,Ao,Aa,upphi,HT,hW,hsigma, ...
        upvarepsilon,lW,ho,0.153); % 无进口堰
else
    yfVsfx=yfVs(do,co,rhoL,rhoV,Ao,Aa,upphi,HT,hW,hsigma, ...
        upvarepsilon,lW,ho,0.2); % 有进口堰
end
if get(handles.popupmenu_YanXing,'value')==1
    yfVsf=yfVsfx(cxhOWf1);% 平直堰
else
    yfVscx=@(lW,hn,yfVsfx)@(Ls)Sawtoothweir(Ls,lW,hn,yfVsfx);
    yfVsf=yfVscx(lW,hn,yfVsfx); % 齿形堰
end
%% 操作线气相负荷 上限 下限
tanVL=Vs/Ls;
ymjdVmaxfx=@(tanVL,ymjdVsf)@(Ls) ...
    tanVL*Ls-ymjdVsf(Ls);
ymjdLmax=fzero(ymjdVmaxfx(tanVL,ymjdVsf),(minLs+maxLs)/2, ...
    optimset('TolX',1e-6));
yfLmaxfx=@(tanVL,yfVsf)@(Ls) tanVL*Ls-yfVsf(Ls);
yfLmaxf=yfLmaxfx(tanVL,yfVsf);
yfLmax=fzero(yfLmaxf,(minLs+maxLs)/2,optimset('TolX',1e-6));
Ind=zeros(1,2);
[opmaxLs,Ind(2)]=min([ymjdLmax yfLmax maxLs]);%操作线液相负荷 上限
opmaxVs=opmaxLs*tanVL;%操作线气相负荷 上限
lyLminfx=@(tanVL,lyVsf)@(Ls) tanVL*Ls-lyVsf(Ls);
lyLminf=lyLminfx(tanVL,lyVsf);
lyLmin=fzero(lyLminf,(minLs+maxLs)/2,optimset('TolX',1e-6));
[opminLs,Ind(1)]=max([lyLmin minLs]);%操作线气相负荷 下限
opminVs=opminLs*tanVL;%操作线液相负荷 下限
flexib=opmaxLs/opminLs; % 操作弹性
%% 作图
LoadperformPlot(handles,Ls,Vs,minLs,maxLs,opminLs,opmaxLs, ...
    opminVs,opmaxVs,ymjdVsf,yfVsf,lyVsf);
% disp(ymjdVsf);
%% 保存结果
stat=0;
svvarx=@(hh)@(nvar,val)setappdata(hh,nvar,val);
svvar=svvarx(handles.figure_main);
cellfun(svvar,{'Ind1','Ind2','minLs','maxLs','opminLs','opmaxLs', ...
    'ymjdVsf','yfVsf','lyVsf'}, ...
    {Ind(1),Ind(2),minLs,maxLs,opminLs,opmaxLs, ...
    ymjdVsf,yfVsf,lyVsf});
cellfun(svvar,{'opmaxVs','opminVs','flexib','fheV'}, ...
    {opmaxVs,opminVs,flexib,fheV});
%% 显示结果
set(handles.text_opmaxVs,'string',sprintf('%.3f',opmaxVs));
set(handles.text_opminVs,'string',sprintf('%.3f',opminVs));
set(handles.text_flexib,'string',sprintf('%.1f',flexib));

function pushbutton_hydro_calc_Callback(~, ~, handles)
if get(handles.popupmenu_BanXing,'value')==1
    stat=HydrodynamicsOfpF([], [], handles);
else
    stat=HydrodynamicsOfpS([], [], handles);
end
if stat < 0
    return
end
if get(handles.popupmenu_BanXing,'value')==1
    stat=LoadperformOfpF([], [], handles);
else
    stat=LoadperformOfpS([], [], handles);
end
if stat < 0
    return
end
%% 显示结果
opstr=get(handles.edit_output,'string');
opstr=vertcat(opstr,{'流体力学计算及绘图成功！'});
set(handles.edit_output,'string',opstr);

function uimenu_Open_Callback(~, ~, handles)
[fin,pan]=uigetfile('*.sav','选择参数文件');
fid=fopen(fullfile(pan,fin),'r');
if fid < 3
    return % 打开错误返回
end
varn={'Vs','Ls','rhoV','rhoL','sigma', ...
    'saf','HT','hL','lWD','uop','Wc','Ws','Fo','do','delta', ...
    'BanXing','C20','hOW','pzyE','YanXing','hn','ho','hwso', ...
    'FuFaJiHe','fkt','D','hW','fkN','lW','Wd','Af', ...
    'upvarepsilon','JinKouYan','upphi','hyK','CF','fhtheta', ...
    'fhFo','fhhOW','g','co','flh','fhfloodL','fheV','hfrL'};
C=textscan(fid,'%[^:]:\t%s');
fvar=C{1}; % 参数名
strv=C{2}; % 参数值
finh=@(sty,hndl)hndl(strcmp(sty,get(hndl,'style')));
for i=1:1:length(varn)
    hndl=findobj(handles.figure_main,'-regexp','tag', ...
        ['[^_]+_',varn{i},'$']);
    hchk=finh('checkbox',hndl);
    tf=strcmp(varn{i},fvar);
    if any(tf)
        hpop=finh('popupmenu',hndl);
        hedit=finh('edit',hndl);
        if ~isempty(hchk)
            set(hchk,'value',1)
            cb=get(hchk,'callback');
            cb(hchk,[]);
        end
        if ~isempty(hpop)
            set(hpop,'value',str2double(strv{tf}));
            cb=get(hpop,'callback');
            cb(hpop,[]);
        elseif ~isempty(hedit)
            set(hedit,'string',strv{tf});
            cb=get(hedit,'callback');
            cb(hedit,[]);
        end
    elseif ~isempty(hchk)
        set(hchk,'value',0)
        cb=get(hchk,'callback');
        cb(hchk,[]);
    end
end
fclose(fid);

function uimenu_Save_Callback(~, ~, ~)
[fin,pan]=uiputfile('*.sav','保存输入参数', datestr(now,'yymmmdd'));
fid=fopen(fullfile(pan,fin),'w');
if fid < 3
    return % 打开错误返回
end
varn={'Vs','Ls','rhoV','rhoL','sigma', ...
    'saf','HT','hL','lWD','uop','Wc','Ws','Fo','do','delta', ...
    'BanXing','C20','hOW','pzyE','YanXing','hn','ho','hwso', ...
    'FuFaJiHe','fkt','D','hW','fkN','lW','Wd','Af', ...
    'upvarepsilon','JinKouYan','upphi','hyK','CF','fhtheta', ...
    'fhFo','fhhOW','g','co','flh','fhfloodL','fheV','hfrL'};
for i=1:1:length(varn)
    hvare=findobj(gcf,'-regexp','tag',['[^_]+_',varn{i},'$'], ...
        'style','edit','enable','on');
    hvarp=findobj(gcf,'-regexp','tag',['[^_]+_',varn{i},'$'], ...
        'style','popupmenu','enable','on');
    if ~isempty(hvare)
        fprintf(fid,'%s:\t%s\n',varn{i},get(hvare,'string'));
    elseif ~isempty(hvarp)
        fprintf(fid,'%s:\t%d\n',varn{i},get(hvarp,'value'));
    end
end
fclose(fid);

function uimenu_Export_Callback(~, ~, handles)
% 结果导出回调函数
varn={'Vs','Ls','rhoV','rhoL','sigma', ...
    'saf','HT','hL','lWD','uop','Wc','Ws','Fo','do','delta', ...
    'BanXing','C20','hOW','pzyE','YanXing','hn','ho','hwso', ...
    'FuFaJiHe','fkt','D','hW','fkN','lW','Wd','Af', ...
    'upvarepsilon','JinKouYan','upphi','hyK','CF','fhtheta', ...
    'fhFo','fhhOW','g','co','flh','fhfloodL','fheV','hfrL'};
dpvarn={'scH','scLV','umax','u','theta','Aa','N','uo','Ao', ...
    'portrate','AT','hc','hl','hsigma','hp','dltpp','dltpc','dltpl', ...
    'dltpsigma','hd','Hd','Ind1','Ind2','minLs','maxLs','opminLs', ...
    'opmaxLs','opminVs','opmaxVs','ELsmax','flexib'};
dpFvn={'uoc','ZL','Ab','floodL','floodS','fktp'};
dpSvn={'co','hf', 'eV','uomin','lyK'};
dscdpvn={'史密斯关联图的沉降高度, m','史密斯关联图的液气动能比, 1', ...
    '泛点气速, m/s','空塔气速, m/s','液体在降液管内停留时间, s', ...
    '塔板鼓泡区面积, m^2','一层塔板上的浮阀或筛孔总数', ...
    '阀孔或筛孔气速, m/s','阀孔或筛孔总面积, m^2','开孔率, 1', ...
    '塔截面积, m^2','干板压降的液柱当量高度, m液柱', ...
    '板上液层阻力的液柱当量高度, m液柱', ...
    '克服表面张力的液柱当量高度, m液柱', ...
    '单板压降液柱当量高度, m液柱','单板压降, Pa','干板阻力压降, Pa', ...
    '板上液层阻力压降, Pa','克服液体表面张力压降, Pa', ...
    '液体经过降液管阻力的液柱当量高度, m液柱', ...
    '降液管清液层高度, m',['气相负荷下限控制数, 1:由漏液控制,', ...
    ' 2:由液相负荷下限控制'],['气相负荷上限控制数, 1:由液沫夹带控制,' ...
    ' 2:由液泛控制, 3:由液相负荷上限控制'],'液相负荷下限线常数值, m^3/s', ...
    '液相负荷上限线常数值, m^3/s','塔板的液相负荷下限, m^3/s', ...
    '塔板的液相负荷上限, m^3/s','塔板的气相负荷下限, m^3/s', ...
    '塔板的气相负荷上限, m^3/s','液流收缩系数推荐线流量, m^3/s', ...
    '操作弹性'};
dscdpF={'临界孔速, m/s','板上液体流径长度, m','板上液流面积, m^2', ...
    '泛点率, 1, 由含有流径长度参数的经验公式算出', ...
    '泛点率, 1, 由不含流径长度参数的经验公式算出', ...
    '液流向平行阀孔中心间距, m'};
dscdpS={'流量系数, 无因次','板上鼓泡层高度, m', ...
    '液沫夹带量, kg(液)/kg(气)','漏液点气速, m/s','稳定系数, 无因次'};
%% 变量分类
varn=sort(varn);
if get(handles.popupmenu_BanXing,'value')==1
    dpvn=horzcat(dpvarn,dpFvn);
    dscdp=horzcat(dscdpvn,dscdpF);
else
    dpvn=horzcat(dpvarn,dpSvn);
    dscdp=horzcat(dscdpvn,dscdpS);
end
[dpvn,idx]=sort(dpvn);
dscdp=dscdp(idx);
%% 输出到文件
[fin,pan]=uiputfile('*.sav','保存输入参数',  ...
    ['Exp-',datestr(now,'yymmmdd')]);
fid=fopen(fullfile(pan,fin),'w');
if fid < 3
    return % 打开错误返回
end
fprintf(fid,'%%可定义参数：带#号为自变量\r\n');
exportAllVarToFile(handles,fid,varn,cell(size(varn)),[]);
fprintf(fid,'%%其它因变量\r\n');
exportAllVarToFile(handles,fid,dpvn,dscdp,'dependent');
fclose(fid);

function exportAllVarToFile(handles,fid,vn,dsc,stat)
%% 参数输出整理函数
for i=1:1:length(vn)
    if strcmp(stat,'dependent')
        val=getappdata(handles.figure_main,vn{i});
        str=sprintf('%.3e',val);
    else
        hpop=findobj(handles.figure_main,'tag',['popupmenu_',vn{i}]);
        hed=findobj(handles.figure_main,'tag',['edit_',vn{i}]);
        if ~isempty(hpop)
            val=get(hpop,'value');
            strc=get(hpop,'string');
            str=strc{val};
        elseif ~isempty(hed) && strcmp(get(hed,'enable'),'on')
            str=get(hed,'string');
            vn{i}=horzcat('#',vn{i});
        else
            val=getappdata(handles.figure_main,vn{i});
            str=sprintf('%.3e',val);
        end
    end
    if ~isempty(dsc{i})
        dsc{i}=sprintf('\t%s\r\n',dsc{i});
    end
    fprintf(fid,'%s:\t%s\r\n%s',vn{i},str,dsc{i});
end

function uimenu_ExportAxesFunc_Callback(~, ~, handles)
% 导出图片到另一窗口，并导出函数句柄到工作空间

%获取参数
[Ls,Vs,minLs,maxLs,opminLs,opmaxLs,opminVs,opmaxVs, ...
    ymjdVsf,yfVsf,lyVsf]=geteditVar({'Ls','Vs','minLs','maxLs', ...
    'opminLs','opmaxLs','opminVs','opmaxVs', ...
    'ymjdVsf','yfVsf','lyVsf'},handles);
if isempty(lyVsf)
    return
end
%% 作图到新窗口并导出参数到工作空间
nf.hf=figure();
nf.axeslp=axes('Parent',nf.hf);
LoadperformPlot(nf,Ls,Vs,minLs,maxLs,opminLs,opmaxLs, ...
    opminVs,opmaxVs,ymjdVsf,yfVsf,lyVsf);
axc=cell2struct({Ls,Vs,minLs,maxLs,opminLs,opmaxLs,opminVs,opmaxVs, ...
    ymjdVsf,yfVsf,lyVsf},{'Ls','Vs','minLs','maxLs','opminLs', ...
    'opmaxLs','opminVs','opmaxVs','ymjdVsf','yfVsf','lyVsf'},2);
assignin('base','axb',axc);

function [ a ] = NewtonIteration( a,epsilon ,ff,ffp)
% Use Newton Iteration to solve the equation
Delt = ff(a)/ffp(a);
a = a- Delt;
if abs(Delt) > epsilon
    [a] = NewtonIteration(a,epsilon,ff,ffp);
else
    a = a- Delt;
end

function Eii=iterationSim(func,Ei,precision)
Eii=func(Ei);
if abs(Eii-Ei) > precision
    Eii=iterationSim(func,Eii,precision);
end