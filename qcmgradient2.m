function varargout = qcmgradient2(varargin)
% QCMGRADIENT2 MATLAB code for qcmgradient2.fig
%      QCMGRADIENT2, by itself, creates a new QCMGRADIENT2 or raises the existing
%      singleton*.
%
%      H = QCMGRADIENT2 returns the handle to a new QCMGRADIENT2 or the handle to
%      the existing singleton*.
%
%      QCMGRADIENT2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QCMGRADIENT2.M with the given input arguments.
%
%      QCMGRADIENT2('Property','Value',...) creates a new QCMGRADIENT2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before qcmgradient2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to qcmgradient2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qcmgradient2

% Last Modified by GUIDE v2.5 18-Jul-2016 13:24:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @qcmgradient2_OpeningFcn, ...
    'gui_OutputFcn',  @qcmgradient2_OutputFcn, ...
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

% --- Executes just before qcmgradient2 is made visible.
function qcmgradient2_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(hObject,'toolbar','figure');
handles.zq=8.84e6;
handles.f1=5e6;
set(hObject,'CloseRequestFcn',@closeGUI)
set(0,'defaultlinelinewidth',2)
addpath('SLMtools')
guidata(hObject, handles);

function varargout = qcmgradient2_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function init_Callback(hObject, eventdata, handles)
% rheological parameters are defined and stored in the handles structure.
% Generally we should not need to do this.
handles.rheology=getrheology(hObject, handles);
plotgradient(hObject,handles)
guidata(hObject, handles)

function ming_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function ming_CreateFcn(hObject, eventdata, handles)

function maxg_Callback(hObject, eventdata, handles)
% we need to set maximum possible value of g to G0
[handles, rheology] = checkrheology(hObject, handles)
[nh, legendtext]=getharmonics(hObject, handles);
if str2num(get(handles.maxg,'string')) > max(rheology.g(max(nh),:))
    str=sprintf('%0.3e', max(rheology.g(max(nh),:)));
%     str=[str(1:4) str(6) str(10)];
    set(handles.maxg,'string',str);
end
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function maxg_CreateFcn(hObject, eventdata, handles)

function nlayers_Callback(hObject, eventdata, handles)
nl = get(handles.nlayers, 'value');
for i = 1:10
    if i <= nl
        set(handles.(['slider' num2str(i)]), 'visible', 'on')
    else
       set(handles.(['slider' num2str(i)]), 'visible', 'off') 
    end
end
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function nlayers_CreateFcn(hObject, eventdata, handles)

function thickness_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function thickness_CreateFcn(hObject, eventdata, handles)

function rho_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
        handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles) 
    guidata(hObject, handles)
end
function rho_CreateFcn(hObject, eventdata, handles)

function g0_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles) 
    guidata(hObject, handles)
end
function g0_CreateFcn(hObject, eventdata, handles)

function gr_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
        handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles) 
    guidata(hObject, handles)
end
function gr_CreateFcn(hObject, eventdata, handles)

function x1dataselect_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function x1dataselect_CreateFcn(hObject, eventdata, handles)

function y1dataselect_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function y1dataselect_CreateFcn(hObject, eventdata, handles)

function x2dataselect_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function x2dataselect_CreateFcn(hObject, eventdata, handles)

function y2dataselect_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function y2dataselect_CreateFcn(hObject, eventdata, handles)

function slider1_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider1_CreateFcn(hObject, eventdata, handles)

function slider2_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider2_CreateFcn(hObject, eventdata, handles)

function slider3_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider3_CreateFcn(hObject, eventdata, handles)

function slider4_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider4_CreateFcn(hObject, eventdata, handles)

function slider5_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider5_CreateFcn(hObject, eventdata, handles)

function slider6_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider6_CreateFcn(hObject, eventdata, handles)

function slider7_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider7_CreateFcn(hObject, eventdata, handles)

function slider8_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider8_CreateFcn(hObject, eventdata, handles)

function slider9_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider9_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>

function slider10_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    plotgradient(hObject,handles)
end
function slider10_CreateFcn(hObject, eventdata, handles)

function delg1_CreateFcn(hObject, eventdata, handles)
function delf1_CreateFcn(hObject, eventdata, handles)

function closeGUI(src,~)
%src is the handle of the object generating the callback (the source of the event)
%evnt is the The event data structure (can be empty for some callbacks)
selection = questdlg('Do you want to save the file?',...
    'Close Request Function',...
    'Yes','No','Cancel','Yes');
switch selection,
    case 'Yes',
        hgsave(src,'qcmgradient2.fig')
        delete(gcf)
    case 'No'
        delete(gcf)
    case 'Cancel'
        return
end

function handles = calc_Callback(hObject, eventdata, handles)
plotgradient(hObject, handles);
handles = calcqcm(hObject, handles);
handles = findsolution(hObject, eventdata, handles);
set(handles.herrorresults, 'string', '')

function handles = calcqcm(hObject, handles)
%	Function to predict the QCM response to a multilayer system of curing
%   layer with a gradient of properties
f1=handles.f1;
zq=handles.zq;
handles = checkrheology(hObject, handles);
nh=handles.rheology.nh;
layerprops=getlayerprops(hObject, handles);
if isempty(layerprops)
    warndlg('No layer data entered. Please enter data')
    return
end
delfstar=calcdelfstar(f1,zq,layerprops,nh);
if get(handles.nh1,'value')
    set(handles.delf1,'string',commanumber(real(delfstar(1))));
    set(handles.delg1,'string',commanumber(imag(delfstar(1))));
else
    set(handles.delf1,'string','-');
    set(handles.delg1,'string','-');
end
if get(handles.nh3,'value')
    set(handles.delf3,'string',commanumber(real(delfstar(3))));
    set(handles.delg3,'string',commanumber(imag(delfstar(3))));
else
    set(handles.delf3,'string','-');
    set(handles.delg3,'string','-');
end
if get(handles.nh5,'value')
    set(handles.delf5,'string',commanumber(real(delfstar(5))));
    set(handles.delg5,'string',commanumber(imag(delfstar(5))));
else
    set(handles.delf5,'string','-');
    set(handles.delg5,'string','-');
end

% add point to property plot
function delfstar=calcdelfstar(f1,zq,layerprops,nhvals)
struct2var(layerprops)
rho=rho*1000;  % convert to kg/m^3
for nh=nhvals
    fn=nh*f1;
    for j=1:nlayers
        gstar=glayer(nh,j)*exp(pi*1i*philayer(nh,j)/180);   % complex shear modulus of the layer of the iteration
        zstar=(rho*gstar)^0.5; % complex impedance of the layer of the iteration
        kstar=2*pi*fn*rho/zstar; % Complex wave number of the layer of the iteration
        
        % lstar is propagator of the shear wave through layer j
        % sstar is propagator of the shear wave through the interface of layer j-1 and j
        lstar=[exp(1i*kstar*subthickness(j)),0;...
            0,exp(-1i*kstar*subthickness(j))];
        % No sstar is needed for first layer, but lstar and sstar are needed for all other layers
        if j==1
            waveMatrix=lstar; % Start the waveMatrix equation which records the amplitude of the shear wave to and away from the surface
            z1star=zstar; % Save the compl impedance of first layer because we need it later
            zm1star=zstar;
        else
            sstar=0.5*[1+zstar/zm1star,1-zstar/zm1star;1-zstar/zm1star,1+zstar/zm1star];
            waveMatrix=waveMatrix*sstar*lstar;
            zm1star=zstar;   % Set the newest value of zstar to the minus 1 for the next loop
        end
    end
    % Termination matrix
    %This is here for generality, and is not really needed for the free surface
    %boundary condition that we are using
    
    tstar=[1,0;0,1];
    waveMatrix=waveMatrix*tstar*[1;1];
    
    % Extract u values and determine overall change in impedance
    
    uPosstar=waveMatrix(1,1);   % Amplitude of the wave traveling toward the QCM surface
    uNegstar=waveMatrix(2,1);   % Amplitude of the wave traveling away from the QCM surface
    rstar=uNegstar/uPosstar;   % Reflection coefficient describing the amplitude ratio for the forward and reflected waves
    zContstar=z1star*(1-rstar)/(1+rstar);   % Total complex impedance of the loading material
    delfstar(nh)=f1*1i*zContstar/pi/zq; % Calculate the complex frequency shift from the complex impedance zContstar
end

function handles = update_Callback(hObject, eventdata, handles) %#ok<*INUSL>
plotgradient(hObject, handles);
handles = findsolution(hObject, eventdata, handles);

function plotgradient(hObject,handles)
[handles, rheology] = checkrheology(hObject, handles);
layerprops=getlayerprops(hObject,handles);
if isempty(layerprops)
    warndlg('No layer data entered. Please enter data')
    return
end
nlayers=layerprops.nlayers;
d=layerprops.d;
subthickness=layerprops.subthickness;
% get the x and y data based on the dropdown menus
[yplotdata]=getaxesdata(handles.y1dataselect,rheology,layerprops);
ylayerdata=yplotdata.layerdata;

for sublayer=1:nlayers
    xdata(2*sublayer-1)=sum(subthickness(1:sublayer-1));
    xdata(2*sublayer)=sum(subthickness(1:sublayer));
%   xdata(2*sublayer-1)=(sublayer-1)*subthickness(sublayer);
%   xdata(2*sublayer)=sublayer*subthickness(sublayer);
    for nh=rheology.nh
        ydata(nh,2*sublayer-1:2*sublayer)=ylayerdata(nh,sublayer);
    end
end
xdata(2*nlayers+1)=d;
ydata(:,2*nlayers+1)=0;

switch get(handles.x1dataselect,'value')
    case 1
        xdata=1e6*xdata;  % convert to microns
        xlabeltext='z (\mum)';
    case 2
        xdata=xdata/d;
        xlabeltext='z/d';
end

% now make the plot
for nh=rheology.nh
    if nh==min(rheology.nh)
        hold(handles.axes1,'off')
    else
        hold(handles.axes1,'on')
    end
    plot(handles.axes1,xdata,ydata(nh,:),...
        'color',rheology.colors{nh},'linestyle','-');
end

% add plot labels
xlabel(handles.axes1,xlabeltext);
ylabel(handles.axes1,yplotdata.labeltext);

% handle linear/log scaling of y axis
if get(handles.y1auto,'value')
    set(handles.axes1,'yscale',yplotdata.scale)
elseif get(handles.y1lin,'value')
    set(handles.axes1,'yscale','linear')
elseif get(handles.y1log,'value')
    set(handles.axes1,'yscale','log')
end

set(handles.axes1,'xscale','linear')
plotproperties(hObject, handles);

function plotproperties(hObject, handles)
rheology=handles.rheology;
layerprops=getlayerprops(hObject,handles);
if isempty(layerprops)
    warndlg('No layer data entered. Please enter data')
    return
end
% get the x and y data based on the dropdown menus
[xplotdata]=getaxesdata(handles.x2dataselect,rheology,layerprops);
xdata=xplotdata.data;
xlayerdata=xplotdata.layerdata;
[yplotdata]=getaxesdata(handles.y2dataselect,rheology,layerprops);
ydata=yplotdata.data;
ylayerdata=yplotdata.layerdata;

% now make the full plot
[r c]=size(ydata);
for nh=rheology.nh
    if nh==min(rheology.nh)
        hold(handles.axes2,'off')
    else
        hold(handles.axes2,'on')
    end
    plotprop(nh)=plot(handles.axes2,xdata(nh,:),ydata(nh,:), 'color',rheology.colors{nh});
    hold(handles.axes2,'on')
    plotlayer(nh)=plot(handles.axes2,xlayerdata(nh,:),ylayerdata(nh,:),...
        'color',rheology.colors{nh},'marker','+','linestyle','none');
    if r==10;  % these means we have two sets of y data - only happens for G', G"
        plotprop2(nh)=plot(handles.axes2,xdata(nh,:),ydata(nh+5,:),...
	  'color', rheology.colors{nh}, 'linestyle', '-');
        hold(handles.axes2,'on')
        plotlayer2(nh)=plot(handles.axes2,xlayerdata(nh,:),ylayerdata(nh+5,:),...
	  'color',rheology.colors{nh},'marker','o','linestyle','none');
    end
end

% now we add a point for the gradient calculation, with the x axis
% corresponding to the properties of layer 1
switch get(handles.y2dataselect,'value')
    case 6
        delfstar=calcdelfstar(rheology.f1,rheology.zq,layerprops,rheology.nh);
        for nh=rheology.nh
            plot(handles.axes2,xlayerdata(nh,1),real(delfstar(1,nh)),...
                'color',rheology.colors{nh},'marker','o','linestyle','none',...
                'markerfacecolor',rheology.colors{nh});
        end
    case 7
        delfstar=calcdelfstar(rheology.f1,rheology.zq,layerprops,rheology.nh);
        for nh=rheology.nh
            plot(handles.axes2,xlayerdata(nh,1),imag(delfstar(1,nh)),...
                'color',rheology.colors{nh},'marker','o','linestyle','none',...
                'markerfacecolor',rheology.colors{nh});
        end
end

% now we reset the axis properties as desired
if get(handles.uselayerlimits,'value')
    % in this case the limits are determined from the getaxislimits
    % function, and are based on Min G. and Max. G
    xlim(handles.axes2,xplotdata.layerlimits)
    ylim(handles.axes2,'auto')
else
    xlim(handles.axes2,xplotdata.limits)
    ylim(handles.axes2,yplotdata.limits)
end
xlabel(handles.axes2,xplotdata.labeltext);
ylabel(handles.axes2,yplotdata.labeltext);

% now we sort out the log or linear scaling
if get(handles.x2auto,'value')
    set(handles.axes2,'xscale',xplotdata.scale)
elseif get(handles.x2lin,'value')
    set(handles.axes2,'xscale','linear')
elseif get(handles.x2log,'value')
    set(handles.axes2,'xscale','log')
end

if get(handles.y2auto,'value')
    set(handles.axes2,'yscale',yplotdata.scale)
elseif get(handles.y2lin,'value')
    set(handles.axes2,'yscale','linear')
elseif get(handles.y2log,'value')
    set(handles.axes2,'yscale','log')
end

if strcmp(yplotdata.legendtext,'nothing')==0
    legend(plotprop(rheology.nh),rheology.legendtext,'location',yplotdata.legendlocation);
end
set(handles.axes2,'nextplot','replacechildren')
set(handles.axes2,'nextplot','replacechildren')
%now update the calculated parameters for the multilayer model
calcqcm(hObject, handles);

% now we invoke the cursor functionality

function [plotdata]= getaxesdata(dropdownhandle,rheology,layerprops)
struct2var(rheology)
struct2var(layerprops)
switch get(dropdownhandle,'value')
    case 1  % magnigude of g
        plotdata.data=g;
        plotdata.layerdata=glayer;
        plotdata.labeltext='|G*| (Pa)';
        plotdata.scale='log';
        plotdata.limits=[1e5 1e10];
        plotdata.layerlimits=layerprops.glimits;
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.legendlocation='southeast';
    case 2  % phase angle
        plotdata.data=phi;
        plotdata.layerdata=philayer;
        plotdata.labeltext='\phi_{n} (deg)';
        plotdata.scale='linear';
        plotdata.limits='auto';
        plotdata.layerlimits=layerprops.philimits;
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.legendlocation='southwest';
    case 3  % storage and loss moduli
        plotdata.data=[storage; loss];
        plotdata.layerdata=[storagelayer; losslayer];
        plotdata.labeltext='G'', G" (Pa)';
        plotdata.scale='log';
        plotdata.limits=[1e4, 1e10];
        plotdata.layerlimits=layerprops.glimits;
        plotdata.legendtext={'G'''; 'G"'};
        plotdata.legendlocation='southwest';
    case 4
        plotdata.data=1e6*lambda; %#ok<*NODEF>
        plotdata.layerdata=1e6*lambdalayer;
        plotdata.labeltext='\lambda (\mum)';
        plotdata.scale='log';
        plotdata.limits='auto';
        plotdata.layerlimits=layerprops.lambdalimits;
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.legendlocation='southwest';
    case 5
        plotdata.data=1e6*delta;
        plotdata.layerdata=1e6*deltalayer;
        plotdata.labeltext='\delta (\mum)';
        plotdata.scale='log';
        plotdata.limits='auto';
        plotdata.layerlimits=layerprops.deltalimits;
        plotdata.legendtext='nothing';
        plotdata.legendlocation='southwest';
    case 6
        plotdata.data=delf;
        plotdata.layerdata=delflayer;
        plotdata.labeltext='\Deltaf (Hz)';
        plotdata.scale='linear';
        plotdata.limits='auto';
        plotdata.layerlimits=layerprops.delflimits;
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.legendlocation='southwest';
    case 7
        plotdata.data=delg;
        plotdata.layerdata=delglayer;
        plotdata.labeltext='\Delta\Gamma (Hz)';
        plotdata.scale='linear';
        plotdata.limits='auto';
        plotdata.layerlimits=layerprops.delglimits;
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.legendlocation='northwest';
    case 8
        plotdata.data=lambda;
        plotdata.layerdata=lambdalayer;
        plotdata.labeltext='d/\lambda';
        plotdata.scale='linear';
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.limits=[0 1];
        plotdata.layerlimits=d./layerprops.lambdalimits;
        plotdata.legendlocation='southeast';
    case 9
        % frequency data is adjusted so that the first layer gives the
        % plotted value of g* at n=1
        reflayer=1; % frequency = f1 for properties corresponding to reflayer
        plotdata.data=f1*taushift./taushiftlayer(1,1);
        plotdata.layerdata=f1*taushiftlayer/taushiftlayer(1,1);
        plotdata.scale='log';
        plotdata.limits='auto';
        plotdata.layerlimits=f1*layerprops.taureflimits/taushiftlayer(1,1);
        plotdata.labeltext='fa_{T}';
        plotdata.legendtext={'5 MHz'; '15 MHz'};
        plotdata.legendlocation='southwest';
end

function space_Callback(hObject, eventdata, handles)
% adjusts sliders so that layers are evenly spaced
nlayers=get(handles.nlayers,'value');
set(handles.slider1,'value',0)
for nl=2:nlayers
    mastername=['master' num2str(nl)];
    set(handles.(mastername),'value',1)
    slidername=['slider' num2str(nl)];
    set(handles.(slidername)','value',(nl-1)/(nlayers-1)) 
end
update_Callback(hObject, eventdata, handles)

function layerprops=getlayerprops(hObject, handles)
refg = 3;

ming=str2num(get(handles.ming,'string'));
maxg=str2num(get(handles.maxg,'string'));

% all the data area actually generated in the getrheology function, which
% can be modified as needed
rheology=handles.rheology;

%maxg should actually be for the 5th harmonic.
% define limits that are only defined for first harmonic
layerprops.glimits=[.5*ming, maxg*2];
layerprops.taureflimits=getlayerlimits(rheology,'taushift',layerprops.glimits);
layerprops.philimits=getlayerlimits(rheology,'phi',layerprops.glimits);
layerprops.losslimits=getlayerlimits(rheology,'loss',layerprops.glimits);
layerprops.lambdalimits=getlayerlimits(rheology,'lambda',layerprops.glimits);
layerprops.deltalimits=getlayerlimits(rheology,'delta',layerprops.glimits);
layerprops.delflimits=getlayerlimits(rheology,'delf',layerprops.glimits);
layerprops.delglimits=getlayerlimits(rheology,'delg',layerprops.glimits);
layerprops.philimits=getlayerlimits(rheology,'phi',layerprops.glimits);

if get(handles.modelparameters, 'value')
    try
        layerrepeats = str2num(get(handles.repeats, 'string'));
        nlayers=get(handles.nlayers, 'value') * layerrepeats;
    catch
        nlayers=get(handles.nlayers, 'value');
        layerrepeats = 1;
    end
else
    layerdata = handles.layerpropertytable.Data;
    nlayers = find(cellfun('isempty', layerdata),1)-1; %This finds the first empty cell, so the number of layers is one less than that.   
    if isempty(nlayers) %If all of the layers are full by some chance...
        nlayers = max(size(layerdata)); %This should be however many rows there are in the table
    end
    layerrepeats = 1;
    if nlayers == 0
        layerprops = [];
        return
    end
end

%add other layer-related properties to the layerprops structure
layerprops.rho=str2num(get(handles.rho,'string'));
layerprops.nlayers=nlayers;
layerprops.d=1e-6*str2num(get(handles.thickness,'string'));

for nlr=1:nlayers
    if get(handles.modelparameters, 'value')
        nl = mod(nlr, get(handles.nlayers, 'value'));
        if nl == 0
            nl = get(handles.nlayers, 'value');
        end
        % use sliderpositions to establish |g*| for different layers
        % note that this is the value of |g*| at the third harmonic
        slidername=['slider' num2str(nl)];
        relvalue(nlr)=get(handles.(slidername),'Value');
    else
        relvalue(nlr) = str2num(layerdata{nlr,2});
        nl = nlr;
    end
    
    if get(handles.linearlayer,'value')
        gbase(nlr)=ming+relvalue(nl)*(maxg-ming);
    elseif get(handles.loglayer,'value')
        loggbase(nlr)=log(ming)+relvalue(nl)*(log(maxg)-log(ming));
        gbase(nlr)=exp(loggbase(nlr));
    end
    
    [~, uidx] = unique(rheology.g(refg,:)); %The DMA data ended up with some modulus duplicates, which this removes.
    % first figure out the proper frequency to evaluate everyting at
    
    tauref(nlr)=interp1(rheology.g(refg,uidx),rheology.taushift(refg,uidx),gbase(nl));
    if isnan(tauref(nlr))
        if gbase(nlr) > max(rheology.g(refg,:))
            set(handles.maxg, 'string', num2str(max(rheology.g(refg,:)), '%10.2e'));
            gbase(nlr) = max(rheology.g(refg,:));
            tauref(nlr)=interp1(rheology.g(refg,uidx),rheology.taushift(refg,uidx),gbase(nl));
        else
            tauref(nlr) = rheology.taushift(refg, find(round(rheology.g(refg,uidx)) == round(gbase(nl))));
        end
    end
    for nh=rheology.nh
        % now interpolate the desired quantities to this frequency
        layerprops.glayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.g(nh,:),tauref(nl));
        layerprops.philayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.phi(nh,:),tauref(nl));
        layerprops.storagelayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.storage(nh,:),tauref(nl));
        layerprops.losslayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.loss(nh,:),tauref(nl));
        layerprops.lambdalayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.lambda(nh,:),tauref(nl));
        layerprops.deltalayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.delta(nh,:),tauref(nl));
        layerprops.delflayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.delf(nh,:),tauref(nl));
        layerprops.delglayer(nh,nlr)=interp1(rheology.taushift(refg,:),rheology.delg(nh,:),tauref(nl));
        layerprops.taushiftlayer(nh,nlr)=nh/refg*tauref(nlr);
    end

    if get(handles.modelparameters, 'value')
        layerprops.subthickness(nlr)=layerprops.d/layerprops.nlayers;
    else
        layerprops.subthickness(nlr) = 1e-6*str2num(layerdata{nlr,1});
    end
end

function layerlimits=getlayerlimits(rheology, property, glimits)
[~, uidx] = unique(rheology.g(1,:)); %The DMA data ended up with some modulus duplicates, which this removes.
for nh=rheology.nh
    limits=interp1(rheology.g(1,uidx),rheology.(property)(nh,uidx), glimits);
    limits = sort(limits);
    if nh==1
        layerlimits=limits;
    else    
        layerlimits=[min(limits(1),layerlimits(1)), max(limits(2), layerlimits(2))];
        layerlimits=sort(layerlimits);
    end
end
if isnan(layerlimits(1)) || isnan(layerlimits(2))
    if isnan(layerlimits(1))
        layerlimits(1) = -Inf;
    end
    if isnan(layerlimits(2))
        layerlimits(2) = Inf;
    end
end

function rheology=getrheology(hObject, handles)
% determine harmonics to consider
[rheology.nh,rheology.legendtext]=getharmonics(hObject, handles);

%set different plot colors to use
rheology.colors(1)={[1 0 0]};  % red
rheology.colors(3)={[0 0.5 0]}; % green
rheology.colors(5)={[0 0 1]}; % blue
    
%   Define rheological properties of the material
rho=1000*str2num(get(handles.rho,'string')); %#ok<*ST2NM>
d=1e-6*str2num(get(handles.thickness,'string'));
% these get defined in the openingfunction
zq=handles.zq;
f1=handles.f1;

rheology.f1=f1;
rheology.zq=zq;

if get(handles.rheologysimon, 'value')
    npoints=300;
    logmin=-16;
    logmax=0;
    
    taushift=logspace(logmin,logmax,npoints); %this is shift in time due to curing
    % the specific values here come from Simon's paper
    gi=[0.0215 0.0215 0.0215 0.0215 0.0267 0.0267 0.0375 0.0405 0.0630 0.0630 ...
        0.1054 0.1160 0.1160 0.1653 0.0561 0.0561 0.0199 0.0119 0.0055 0.0028 0.0008 ...
        0.0002 0.0003 0.0003];
    taui=logspace(-7, 5, 24);
    gg=str2num(get(handles.g0,'string')); % unrelaxed modulus
    gr=str2num(get(handles.gr,'string'));  % relaxed modulus
    
    OMEGA=@(n,shift) 2*pi*n*f1*shift.*taui; % here OMEGA=2*pi*n*f1*tau
    storage=@(n,shift) gr+(gg-gr)*sum(gi.*OMEGA(n,shift).^2./(1+OMEGA(n,shift)).^2);
    loss=@(n,shift) (gg-gr)*sum(gi.*OMEGA(n,shift)./(1+OMEGA(n,shift)).^2);
    phi=@(n,shift) atan(loss(n,shift)/storage(n,shift));
    g=@(n,shift) (storage(n,shift)^2+(loss(n,shift))^2).^0.5;
    lambda=@(n,shift) (1/(n*f1))*(g(n,shift)/rho)^0.5*(1/cos(phi(n,shift)/2));
    delta=@(n,shift) lambda(n,shift)*cot(phi(n,shift)/2)/(2*pi);
    fs=@(n) -2*n*f1^2*rho*d/zq;
    delfstar=@(n,shift) fs(n)*tan((2*pi*d/lambda(n,shift))*(1-1i*tan(phi(n,shift)/2)))/...
        ((2*pi*d/lambda(n,shift))*(1-1i*tan(phi(n,shift)/2)));
    
    
    % now generate values to plot
    for nh=rheology.nh
        for shiftindex=1:npoints
            rheology.phi(nh,shiftindex)=(180/pi)*phi(nh,taushift(shiftindex));
            rheology.g(nh,shiftindex)=g(nh,taushift(shiftindex));
            rheology.storage(nh,shiftindex)=storage(nh,taushift(shiftindex));
            rheology.loss(nh,shiftindex)=loss(nh,taushift(shiftindex));
            rheology.lambda(nh,shiftindex)=lambda(nh,taushift(shiftindex));
            rheology.delta(nh,shiftindex)=delta(nh,taushift(shiftindex));
            rheology.delf(nh,shiftindex)=real(delfstar(nh,taushift(shiftindex)));
            rheology.delg(nh,shiftindex)=imag(delfstar(nh,taushift(shiftindex)));
            rheology.taushift(nh,shiftindex)=nh*taushift(shiftindex);
        end
    end
elseif get(handles.rheologyrubber, 'value')
    if get(handles.rubberpbd, 'value')
        rubberidx = 23;
    elseif get(handles.rubbersbr, 'value')
        rubberidx = 24;
    else
        error('qcmgradient2:getrheology:invalidRubber','Valid rubber not identified') 
    end
    
    % Load in the DMA data
    filename=['data/S' num2str(rubberidx) '.txt'];
    data=dlmread(filename,',',1,0);
    freq=data(:,1);
    stord=data(:,2);
    lossd=data(:,3);
    if max(freq) > 10^20
        freq(freq>10^20) = NaN;
        stord(freq>10^20) = NaN;
        lossd(freq>10^20) = NaN;
    end
        
    if rubberidx == 23 %This rubber has a point that I want excluded. This is easier than editing the data file for now.
        stord(145) = NaN;
        lossd(145) = NaN;
    end
      
    npoints=300;
    logmin=log(min(freq));
    logmax=log(max(freq));
    
    taushift=logspace(logmin,logmax,npoints); %this is shift in time due to curing
    %These functions turn the data points into a smoother function.
    slmlogstorage = slmengine(log10(freq),log10(stord), 'knots', 8, 'increasing', 'on', 'plot', 'off');
    slmlogloss = slmengine(log10(freq),log10(lossd), 'knots', 8, 'plot', 'off');
    
    storage=@(n,shift) 10.^slmeval(log10(shift*n), slmlogstorage,0);
    loss=@(n,shift) 10.^slmeval(log10(shift*n), slmlogloss,0);
    phi=@(n,shift) atan(loss(n,shift)/storage(n,shift));
    g=@(n,shift) (storage(n,shift)^2+(loss(n,shift))^2).^0.5;
    lambda=@(n,shift) (1/(n*f1))*(g(n,shift)/rho)^0.5*(1/cos(phi(n,shift)/2));
    delta=@(n,shift) lambda(n,shift)*cot(phi(n,shift)/2)/(2*pi);
    fs=@(n) -2*n*f1^2*rho*d/zq;
    delfstar=@(n,shift) fs(n)*tan((2*pi*d/lambda(n,shift))*(1-1i*tan(phi(n,shift)/2)))/...
        ((2*pi*d/lambda(n,shift))*(1-1i*tan(phi(n,shift)/2)));
    
    for nh=rheology.nh
        for shiftindex=1:npoints
            rheology.phi(nh,shiftindex)=(180/pi)*phi(nh,taushift(shiftindex));
            rheology.g(nh,shiftindex)=g(nh,taushift(shiftindex));
            rheology.storage(nh,shiftindex)=storage(nh,taushift(shiftindex));
            rheology.loss(nh,shiftindex)=loss(nh,taushift(shiftindex));
            rheology.lambda(nh,shiftindex)=lambda(nh,taushift(shiftindex));
            rheology.delta(nh,shiftindex)=delta(nh,taushift(shiftindex));
            rheology.delf(nh,shiftindex)=real(delfstar(nh,taushift(shiftindex)));
            rheology.delg(nh,shiftindex)=imag(delfstar(nh,taushift(shiftindex)));
            rheology.taushift(nh,shiftindex)=nh*taushift(shiftindex);
        end
    end
elseif get(handles.rheologycustom, 'value')
    set(handles.ming, 'string', num2str(str2num(get(handles.G1, 'string'))./(rho/1000),'%0.2g'));
    set(handles.maxg, 'string', num2str(str2num(get(handles.G2, 'string'))./(rho/1000),'%0.2g'));
    npoints=2;
    
    try
        grho(3,1) = str2num(get(handles.G1, 'string'))*1000;
        grho(3,2) = str2num(get(handles.G2, 'string'))*1000;
        phi(3,1) = str2num(get(handles.phi1, 'string'));
        phi(3,2) = str2num(get(handles.phi2, 'string'));
    catch
        warndlg('No values entered. Please enter values before continuing')
        return
    end
    
    taushift=[1,2]; %this is shift in time due to curing
    
    g = @(n, shift) grho(3,shift)/rho*(n^(phi(3,shift)/90))/(3^(phi(3,shift)/90));
    phi(1,1:2) = phi(3,1:2);
    phi(5,1:2) = phi(3,1:2);
    storage = @(n,shift) sqrt((g(n,shift)^2)/(1+tand(phi(n,shift))^2));
    loss = @(n, shift) sqrt((g(n,shift)^2)/(1+(1/tand(phi(n,shift))^2)))  ;
    
    lambda=@(n,shift) (1/(n*f1))*(g(n,shift)*rho)^0.5*(1/cosd(phi(n,shift)/2));
    delta=@(n,shift) lambda(n,shift)*cotd(phi(n,shift)/2)/(2*pi);
    fs=@(n) -2*n*f1^2*rho*d/zq;
    delfstar=@(n,shift) fs(n)*tan((2*pi*d/lambda(n,shift))*(1-1i*tand(phi(n,shift)/2)))/...
        ((2*pi*d/lambda(n,shift))*(1-1i*tand(phi(n,shift)/2)));
    
    for nh=rheology.nh
        for shiftindex=1:npoints
            rheology.phi(nh,shiftindex)=phi(nh,taushift(shiftindex));
            rheology.g(nh,shiftindex)=g(nh,taushift(shiftindex));
            rheology.storage(nh,shiftindex)=storage(nh,taushift(shiftindex));
            rheology.loss(nh,shiftindex)=loss(nh,taushift(shiftindex));
            rheology.lambda(nh,shiftindex)=lambda(nh,taushift(shiftindex));
            rheology.delta(nh,shiftindex)=delta(nh,taushift(shiftindex));
            rheology.delf(nh,shiftindex)=real(delfstar(nh,taushift(shiftindex)));
            rheology.delg(nh,shiftindex)=imag(delfstar(nh,taushift(shiftindex)));
            rheology.taushift(nh,shiftindex)=nh*taushift(shiftindex);
        end
    end  
else
    error('qcmgradient2:getrheology:invalidModel','Valid model not identified')
end


function struct2var(s)
%STRUCT2VAR Convert structure array to workspace variables.
%   STRUCT2VAR(S) converts the M-by-N structure S (with P fields)
%   into P variables defined by fieldnames with dimensions M-by-N.  P
%   variables are placed in the calling workspace.
%
%   Example:
%     clear s, s.category = 'tree'; s.height = 37.4; s.name = 'birch';
%     c = struct2cell(s); f = fieldnames(s);
%
%   See also STRUCT2CELL, FIELDNAMES.

% Copyright 2010 The MathWorks, Inc.
if nargin < 1
    error('struct2var:invalidInput','No input structure found')
elseif nargin > 1
    error('struct2var:invalidInput','Too many inputs')
elseif ~isstruct(s)
    error('struct2var:invalidInput','Input must be a structure data type')
end

%[~,c] = size(s);
names = fieldnames(s);

for i=1:length(names)
    assignin('caller',names{i},s.(names{i}))
end

function fontsize_CreateFcn(hObject, eventdata, handles)

function fontsize_Callback(hObject, eventdata, handles)  % change font size of gui
switch get(handles.fontsize,'Value')
    case 1
        fontsize=8;
    case 2
        fontsize=10;
    case 3
        fontsize=12;
    case 4
        fontsize=14;
    case 5
        fontsize=16;
    case 6
        fontsize=18;
    case 7
        fontsize=20;
    case 8
        fontsize=22;
    case 9
        fontsize=24;
    otherwise
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize,'FontUnits','points')
guidata(hObject, handles)

function uselayerlimits_Callback(hObject, eventdata, handles)

function autoupdate_Callback(hObject, eventdata, handles)

function scaley2_SelectionChangeFcn(hObject, eventdata, handles)
if get(handles.y2lin,'value')
    set(handles.axes2,'yscale','linear')
elseif get(handles.y2log,'value')
    set(handles.axes2,'yscale','log')
end

function scalex2_SelectionChangeFcn(hObject, eventdata, handles)
if get(handles.x2lin,'value')
    set(handles.axes2,'xscale','linear');
elseif get(handles.x2log,'value')
    set(handles.axes2,'xscale','log');
end

function scaley1_SelectionChangeFcn(hObject, eventdata, handles)
if get(handles.y1lin,'value')
    set(handles.axes1,'yscale','linear');
elseif get(handles.y1log,'value')
    set(handles.axes1,'yscale','log');
end

function cursorstart_Callback(hObject, eventdata, handles)
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [x,y,but] = ginput(1);
    set(handles.xvalue,'string',num2str(x,'%6.2e'));
    % use commanumber for delf, delg (positions 6,7 in y2dataselect
    if sum(ismember([6 7],get(handles.y2dataselect,'value')))
        set(handles.yvalue,'string',commanumber(y));
    % format to tenth of a degree for phi
    elseif get(handles.y2dataselect,'value')==2
        set(handles.yvalue,'string',num2str(y,'%0.1f'));
    else
        set(handles.yvalue,'string',num2str(y,'%0.2e'));
    end
end

function [str]=commanumber(num)
num=round(num);
str = num2str(num);
fin = length(str);
for i = fin-2:-3:2
    str(i+1:end+1) = str(i:end);
    str(i) = ',';
end
newfin=length(str);
% avoid -, output at beginning of string
if strcmp(str(1),'-') && strcmp(str(2),',')
    str=[str(1) str(3:newfin)];
end

function nh1_Callback(hObject, eventdata, handles)
getharmonics(hObject, handles);

function nh3_Callback(hObject, eventdata, handles)
getharmonics(hObject, handles);

function nh5_Callback(hObject, eventdata, handles)
getharmonics(hObject, handles);

function [nh, legendtext]=getharmonics(hObject, handles)
nh=[];
legendtext={};
if get(handles.nh1,'value')
    nh=[nh 1];
    legendtext=[legendtext; '5MHz'];
end
if get(handles.nh3,'value')
    nh=[nh 3];
    legendtext=[legendtext; '15MHz'];
end
if get(handles.nh5,'value')
    nh=[nh 5];
    legendtext=[legendtext; '25MHz'];
end


guidata(hObject, handles)

function handles = findsolution(hObject, eventdata, handles)
switch get(handles.onelayernh, 'value')
    case 1
        nhset = [1 3 1 5];
    case 2
        nhset = [1 3 3 5];
    case 3
        nhset = [1 5 1 3];
    case 4
        nhset = [1 5 5 3];
    case 5
        nhset = [3 5 3 1];
    case 6
        nhset = [3 5 5 1];
end

%Reads the values of df and dg off the screen
for nh = 1:2:5
    fname = ['delf',num2str(nh)];
    gname = ['delg',num2str(nh)];
    df(nh) = str2double(get(handles.(fname),'string'));
    dg(nh) = str2double(get(handles.(gname),'string'));
end

f1 = handles.f1;
zq = handles.zq;

dissratio = dg(nhset(3))/df(nhset(3));
harmratio = nhset(1)*df(nhset(2))/(nhset(2)*df(nhset(1)));

ming = str2num(get(handles.ming, 'string'));
maxg = str2num(get(handles.maxg, 'string'));

% All properties defined for 15MHz.
grho = 1e3 * (1/ming + 1/maxg)^(-1) * str2num(get(handles.rho, 'string'));
drho = 1e-3*str2num(get(handles.rho, 'string')) * str2num(get(handles.thickness, 'string'));
phi = 45; % I can't at the moment figure out how to get the 5MHz value.
phir = phi*pi/180;
lambdarho=1/(3*5e6)*(grho^0.5)/cosd(phi/2);
d=drho/lambdarho;

x0(1) = phi;
x0(2) = d;

ftosolve = @(x) fstar(x,nhset,harmratio, dissratio, 3);
options = optimset('Display','off');
[x,fval, exitflag,output, jacobian] = fsolve(ftosolve,x0,options); % Call solver

%If the function did not solve satisfactorily (exitflag ~= 1), set data
%to NaNs also.
if exitflag ~= 1 || x(1) > 90
    set(handles.drho1,'string', '--')
    set(handles.grho1,'string', '--')
    set(handles.phi,'String', '--')
    
    disp(['Did not solve using ratio '...
        num2str(nhset(1)) ':' num2str(nhset(2)) ',' num2str(nhset(3))...
        ' with error ' num2str(exitflag)])
    % Also discard data if the phi value is impossible.
else
    phi=x(1); %#ok<*SAGROW>
    dref=x(2);
    phir=phi*pi/180;  % phase angle in radians
    nh = 3;
    refG = 3;
    
    d(nh) = dref;
    drho=(df(nh)/real(delfstar(d(nh),phi)))*zq/(2*nh*f1^2);  %#ok<*AGROW>
    grhoref=((drho/d(nh))*nh*f1*cos(phir/2))^2;
        
    for n = 1:2:5
        lambdarho(n)=lambdarhof(refG, n, grhoref, phi);
        deltarho(n)=lambdarho(n)*(1/(2*pi))*cotd(phi/2);
        grho(n)=grhof(refG, n, grhoref, phi);
        d(n)=drho/lambdarho(n);
        dfds(n)=real(delfstar(d(n),phi));
        dfcalc(n)=sauerbrey(n,drho)*dfds(n);
        dgcalc(n)=sauerbrey(n,drho)*imag(delfstar(d(n),phi));
        if get(handles.modelparameters, 'value')
            set(handles.(['delfcalc' num2str(n)]), 'string', commanumber(real(dfcalc(n))));
            set(handles.(['delgcalc' num2str(n)]), 'string', commanumber(real(dgcalc(n))));
        end
    end
    
    
    % update the data table on the gui with the updated values
    drho=drho*1e3; % convert to g/m^2
    grho=grho*1e-3;  % convert to Pa-g/cm^3
    lambdarho=1e6*lambdarho;  % convert to microns
    set(handles.drho1,'string',num2str(drho,'%8.3f'))
    set(handles.grho1,'string',num2str(grho(nh),'%8.2e'))
    set(handles.phi,'String',num2str(phi,'%8.2f'))
end

function F=sauerbrey(n,drho)
% Calculates the sauerbry shift based on the harmonic and areal density.
F=2*n*5e6^2*drho/8.84e6;

function F = fstar(x,nh,harmratio,dissratio, refn)
% Compares ideal and measured values for fstar for use in the solve
% function. The ideal output of the function is [0 0]. x(1) is the phase
% angle and x(2) is d/lambda.
phi = x(1);  % phase angle in degrees
drefn = x(2);
d(1)=drefn*(nh(1)/refn)^(1-phi/180); % d/lambda at nh(1)
d(2)=drefn*(nh(2)/refn)^(1-phi/180); % d/lambda at nh(2)
d(3)=drefn*(nh(3)/refn)^(1-phi/180); % d/lambda at nh(3)

for n=1:3 
    delf(n)=delfstar(d(n),phi);
end
harmratiocalc=real(delf(2))/real(delf(1));
dissratiocalc=imag(delf(3))/real(delf(3));
F=[dissratio-dissratiocalc;
    harmratio-harmratiocalc]; %difference in ratios.

function F=delfstar(d,phi)  % input phase angle is in degrees
% Calculates delfstar (rhs of delf/delfsn equation) with input of d/lambda 
% and phi
phir=phi*pi/180;
F=-(1/((2*pi*d)*(1-1i*tan(phir/2))))* ...
    tan(2*pi*d*(1-1i*tan(phir/2)));

function onelayernh_Callback(hObject, eventdata, handles)
findsolution(hObject, eventdata, handles);

function onelayernh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function openhandles_Callback(hObject, eventdata, handles)
% Query handles structure. Type return in the command line to return to
% normal operation again.
keyboard;

function rheologymodel_SelectionChangedFcn(hObject, eventdata, handles)
modeltype = eventdata.NewValue.String;
if strcmp(modeltype, eventdata.OldValue.String)
    return
end
switch modeltype
    case 'Simon'
        set(handles.choosesimon, 'visible', 'on')        
        set(handles.chooserubber, 'visible', 'off')
        set(handles.choosecustom, 'visible', 'off') 
    case 'Rubber'
        set(handles.choosesimon, 'visible', 'off')
        set(handles.choosecustom, 'visible', 'off') 
        set(handles.chooserubber, 'visible', 'on')
    case 'Custom'
        set(handles.choosecustom, 'visible', 'on')
        set(handles.chooserubber, 'visible', 'off') 
        set(handles.choosesimon, 'visible', 'off') 
end
handles.rheology=getrheology(hObject, handles);
plotgradient(hObject,handles);
guidata(hObject, handles)

function chooserubber_SelectionChangedFcn(hObject, eventdata, handles)
handles.rheology=getrheology(hObject, handles);
plotgradient(hObject,handles);
guidata(hObject, handles)

function repeats_Callback(hObject, eventdata, handles)
plotgradient(hObject,handles);
findsolution(hObject, eventdata, handles);

function repeats_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G1_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles);
    findsolution(hObject, eventdata, handles);
    guidata(hObject, handles)
end

function G1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function phi1_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles) ;
    findsolution(hObject, eventdata, handles);
    guidata(hObject, handles)
end

function phi1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function phi2_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles) ;
    findsolution(hObject, eventdata, handles);
    guidata(hObject, handles)
end

function phi2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G2_Callback(hObject, eventdata, handles)
if get(handles.autoupdate,'value')
    handles.rheology=getrheology(hObject, handles); 
    plotgradient(hObject,handles); 
    findsolution(hObject, eventdata, handles);
    guidata(hObject, handles)
end

function G2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lrho = lambdarhof(refn, n, grho, phi)
if refn == n
    lrho = 1/(n*5e6)*(grho^0.5)/cosd(phi/2);
else
    lrho = 1/(n*5e6)*((grho*(n^(phi/90))/(refn^(phi/90)))^0.5)/cosd(phi/2);
end

function grho = grhof(refn, n, grhoref, phi)
grho = grhoref*(n^(phi/90))/(refn^(phi/90));

function dl = dlf(refn, refnh, dlin, phi);
dl = dlin*(refnh/refn)^(1-(phi/180));

function layersweep_Callback(hObject, eventdata, handles)
ntimes = str2num(get(handles.repeats, 'string'));
phi = []; grho = [];
figure
h(1) = subplot(1,2,1);
hold on
xlim([0 2*ntimes]);
xlabel('Number of layers')
ylabel('G_3')
h(2) = subplot(1,2,2);
hold on
xlim([0 2*ntimes]);
xlabel('Number of layers')
ylabel('\phi')
for i = 1:ntimes
   set(handles.repeats, 'string', num2str(i)); %Set repeats to the correct number
   plotgradient(hObject, handles);
   handles = calcqcm(hObject, handles);
   handles = findsolution(hObject, eventdata, handles);
   drawnow;
   pause(.5);
   grho = [grho str2num(get(handles.grho1, 'string'))];
   phi = [phi str2num(get(handles.phi, 'string'))];
   layers = 2:2:2*i;
   plot(h(1), layers, grho, 'bx');
   plot(h(2), layers, phi, 'bx');
end
assignin('base', 'grho', grho)
assignin('base', 'phi', phi)

% --- Executes on button press in modelparameters.
function modelparameters_Callback(hObject, eventdata, handles)
%Switches which panel is showing for entering the layer properties and
%thicknesses.
if get(handles.modelparameters, 'value')
    set(handles.modelpanel, 'visible', 'on')
    set(handles.fitpanel, 'visible', 'off')
else
    set(handles.modelpanel, 'visible', 'off')
    set(handles.fitpanel, 'visible', 'on')
end

for n = 1:2:5 
%Reset all of the top values to blank to indicate nothing has been solved yet
    set(handles.(['delfcalc' num2str(n)]), 'string', '-');
    set(handles.(['delgcalc' num2str(n)]), 'string', '-');
    set(handles.(['delf' num2str(n)]), 'string', '-');
    set(handles.(['delg' num2str(n)]), 'string', '-');
end

% --- Executes when entered data in editable cell(s) in layerpropertytable.
function layerpropertytable_CellEditCallback(hObject, eventdata, handles)
%Whenever a layer is added or changed, the total thickness changes. This
%updates the total thickness of the film to be equal to the current total
%thickness of all of the layers.
layerdata = handles.layerpropertytable.Data;
nlayers = find(cellfun('isempty', layerdata),1)-1;
if isempty(nlayers) %If all of the layers are full by some chance...
        nlayers = max(size(layerdata)); %This should be however many rows there are in the table
end
thickness = 0;
for i = 1:nlayers
   thickness = thickness + str2num(layerdata{i,1});
end
set(handles.thickness, 'string', num2str(thickness, 3));
% If there is enough data in the table, it should update the figure. If the
% thickness that was just changed has a corresponding modulus value, update
% the figure.
where = eventdata.Indices;
if where(2) == 1 && ~isempty(layerdata{where(1), 2}) %The thickness was updated and there is a corresponding layer data
    update_Callback(hObject, eventdata, handles);
elseif where(2) == 1
    return
else %Layer properties were changed, check that they are valid
   if str2num(layerdata{where(1), where(2)}) < 0 || str2num(layerdata{where(1), where(2)}) > 1
       warndlg(['You have entered a relative modulus value outside the '...
           'limits of 0 to 1. Please enter a value between 0 and 1.']);
       handles.layerpropertytable.Data{where(1), where(2)} = '-';
   else
       update_Callback(hObject, eventdata, handles);
   end
end

function [lowdelfstar, highdelfstar] = calcerrorrange(hObject, handles)
f1=handles.f1;
zq=handles.zq;
handles = checkrheology(hObject, handles);
nh=handles.rheology.nh;
layerprops=getlayerprops(hObject, handles);
if isempty(layerprops)
    warndlg('No layer data entered. Please enter data')
    return
end
highlayerprops = layerprops;
lowlayerprops = layerprops;
pererror = str2num(get(handles.herrorpercent, 'string'));
if isempty(pererror)
    warndlg('You need to enter an error value')
    lowdelfstar = [];
    highdelfstar = [];
    return
elseif pererror > 1 && pererror < 100
    pererror = pererror/100;
elseif pererror > 100
    warndlg('Your error is too big')
    lowdelfstar = [];
    highdelfstar = [];
    return
end

for i = 1:layerprops.nlayers
    highlayerprops.subthickness(i) = layerprops.subthickness(i)*(1+pererror);
    lowlayerprops.subthickness(i) = layerprops.subthickness(i)*(1-pererror);
end

highdelfstar=calcdelfstar(f1,zq,highlayerprops,nh);
lowdelfstar=calcdelfstar(f1,zq,lowlayerprops,nh);
foutputstring = '';
goutputstring = [];
for i = 1:2:5
    foutputstring = sprintf('%s%2.0f - %2.0f\n', foutputstring, real(lowdelfstar(i)), real(highdelfstar(i)));
    goutputstring = sprintf('%s%1.0f - %1.0f\n', goutputstring, imag(lowdelfstar(i)), imag(highdelfstar(i)));
end
set(handles.herrorresults, 'string', [foutputstring goutputstring])

% --- Executes on button press in herror.
function herror_Callback(hObject, eventdata, handles)
[lowdelfstar, highdelfstar] = calcerrorrange(hObject, handles);

function herrorpercent_Callback(hObject, eventdata, handles)

function herrorpercent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [handles, rheology] = checkrheology(hObject, handles)
try
    rheology = handles.rheology;
catch
    rheology = getrheology(hObject, handles);
    handles.rheology = rheology;
end
