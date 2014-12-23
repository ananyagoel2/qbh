function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 12-Dec-2014 02:26:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

recObj = audiorecorder;

% --- Executes on button press in record.
function record_Callback(hObject, eventdata, handles)
% hObject    handle to record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
recordblocking(recObj, 5);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
[FileName,PathName] = uigetfile('*.wav','Select the hummed audio file');
path=strcat(PathName,'\',FileName);
[data1,sr1]=wavread(path);
plot(data1)
axes(handles.axes2);
x=data1;
q=sr1;
p.maxprd=150;
p.minprd=2;
p.wsize=p.maxprd;
 p.hop=p.wsize;
 p.thresh=0.1;
 p.smooth=p.minprd/2;
if min(size(x)) ~= 1; error('data should be 1D'); end
x=x(:);
nsamples=numel(x);
 
nframes=floor((nsamples-p.maxprd-p.wsize)/p.hop);
pwr=zeros(1,nframes);
prd=zeros(1,nframes);
ap=zeros(1,nframes);
 
% shifted data
x=convmtx(x,p.maxprd+1);
x=x(p.maxprd:end-p.maxprd,:);
 nframes
for k=1:nframes
   
    start=(k-1)*p.hop; % offset of frame
    xx=x(start+1:start+p.wsize,:);
   
    d=mean( (xx - repmat(xx(:,1),1,p.maxprd+1)).^2 )/2;  
     % squared difference function
    dd= d(2:end) ./ (cumsum(d(2:end)) ./ (1:(p.maxprd)));   % cumulative mean - normalized
    
    % parabolic interpolation of all triplets to refine local minima
    min_pos=1:numel(dd);    % nominal position of each sample
    x1=dd(1:end-2);
    x2=dd(2:end-1);
    x3=dd(3:end);
    a=(x1+x3-2*x2)/2;
    b=(x3-x1)/2;
    shift=-b./(2*a);        % offset of interpolated minimum re current sample
    val=x2-b.^2./(4*a);     % value of interpolated minimum
    
    % replace all local minima by their interpolated value, 
    idx= 1 + find(x2<x1 & x2<x3);
    dd(idx)=val(idx-1);
    min_pos(idx)=min_pos(idx-1)+shift(idx-1);
    
    % find index of first min below threshold
    a=dd<p.thresh;
    if isempty(find(a))
        [~,prd0]=min(dd); % none below threshold, take global min instead
    else
        b=min(find(a)); % left edge
        c=min(b*2,numel(a));
        [~,prd0]=min(dd(b:(c-1))); 
        prd0=b+prd0-1;
    end
    
    prd=min_pos(prd0)+1;
        
    if prd>2 & prd<numel(dd) & d(prd0)<d(prd0-1) & d(prd0)<d(prd0+1)
        
        % refine by parabolic interpolation of raw difference function 
        k
        
        x1=d(prd-1);
        x2=d(prd);
        x3=d(prd+1);
        a=(x1+x3-2*x2)/2;
        b=(x3-x1)/2;
        shift=-b./(2*a);        % offset of interpolated minimum re current sample
        val=x2-b.^2./(4*a);     % value of interpolated minimum
        prd=prd+shift-1;
    end
%     
    % aperiodicity
    frac=prd-floor(prd);
    if frac==0
        yy=xx(:,prd);
    else
        yy1=(1-frac)*xx(:,floor(real(prd+1)))+frac*xx(:,floor(real(prd+1))+1); % linear interpolation
       yy=yy1;
    end
    pwr=(mean(xx(:,1).^2) + mean(yy.^2))/2; % average power over fixed and shifted windows
    res=mean(((xx(:,1) - yy)).^2) / 2;
    ap=res/pwr;
            
    r.prd(k)=prd;
    r.ap(k)=ap;
    r.prd(k);

   
 
end
C = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890';
for i=1:nframes
     y(i)=1000/r.prd(i);
     %g(i)=music.freq2note(y(i))
     t(i)=music.freq2tone(y(i));
        if (t(i)>=-49&&t(i)<=12)
            
            if(t(i)<0)
             
             n(i)=C(50+t(i));
            else 
            n(i)=C(50+t(i));
            end
        
            
        end
     end

    distance=n;
    
    [string1, string2,string3, string4, string5] = textread('mydata.dat', ...
'%q %q %q %q %q', 1);
    string11=char(string1) ;
    string22=char(string2) ;
    string33=char(string3) ;
    string44=char(string4) ;
    string55=char(string5) ;
   %d1=LCS('abc',string11)
    chord1=LCS(distance,string11)
    chord2=LCS(distance,string22)
    chord3=LCS(distance,string33)
    chord4=LCS(distance,string44)
    chord5=LCS(distance,string55)
    compare_array=[chord1,chord2,chord3,chord4,chord5];
    [g,i]=max(compare_array)
    set(handles.text2,'string',{'Chord' num2str(i) ' ' 'Chord1 = A' 'Chord2 = B' 'Chord3 = C' 'Chord4 = D' 'Chord5 = E'})


    %subplot 211; plot(r.prd); title('period'); xlabel('frame'); ylabel('samples');
    %subplot 212; plot(r.ap); title('periodicity'); xlabel('frame');
%     subplot 313; plot(sort(r.prd)); title('period'); xlabel('frame'); ylabel('samples');

   
plot(r.prd); title('period'); xlabel('frame'); ylabel('samples');   
axes(handles.axes3);
plot(r.ap); title('periodicity'); xlabel('frame');
axes(handles.axes4);
r=[];

 



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot(data1);
axes(handles.axes3);
