function varargout = test(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_OpeningFcn, ...
                   'gui_OutputFcn',  @test_OutputFcn, ...
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


function test_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

function varargout = test_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%function pushbutton2_Callback(hObject, eventdata, handles)
function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
N=get(handles.edit3,'string');
N=str2num(N);


function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton12_Callback(hObject, eventdata, handles)
close;


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
luachon = get(handles.luachon,'value');
if luachon==1
    N=get(handles.edit3,'string');
    N=str2num(N);
    dv=randi([0 1],1,N*4);
    dvd1=num2str(dv);
    set(handles.edit1,'string',dvd1);
    sobitvao=num2str(N*4);
    set(handles.sobitvao,'string',sobitvao)
else
    dv=eval(get(handles.edit1,'string'));
    NN=length(dv);
    sobitvao=num2str(NN);
    set(handles.sobitvao,'string',sobitvao)
    if rem(NN,4)==1
        N=(NN+3)/4;
        dv=[dv 0 0 0 ];
    elseif rem(NN,4)==2
        N=(NN+2)/4;
        dv=[dv 0 0];
    elseif rem(NN,4)==3
        N=(NN+1)/4;
        dv=[dv 0];
    else
        N=NN/4;
        dv=dv;
    end
end


Frame=2;
 
H1=0.971960;
H2=sqrt(2)/2;
H3=0.235147;
factech=1+2*(H1+H2+H3);
hef(1:4*N)=0;
for i=1:4*N-1
    hef(1+i)=1-2*H1*cos(pi*i/(2*N))+2*H2*cos(pi*i*2/(2*N))-2*H3*cos(pi*i*3/(2*N));
end
hef=hef/factech;
h=hef;
y=zeros(1,4*Frame*N);
sdemod=zeros(N,Frame);
 
 
 
c=reshape(dv,N*2,2);
 for e=1:(N*2)
    if c(e,:)==[0 0]
        g(e)=-1;
    elseif c(e,:)==[1 1]
        g(e)=1; 
    elseif c(e,:)==[0 1]
        g(e)=j;
    else
        g(e)=(-j);
    end
 end
 gnv=g;
 
 
 s=reshape(g,N,2);
 
 
 
 for ntrame=1:Frame
 aa=ones(N,N);
 [m n] =size(aa);
 for i=1:m
    for ii=1:n
        W(i,ii)=exp(-2*pi*j*(ii-1)*(i-1)/N)*aa(i,ii);
    end
end
x(:,ntrame)=(1/N)*W*s(:,ntrame);
x4=[x(:,ntrame).' x(:,ntrame).' x(:,ntrame).' x(:,ntrame).'];
x5=[x(:,1).' x(:,1).' x(:,1).' x(:,1).'];
signal=x4.*h;
signal2=x5.*h;
 
 
 
%Parallel to Serial Conversion 
y(1+(ntrame-1)*4*N:4*ntrame*N)=y(1+(ntrame-1)*4*N:4*ntrame*N)+signal;
y2(1+(ntrame-1)*4*N:4*ntrame*N)=y(1+(ntrame-1)*4*N:4*ntrame*N)+signal2;
 end
 
 
%%ChANNEL
channel =y
 
%%%Serial to Paralll Conversion
for ntrame=1:Frame
    r=channel(1+(ntrame-1)*4*N:4*ntrame*N).*h;
    u=zeros(2,N);
    for k=1:4
        u(ntrame,:)=u(ntrame,:)+r(1,1+(k-1)*N:k*N);
    end
    
    %%%ifft 
    for e= 1:m
        for ee=1:n
            W2(e,ee)=exp(2*pi*j*(ee-1)*(e-1)/N)*aa(e,ee);
        end
    end
    sest(ntrame,:)=W2*u(ntrame,:).';
    sest(ntrame,:)=sest(ntrame,:)/0.6863;
    sest(ntrame,:)=round(sest(ntrame,:));
end
gg=sest.';
sdemod=gg;
 
%%Parallel to serial conversion
gnr=reshape(gg,1,N*2);
 
ra=zeros(N*2,2);
for rr=1:(N*2)
    if gnr(rr)==-1
        ra(rr,:)=[0 0];
    elseif gnr(rr)==1
        ra(rr,:)=[1 1];
    elseif gnr(rr)==j
        ra(rr,:)=[0 1];
    else
        ra(rr,:)=[1 0];
    end
end
ra2=reshape(ra,1,N*4);
 
if luachon==2
    if rem(NN,4)==1
    for i=1:(4*N-3)
        ngora(i)=ra2(i);
    end
    elseif rem(NN,4)==2
    for i=1:(4*N-2)
        ngora(i)=ra2(i);
    end
    elseif rem(NN,4)==3
    for i=1:(4*N-1)
        ngora(i)=ra2(i);
    end
    else
        N=NN/4;
        ngora=ra2;
    end
else
    ngora=ra2;
end
 
%%-----------PHAN CODE DE VE DO THI-------------%%
%%------------Symbol----------%%
for ntrame=1:Frame
    r2=y2(1+(ntrame-1)*4*N:4*ntrame*N).*h;
    
    
    u2=zeros(1,N);
    for k=1:4
        u2=u2+r2(1,1+(k-1)*N:k*N);
    end
    
     u2=u2.';
end
 
%%--------------KET THUC PHAN CODE DE VE DOTHI--------------------------%%
 
 
%%-----------------------------------------------------------%%
 
%%--------------------------------------------%%
%%--------------------------------------------
 
luachon2=get(handles.luachon2,'value');
if luachon2==1
    axes(handles.axes17);
    plot(s,'or');
    grid on;
    title('CHOM SAO PHIA PHAT');
    xlabel('amplitude(I)');
    ylabel('amplitude(Q)');
    xlim([-1.2 1.2]);
    ylim([-1.2 1.2]);
else
    axes(handles.axes17);
    hold off;
    stem(real(gnv),'or');
    hold on
    stem(imag(gnv),'xb');
    grid on;
    title('NGO VAO OQAM');
    xlabel('Frequency');
    ylabel('Amplitude');
    xlim([1 length(gnv)]);
    ylim([-1.2 1.2]);
    hold off 
end
 
%%----------------VE BIEU DO BO LOC--------%%
 
 axes(handles.axes18);
 plot(h);
 grid on;
 title('BO LOC PDF');
 xlabel('Time');
 ylabel('Amplitude');
 xlim([1 length(h)]);
 ylim([-1 1.2]);
 
 %%-----------VE BIEU DO IO OQAM PHIA PHAT------------%%
 
 axes(handles.axes1);
 hold off;
 stem(real(s(:,1)),'b');
 hold on;
 stem(real(s(:,2)),'xr');
 xlim([1 length(s)]);
 ylim([-1.2 1.2]);
 xlabel('Data points I (Frequncy)');
 ylabel('Amplitude');
 grid on;
 hold off;
 
 
 axes(handles.axes2);
 hold off;
 stem(imag(s(:,1)));
 hold on;
 stem(imag(s(:,2)),'xr');
 xlim([1 length(s)]);
 ylim([-1.2 1.2]);
 xlabel('Data points Q (Frequncy)');
 ylabel('Amplitude');
 grid on;
 hold off;
 
 %%------------VE BIEU DO IO IFFT--------------%%
 
 axes(handles.axes3);
 hold off;
 plot(real(x(:,1)),'r')
 hold on;
 plot(real(x(:,2)));
 grid on;
 xlabel('Data points I (Time)');
 ylabel('Amplitude');
 xlim([1 length(x)]);
 ylim([-1 1]);
 grid on;
 hold off;
 
 axes(handles.axes4);
 hold off;
 plot(imag(x(:,1)),'r')
 hold on;
 plot(imag(x(:,2)));
 grid on;
 xlabel('Data points 0 (Time)');
 ylabel('Amplitude');
 xlim([1 length(x)]);
 ylim([-1 1]);
 hold off;
 %%-----VE BIEU DO IO SAU BO LOC-------%%
 axes(handles.axes5);
 hold off;
 plot(real(signal),'b');
 grid on;
 hold on;
 plot(real(signal2),'r');
 legend('symbol 1','Symbol 2')
 xlabel('Time (I)');
 ylabel('Amplitude');
 ylim([-1 1]);
 xlim([1 length(signal)]);
 hold off;
 
 
 axes(handles.axes6);
 hold off;
 plot(imag(signal))
 grid on;
 hold on;
 plot(imag(signal2),'r');
 xlabel('Time (Q)');
 ylabel('Amplitude');
 ylim([-1 1]);
 xlim([1 length(signal)]);
 hold off;
 
 %%---------VE BIEU DO IQ CHANNEL-----------%%
 axes(handles.axes7);
 plot(real(channel))
 grid on;
 xlabel('Time (I)');
 ylabel('Amplitude');
 ylim([-1 1]);
 xlim([1 length(channel)]);
 hold off;
 axes(handles.axes8);
 plot(imag(channel))
 grid on;
 xlabel('Time (Q)');
 ylabel('Amplitude');
 ylim([-1 1]);
 xlim([1 length(channel)]);
 hold off;
 
 %%-------VE BIEU DO IQ BO LOC PHIA THU----------------------%%
 axes(handles.axes11);
 hold off;
 plot(real(r))
 hold on;
 plot(real(r2),'r')
 grid on
 xlabel('Time (I)');
 ylabel('Amplitude');
 ylim([-1 1]);
 xlim([1 length(r2)]);
 hold off;
 
 axes(handles.axes12);
 hold off;
 plot(imag(r))
 hold on;
 plot(imag(r2),'r')
 grid on
 xlabel('Time (Q)');
 ylabel('Amplitude');
 ylim([-1 1]);
 xlim([1 length(r2)]);
 hold off;
 
 
 %%-------VE BIEU DO IQ BO FFT----------------------%%
 axes(handles.axes13);
 hold off;
 u=u(2,:).';
 plot(real(u),'b');
 hold on;
 plot(real(u2),'r');
 grid on;
 xlabel('Data Points I (TIME)');
 ylabel('Amplitude');
 xlim([1 length(u)]);
 ylim([-1 1]);
 hold off;
 
 
 axes(handles.axes14);
 hold off;
 plot(imag(u),'b');
 hold on;
 plot(imag(u2),'r');
 grid on;
 xlabel('Data Points 0 (TIME)');
 ylabel('Amplitude');
 xlim([1 length(u)]);
 ylim([-1 1]);
 hold off;
 
 %%---------VE BIEU DO IQ OQAM THU--------%%
 axes(handles.axes15);
 hold off;
 stem(real(sdemod(:,1)));
 hold on;
 stem(real(sdemod(:,2)),'xr');
 grid on;
 xlabel('Data Points I (frequency)');
 ylabel('Amplitude');
 xlim([1 length(sdemod)]);
 ylim([-1.2 1.2]);
 hold off;
 
 
 axes(handles.axes16);
 hold off;
 stem(imag(sdemod(:,1)));
 hold on;
 stem(imag(sdemod(:,2)),'xr');
 grid on;
 xlabel('Data Points 0 (frequency)');
 ylabel('Amplitude');
 xlim([1 length(sdemod)]);
 ylim([-1.2 1.2]);
 hold off;
 
 %%---------VE BIEU DO CHOM PHIA THU---------%%
 luachon3=get(handles.luachon3,'value');
 if luachon3==1
     axes(handles.axes19);
     plot(sdemod,'or');
     grid on;
     title('CHOM SAO PHIA THU')
     xlabel('Amplitude(I)');
     ylabel('Amplitude(Q)');
     xlim([-1.2  1.2]);
     ylim([-1.2 1.2]);
 else   
      axes(handles.axes19);
      hold off;
      stem(real(gnr),'or');
      hold on;
      stem(imag(gnr),'xb');
      grid on;
      title('NGO RA OQAM')
      xlabel('frequency');
      ylabel('Amplitude');
      xlim([0 length(gnr)+1]);
      ylim([-1.2 1.2]);
      hold off;
 end
 %%HIEN THU DU LIEU RA PHIA THU-----------------%%
 dulieuraall=ngora;
 sobitra=length(ngora);
 sobitra=num2str(sobitra);
 set(handles.sobitra,'string',sobitra);
 %%------hien thi ra Guide-----%%
 dlr=num2str(dulieuraall);
 set(handles.edit2,'string',dlr)
 
 %------hien thi SO DIEM IFFT/FFT-----%%
 N=num2str(N);
 set(handles.edit3,'string',N);


% --- Executes on selection change in luachon2.
function luachon2_Callback(hObject, eventdata, handles)
% hObject    handle to luachon2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns luachon2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from luachon2


% --- Executes during object creation, after setting all properties.
function luachon2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to luachon2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in luachon3.
function luachon3_Callback(hObject, eventdata, handles)
% hObject    handle to luachon3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns luachon3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from luachon3


% --- Executes during object creation, after setting all properties.
function luachon3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to luachon3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in luachon.
function luachon_Callback(hObject, eventdata, handles)
% hObject    handle to luachon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns luachon contents as cell array
%        contents{get(hObject,'Value')} returns selected item from luachon


% --- Executes during object creation, after setting all properties.
function luachon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to luachon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function sobitvao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sobitvao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in sobitra.
function sobitra_Callback(hObject, eventdata, handles)
% hObject    handle to sobitra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function sobitra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sobitra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
