function [ DLP ] = initializeDLP(varargin)
%this script sets all the projectors to HDMI Video Input
%then sets their video mode settings to 60 hz 8 bit and monochrome green

list = 1:5;
current = 1;
bitDepth = 7;
color = 3;

for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

if current > 1
    current = 1;
end

current = round(current*274);

numDLP = length(list);
DLP = LightCrafter();

%third byte of the IP address
IP = {'1','2','3','4','5'};

for ii = 1:numDLP
    disp(list(ii));
    try
        tcpObject = tcpip(['192.168.',IP{list(ii)},'.110'],21845);
        fopen(tcpObject);
    catch err
        error(['rig ',num2str(list(ii)),' didn''t connect']);
    end

    DLP.setDisplayModeHDMIVideoInput(tcpObject);
    pause(0.01);
    DLP.setVideoModeSetting(tcpObject,bitDepth,color);
    pause(0.01);
    DLP.setLEDcurrent(tcpObject,current);
    pause(0.01);
    fclose(tcpObject);
end

end