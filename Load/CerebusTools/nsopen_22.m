function [r,SamplingRateInKHZ,nChans,tempNS] = nsopen_22(filename)%,d)

% 13-Aug-2014 BH created to deal with opening filespc 2.2 ns5 files
% builds on nsopen.m

% Created DLR who knows when
% 2010-03-16 AZ Tweaked some code, added memmapfile at end


% global pepNEV;
% 
% if ~isempty(pepNEV)
%     nevclose;
% end
% filename = sprintf('%s.ns%u',filename,d);
fid = fopen(filename, 'r', 'ieee-le');
if fid == -1, 
    r = -1;
    return;
end
% AZ 2009-03-04
% pepNEV.index.fid = fid;

%%%%%%%%%%%  READ HEADER (see Cerebus manual for info/naming)  %%%%%%%%%%%

if fseek(fid,0,'bof')
   r = -1;
   return;
end

fileid = fread(fid,8,'*char')';
% if ~strcmp(fileid,'NEURALSG')
if ~strcmp(fileid,'NEURALCD') %filespec 2.2 format
   r = -1;
   return;
end

ExtHeaderLength = 66; %hard-coded from cerebus documentation
BasicHeader = fread(fid, 306, '*uint8');
HeaderBytes = double(typecast(BasicHeader(3:6), 'uint32'));
TimeRes = double(typecast(BasicHeader(283:286), 'uint32'));
SamplingFreq = TimeRes / double(typecast(BasicHeader(279:282), 'uint32'));
SamplingRateInKHZ = SamplingFreq;
ChannelCount = double(typecast(BasicHeader(303:306), 'uint32'));
readSize = double(ChannelCount * ExtHeaderLength);
ExtendedHeader = fread(fid, readSize, '*uint8');

EOexH = double(ftell(fid));
fseek(fid, 0, 'eof');
EOF = double(ftell(fid));

% label = fread(fid,16,'*char')';
% a = findstr(label,' kS/s');
% label = str2double(label(1:a-1));
% if label == 0 or isempty(a), "LFP Low" sampling rate
% if not LFP, label and SamplingRateInKHZ should be equal

disp([filename,sprintf(' sampled at %ukHz',SamplingRateInKHZ)]);

% nChans = fread(fid,1,'uint32=>double');
nChans = double(typecast(BasicHeader(303:306), 'uint32'));
% pepNEV.nchan = nChans;

% pepNEV.ids = fread(fid,nChans,'uint32=>double');
% pepNEV.ids = EOexH;
% This is the end of the header / beginning of the data section

% nbytesInHeader = 32 + 4*nChans;
nbytesInHeader = double(typecast(BasicHeader(3:6), 'uint32'));
% tempNS = openNSx(filename,'noread','report');
tempNS = openNSx(filename);
% DataPoints = double((EOF-EOexH)/(nChans*2)); %filespec 2.2 has 9 preceding bytes of extra header
DataPoints = tempNS.MetaTags.DataPoints;
% Specify DATA file format
% fseek(fid,0,'eof');
nbytes = EOF;
fclose(fid);
% fseek(fid,nbytesInHeader,'bof');

% nsformat = { 'int16' [nChans (nbytes-nbytesInHeader)/(nChans*2)] 'data';};
nsformat = { 'int16' [nChans DataPoints] 'data';};
% Open file
% pepNEV.ns = memmapfile(filename,'Format',nsformat,'Offset',nbytesInHeader);

% clear tempNS

r = 1;