% nevopen_memmap
% Created by AZ 2009-04-15

filename = ['..' filesep '..'      filesep 'Data' filesep 'Cerebus' ...
                 filesep 'CATZ077' filesep 'u004_022.nev'];
elec = 1;

% Specify STANDARD HEADER file format (336 BYTES)
nevformat = {'uint8'  [1   8]             'identifier'; ...
             'uint8'  [1   2]               'filespec'; ...
             'uint16' [1   1]              'add_flags'; ...
             'uint32' [1   1]          'nbytes_header'; ...
             'uint32' [1   1] 'nbytes_per_data_packet'; ...
             'uint32' [1   1]        'timeres_tstamps'; ...
             'uint32' [1   1]        'timeres_samples'; ...
             'uint16' [1   8]             'timeorigin'; ...
             'uint8'  [1  32]               'app_used'; ...
             'uint8'  [1 256]                'comment'; ...
             'uint32' [1   1]        'num_ext_headers'     };
% Open file
nev = memmapfile(filename,'Format',nevformat,'Repeat',1);

% SPEW OUT STANDARD HEADER DATA
nev.Data
% char(nev.Data.identifier)
% nev.Data.filespec
waveform_length_flag   = str2double(dec2bin(nev.Data.add_flags));
nbytes_header          = double(nev.Data.nbytes_header);
nbytes_per_data_packet = double(nev.Data.nbytes_per_data_packet);
timeres_tstamps        = double(nev.Data.timeres_tstamps);
timeres_samples        = double(nev.Data.timeres_samples);
% nev.Data.timeorigin
% char(nev.Data.app_used)
% char(nev.Data.comment)
num_ext_headers        = double(nev.Data.num_ext_headers);

if (nbytes_header-336 ~= 32*num_ext_headers)
   warning('Unexpected extended header size'); %#ok<WNTAG>
end

% Specify EXTENDED HEADER file format
nevformat = {'uint8'  [1   1]             'identifier'; ...
             'uint16' [1   1]             'info_field'; ...
             'uint8'  [1   1]                  'empty';    };
% Open file
nev = memmapfile(filename,'Format',nevformat,'Repeat',6,'Offset',336+32*7);

% SPEW OUT EXTENDED HEADER DATA
% dig_scaling_factor_nV = double( nev.Data(1).info_field);
% energy_thrshold       = double( nev.Data(2).info_field);
% amp_hi_thresh_in_uV   = double( nev.Data(3).info_field);
% amp_lo_thresh_in_uV   = double(-nev.Data(4).info_field);
% nsorted_units         = double( nev.Data(5).info_field);
nbytes_per_waveform   = double( nev.Data(6).info_field);

if waveform_length_flag && (nbytes_per_waveform ~= 2)
   warning('Unexpected waveform length'); %#ok<WNTAG>
elseif ~waveform_length_flag && (nbytes_per_waveform ~= 2)
   warning(['Nonstandard waveform length: ',...
      num2str(nbytes_per_waveform),' bytes']); %#ok<WNTAG>
end

% Specify DATA file format
nevformat = {'uint32' [1   1]              'timestamp'; ...
             'uint16' [1   1]               'packetid'; ...
             'uint8'  [1   1]         'unit_class_num'; ...
             'uint8'  [1   1]               'reserved'; ...
             'int16'  [1 (nbytes_per_data_packet-8)/2]  ...
                                            'waveform';    };

% Calculate .nev file size, initialize variables
nevstats = dir(filename);
num_data_samples = (nevstats.bytes - nbytes_header)/nbytes_per_data_packet;
packetid         = zeros(num_data_samples,1);

% Open file, populate variable arrays
N = 500000;
N_end = N;
N_last = mod((nevstats.bytes - nbytes_header),N*nbytes_per_data_packet);
i = -1;
j = 0;

% Initialize data structs/arrays
data = repmat(struct('timestamp', [], 'waveform', []),96,1);

while i < floor(num_data_samples/N)
   i = i + 1;
   if      (nevstats.bytes - (nbytes_header+i*N*nbytes_per_data_packet)) < ...
                                              N*nbytes_per_data_packet
   N_end = (nevstats.bytes - (nbytes_header+i*N*nbytes_per_data_packet)) / ...
                                                nbytes_per_data_packet;
   end
   nev = memmapfile(filename,'Format',nevformat,...
      'Offset',nbytes_header+i*N*nbytes_per_data_packet,'Repeat',N_end);
   nevdata = rmfield(nev.Data,{'reserved';'unit_class_num'});
   
   packetid(i*N+1:i*N+N_end) = double([nevdata.packetid]');
   
   % SAVE DATA
%    nonzero_ids = find(packetid(i*N+1:i*N+N_end));
   for j = elec
      this_elec = find(packetid(i*N+1:i*N+N_end) == j);
      data(j).timestamp = [data(j).timestamp; ...
                   double([nevdata(this_elec).timestamp]')                        ];
      data(j).waveform  = [data(j).waveform ; ...
           reshape(double([nevdata(this_elec).waveform ]'),48,size(this_elec,1))' ];
   end
%    data(packetid(j)).timestamp = [data(packetid(j)).timestamp; ...
%                        double([nevdata(nonzero_ids).timestamp]')                         ];
%    data(packetid(j)).waveform  = [data(packetid(j)).waveform ; ...
%                reshape(double([nevdata(nonzero_ids).waveform]'),48,size(nonzero_ids,1))' ];
   
   disp([num2str((i*N+N_end)*100/num_data_samples),'% done, at sample # ', ...
      num2str(i*N+N_end)]);
end


% nevdata  = nev.Data;
% packetid = zeros(size(nevdata,1),1);
% for i = 1:size(nevdata,1)
%    packetid(i) = nevdata(i).packetid;
% end
% packetid = double(packetid);
% % figure; hold on;
% % for i = 1:N
% %    plot(double(nev.Data(i).waveform));
% % end