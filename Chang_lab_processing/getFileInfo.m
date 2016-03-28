folders = {'EC108_2fd7f1ac';'EC108_9ae22432';'EC108_e6336384'; ...
    'EC108_03c06afe';'EC108_ff348c69';'EC108_8337e405';'EC108_cd686c40'};
timelist=struct([]);

for d=1:length(folders)
  
    h5files=dir(sprintf('%s/*.h5', char(folders(d))));
    
    for f=1:length(h5files)
        filename=sprintf('%s/%s',char(folders(d)), char(h5files(f).name));
        timestamps = h5read(filename, '/timestamp vector'); % These are POSIX (UNIX) timestamps
        fs = h5readatt(filename, '/ECoG Array', 'Sampling Rate');
        start_datetime = datetime(timestamps(1), 'TimeZone', 'America/Los_Angeles', 'ConvertFrom', 'posixtime');
        end_datetime = datetime(timestamps(end), 'TimeZone', 'America/Los_Angeles', 'ConvertFrom', 'posixtime');
        
        this.dir = char(folders(d));
        this.fname = char(h5files(f).name);
        this.fs = fs;
        this.start = start_datetime;
        this.end = end_datetime;
        this.startPosix = timestamps(1);
        this.endPosix = timestamps(end);
        timelist = [timelist, this];
    end
    fprintf('done with %s\n', char(folders(d)));
    
end

clear d h5files filename timestamps start_datetime end_datetime this f
