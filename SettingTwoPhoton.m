function analog = SettingTwoPhoton(fileName)
% here a script to extract all the info from the 2P scope file .MDF
% analog signals, objectives, wavelengths, etc.....
     f = figure(1);%                                                                            
    % Read the .MDF file generated by the 2P MScan Software
    mfile = actxcontrol('MCSX.Data', [0, 0, 500, 500], f);
 %
%%

F = 1;
%for i = 1:size(list,2)
                                                                             
    % Read the .MDF file generated by the 2P MScan Software
    mfile = actxcontrol('MCSX.Data', [0, 0, 500, 500], f);
    %
    openResult = mfile.invoke('OpenMCSFile',([fileName '.MDF']));
    %
    % flag to check if the trial that you want to work with exist!
    if openResult > 0
        close all;
        %clear data;
        warndlg('You are an asshole!')
        return        
    end
    %
    % Extract infos and miscellaneous infos about scanning setting and
    NumOfEpisode      = invoke(mfile, 'ReadParameter', 'episode Count')                       ;
    frameHeight       = invoke(mfile, 'ReadParameter', 'Frame Height')                       ;
    frameWidth        = invoke(mfile, 'ReadParameter', 'Frame Width')                        ;
    frameCount        = invoke(mfile, 'ReadParameter', 'Frame Count')                        ;
    frameInterval     = invoke(mfile, 'ReadParameter', 'Frame Interval (ms)')                ;
    frameDuration     = invoke(mfile, 'ReadParameter', 'Frame Duration (s)')                 ;
    %LaserIntensity    = invoke(mfile, 'ReadParameter', 'Laser Intensity ')                   ;
    LaserWaveLength   = '940nm'                                                                ;
    PixelClock        = invoke(mfile, 'ReadParameter', 'Pixel Clock')                        ;
    %PointsInRegion    = invoke(mfile, 'ReadParameter', 'Points In Region')                   ;
    ScanMode          = invoke(mfile, 'ReadParameter', 'Scan Mode')                          ;
    %ObjectiveLens     = invoke(mfile, 'ReadParameter', 'Objective Lens')                     ;
    MicronsPerPixe    = invoke(mfile, 'ReadParameter', 'Microns Per Pixel')                  ;
    %Objective         = '18X'                                                                ;
    Trigger           = invoke(mfile, 'ReadParameter', 'Imaging Trigger')                    ;
    Magnification     = invoke(mfile, 'ReadParameter', 'Magnification')                      ;
    Created           = invoke(mfile, 'ReadParameter', 'Created On')                         ;
    %Dude              = invoke(mfile, 'ReadParameter', 'Created By')                         ;
    AnalogResolution  = invoke(mfile, 'ReadParameter', 'Analog Resolution')                  ;
    AnalogSampling    = invoke(mfile, 'ReadParameter', 'Analog Acquisition Frequency (Hz)')  ;
    AnalogCh0InRan    = invoke(mfile, 'ReadParameter', 'Analog Ch 0 Input Range')            ;
    AnalogCh0Name     = invoke(mfile, 'ReadParameter', 'Analog Ch 0 Name')                   ;
    AnalogCh1InRan    = invoke(mfile, 'ReadParameter', 'Analog Ch 1 Input Range')            ;
    AnalogCh1Name     = invoke(mfile, 'ReadParameter', 'Analog Ch 1 Name')                   ;
    AnalogCh2InRan    = invoke(mfile, 'ReadParameter', 'Analog Ch 2 Input Range')            ;
    AnalogCh2Name     = invoke(mfile, 'ReadParameter', 'Analog Ch 2 Name')                   ;
    AnalogCh3InRan    = invoke(mfile, 'ReadParameter', 'Analog Ch 3 Input Range')            ;
    AnalogCh3Name     = invoke(mfile, 'ReadParameter', 'Analog Ch 3 Name')                   ;
    AnalogCh4InRan    = invoke(mfile, 'ReadParameter', 'Analog Ch 4 Input Range')            ;
    AnalogCh4Name     = invoke(mfile, 'ReadParameter', 'Analog Ch 4 Name')                   ;
    AnalogCh5InRan    = invoke(mfile, 'ReadParameter', 'Analog Ch 5 Input Range')            ;
    AnalogCh5Name     = invoke(mfile, 'ReadParameter', 'Analog Ch 5 Name')                   ;
    AnalogSam        = invoke(mfile, 'ReadParameter', 'Analog Sample Count')                ;
    TrialLenght       = ceil(str2double(frameDuration(1:8))*str2double(frameCount))          ;
    %invoke(mfile, 'SetEpisode',(5))
    %invoke(mfile, 'ReadParameter', 'Episode Start Time')                                              ;
    Points       = (TrialLenght*(str2double(AnalogSampling(1:5))))-10000;
    
    for j = 1:str2num(NumOfEpisode)
        invoke(mfile, 'SetEpisode',(j));
        for h = 1:6
            analog(j).data(:,h) = double(invoke(mfile, 'ReadAnalog',h, str2num(AnalogSam) ,0))/2^15 *10;
        end
        analog(j).startTime = invoke(mfile, 'ReadParameter', 'Episode Start Time');
    end
close all;
