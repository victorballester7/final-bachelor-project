function [TLE] = getSatelliteTLE(ID,inEpochDatenum)
%% *purpose*
% return the TLE for Satellite based on epoch
%% *inputs*
%  ID - Spacecraft ID
%       73027 = Skylab
%  inEpochDatenum - TLEs can be selected based on epoch 
%% *outputs*
%  TLE - the two line element set corresponding to the satellite at that
%        epoch
%% *history*
%  When       Who    What
%  ---------- ------ --------------------------------------------------
%  2019/07/17 mnoah  original code
%  2020/01/19 mnoah  placeholder - edit to put your own TLE parser
%                    depending on your TLE source

if (ID == 25544)
    TLE = { ...
        % 'ISS (ZARYA)'; ...
        % '1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927'; ...
        % '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'};
        'NUTSAT'; ...
        '1 55124U 98067UR  23086.16778994  .01540407  14102-2  26634-2 0  9991'; ...
        '2 55124  51.6237   5.0017 0010757 188.3869 171.6963 16.01878199 13792'};
else
    error('ISS Placeholder only - modify code for getting TLE of other spacecraft');
end

end




