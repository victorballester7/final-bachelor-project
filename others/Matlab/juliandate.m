function jd = juliandate(year, month, day, hour, minute, second, action)
% Calculate Julian date from year, month, day, hour, minute, and second
% Usage: jd = juliandate(year, month, day, hour, minute, second, 'action', value)

if nargin < 7
    action = 'none'; % Default action is 'none'
end

% Calculate the decimal day
decimal_day = day + (hour + (minute + second / 60) / 60) / 24;

% Calculate the Julian date for the given date
if month <= 2
    year = year - 1;
    month = month + 12;
end

A = floor(year / 100);
B = 2 - A + floor(A / 4);
jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + decimal_day + B - 1524.5;

% Apply the specified action, if any
switch lower(action)
    case 'add'
        jd = jd + value;
    case 'subtract'
        jd = jd - value;
    otherwise
        % Do nothing
end

end
