function [output] = getUniqueName(str)
    currentDateTime = datetime('now');
    currentDateTimeString = char(currentDateTime);
    currentDateTimeString = strrep(currentDateTimeString, ':', '_');
    currentDateTimeString = strrep(currentDateTimeString, ' ', '_');
    output = strcat('data/',currentDateTimeString,'__',str);
end