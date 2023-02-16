%% comma2point.m
% replaces all occurences of comma (",") with point (".") in a text-file.
% Note that the file is overwritten, which is the price for high speed.

function comma2point(filespec)

    file = memmapfile( filespec, 'writable', true );
    comma = uint8(',');
    point = uint8('.');
    file.Data( transpose( file.Data==comma) ) = point;
    delete(file)
end