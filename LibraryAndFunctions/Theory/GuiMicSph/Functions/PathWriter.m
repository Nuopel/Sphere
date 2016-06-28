function [ ] = PathWriter( filename,path )

fileID = fopen(filename,'w');
fprintf(fileID,path);
fclose(fileID);

end

