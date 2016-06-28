function [ VECTOR ] = ResizeColumn( VECTOR )
% Check VECTOR and if needed change is dimension from line to column

[a, b ] = size(VECTOR);
if b>a
    VECTOR = VECTOR.';
end

end

