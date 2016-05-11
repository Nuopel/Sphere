function [ VECTOR ] = ResizeLine( VECTOR )
% Check VECTOR and if needed change is dimension from column to line

[a, b ] = size(VECTOR);
if b>a
    VECTOR = VECTOR.';
end

end

