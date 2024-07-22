function [capitalized] = capitalize(name)
%CAPITALIZE Change the first letter of the provided string to uppercase
    capitalized = name;
    capitalized(1) = upper(capitalized(1));
end

