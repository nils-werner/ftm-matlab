function [ result ] = findfigure( name )
%FINDFIGURE Finds a figure for a given name. Creates one if none is found.
%   Detailed explanation goes here

	if(findobj('type','figure','name',name))
		result = figure(findobj('type','figure','name',name));
	else
		result = figure('name',name);
	end

end

