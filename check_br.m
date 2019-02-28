% Copyright 2018 Josef Pánek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function flag = check_br(fn)

id=fopen(fn); t = textscan(id, '%s', 'delimiter', char(10)); fclose(id);
BAD_CHARS = ':/''()"';
flag = 1;
for ibr = 1 : 3 : length(t{1})
	if ibr + 2 <= length(t{1}) & ~isempty(t{1}{ibr}) & ~isempty(t{1}{ibr + 1}) & ~isempty(t{1}{ibr + 2})
	% CHECK SEQUENCE
		j = find(upper(t{1}{ibr + 1}) ~= 'A' & upper(t{1}{ibr + 1}) ~= 'C' & upper(t{1}{ibr + 1}) ~= 'G' & upper(t{1}{ibr + 1})~= 'U');
		if ~isempty(j)
			disp(['Non-sequence chars in ''' t{1}{ibr} ''' (' t{1}{ibr + 1}(j) ') in ' fn '.'])
			flag = 0;
			return
		end
	% CHECK HEADER
		j = fvi(t{1}{ibr}, BAD_CHARS);
		if ~isempty(j)
			disp(['Bad chars (' t{1}{ibr}(j) ')in ''' t{1}{ibr} ''' header in ' fn '.'])
			flag = 0;
			return
		end
	% CHECK SPECIES
		j = strmatch(t{1}{ibr}, t{1}([1 : ibr - 1 ibr + 1 : length(t{1})]), 'exact');
		if ~isempty(j)
			disp(['Ambiguous sequences: ' t{1}{ibr} ' in ' fn '.'])
			flag = 0;
			return
		end
	% CHECK STRUCTURE
		j = find(t{1}{ibr + 2} ~= '.' & t{1}{ibr + 2} ~= '(' & t{1}{ibr + 2} ~= ')');
		if ~isempty(j)
			disp(['Non-structure chars in ''' t{1}{ibr} ''' (''' t{1}{ibr + 2}(j) ''') in ''' fn '''.'])
			flag = 0;
			return
		end	
	end
end


