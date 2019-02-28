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

function [lun run] = getunbalancedparsCAN(br)
% function [lun run] = getunbalancedparsCAN(br)
lun = [];
% p = [];
run = [];
for i = 1 : length(br)
	if br(i) == '('
		lun = [lun i]; 
	else
		if br(i) == ')'	
			if isempty(lun)
				run = [run i];
			else
% 				p = [p [l(length(l)) i]'];
				lun(length(lun)) = [];
			end
		end
	end
end
% if ~isempty(run), disp([''')''s at ' num2str(run) ' are unmatched.']), end
% br(run) = '.';
% if ~isempty(l), disp(['''(''s at ' num2str(l) ' are unmatched.']), end
% br(l) = '.';
