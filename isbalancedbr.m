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

function f = isbalancedbr(br)
if isempty(br), f = 0; return, end
% f = 0;
% if (sum(br=='[') == sum(br==']')) & (sum(br=='{') == sum(br=='}')) & (sum(br=='<') == sum(br=='>'))
%     f = 1;
% end
f=1;
l = [];
run = [];
for i = 1 : length(br)
	if br(i) == '('
		l = [l i]; 
	else
		if br(i) == ')'	
			if isempty(l)
				run = [run i];
			else
				l(length(l)) = [];
			end
		end
	end
end
if ~isempty(run) | ~isempty(l), f = 0; return, end

l = [];
run = [];
for i = 1 : length(br)
	if br(i) == '['
		l = [l i]; 
	else
		if br(i) == ']'
			if isempty(l)
				run = [run i];
			else
				l(length(l)) = [];
			end
		end
	end
end
if ~isempty(run) | ~isempty(l), f = 0; return, end

l = [];
run = [];
for i = 1 : length(br)
	if br(i) == '{'
		l = [l i]; 
	else
		if br(i) == '}'
			if isempty(l)
				run = [run i];
			else
				l(length(l)) = [];
			end
		end
	end
end
if ~isempty(run) | ~isempty(l), f = 0; return, end

l = [];
run = [];
for i = 1 : length(br)
	if br(i) == '<'
		l = [l i]; 
	else
		if br(i) == '>'	
			if isempty(l)
				run = [run i];
			else
				l(length(l)) = [];
			end
		end
	end
end
if ~isempty(run) | ~isempty(l), f = 0; return, end
