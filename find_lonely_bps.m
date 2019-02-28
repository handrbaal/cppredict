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

function lp = find_lonely_bps(br)

p = br2inCAN(br);
lp = [];
if size(p, 2) > 1
% 	a = abs(p(1, 2 : length(p)) - p(1, 1 : length(p) - 1));
% 	b = abs(p(2, 2 : length(p)) - p(2, 1 : length(p) - 1));
% 	v = f_un_gr_ve(find(a > 1 & b > 1));
% 	for iv = 1 : length(v)
% 		for iiv = 2 : length(v{iv})
% 			plp = [plp p(:, v{iv}(iiv))];
% 		end
% 	end
	l = strfind(br, '.(.') + 1;
	r = strfind(br, '.).') + 1;
	for il = 1 : length(l)
		i = find(p(1,:)==l(il));
		if ~isempty(i)
			lp = [lp [l(il); p(2, find(p(1,:)==l(il)))]]; 
		end
	end
end
