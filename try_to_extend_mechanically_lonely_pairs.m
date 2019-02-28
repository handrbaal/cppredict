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

function br = try_to_extend_mechanically_lonely_pairs(sq, br)
% function br = try_to_extend_mechanically_lonely_pairs(sq, br)

lp = find_lonely_bps(br);
for ip = 1 : size(lp, 2)
% TRY TO EXTEND A LONELY PAIR TO OUTSIDE 
	i = 1;
	while lp(1, ip) - i & lp(2, ip) + i <= length(br) & br(lp(1, ip) - i) == '.' & br(lp(2, ip) + i) == '.'
		if (sq(lp(1, ip) - i) == 'G' & sq(lp(2, ip) + i) == 'C') | (sq(lp(1, ip) - i) == 'C' & sq(lp(2, ip) + i) == 'G') | ...
		   (sq(lp(1, ip) - i) == 'A' & sq(lp(2, ip) + i) == 'U') | (sq(lp(1, ip) - i) == 'U' & sq(lp(2, ip) + i) == 'A') | ...
		   (sq(lp(1, ip) - i) == 'G' & sq(lp(2, ip) + i) == 'U') | (sq(lp(1, ip) - i) == 'U' & sq(lp(2, ip) + i) == 'G')
			br([lp(1, ip) - i lp(2, ip) + i]) = '()';
% IF THIS ENABLED, ADD ONLY NEIGBOURING POTENTIAL PAIRS; INTRODUCE NO UNPAIRED NUCLEOTIDES BETWEEN ORIGINAL AND NEW BPS
% 		else
% 			break
		end
 		i = i + 1;
	end
	
% TRY TO EXTEND A LONELY PAIR INTO INSIDE 
	i = 1;
	while abs((lp(1, ip) + i) - (lp(2, ip) - i)) > 3 & br(lp(1, ip) + i) == '.' &  br(lp(2, ip) - i) == '.'
		if (sq(lp(1, ip) + i) == 'G' & sq(lp(2, ip) - i) == 'C') | (sq(lp(1, ip) + i) == 'C' & sq(lp(2, ip) - i) == 'G') | ...
		   (sq(lp(1, ip) + i) == 'A' & sq(lp(2, ip) - i) == 'U') | (sq(lp(1, ip) + i) == 'U' & sq(lp(2, ip) - i) == 'A') | ...
		   (sq(lp(1, ip) + i) == 'G' & sq(lp(2, ip) - i) == 'U') | (sq(lp(1, ip) + i) == 'U' & sq(lp(2, ip) - i) == 'G')
			br([lp(1, ip) + i lp(2, ip) - i]) = '()';
% IF THIS ENABLED, ADD ONLY NEIGBOURING POTENTIAL PAIRS; INTRODUCE NO UNPAIRED NUCLEOTIDES BETWEEN ORIGINAL AND NEW BPS
% 		else
% 			break
		end
		i = i + 1;
	end
end
