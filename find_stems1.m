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

function stems = find_stems1(br, hairpins)
% function stems = find_stems1(br, hairpins)

stems = [];
 
if isempty(hairpins), return, end

% MARK HAIRPINS BY 'x'S
for i = 1 : size(hairpins, 1), br(hairpins(i, 1) : hairpins(i, 2)) = 'x'; end

x = find(br == 'x');
p = br2inCAN(br);
if ~isempty(p)
% FIND THOSE PAIRS IN p WITH NO 'x'S BETWEEN THEM IN NONE OF THE 2 STRANDS. SUCH PAIRS FORM STEMS.
	stem = 1;
	for ip = 2 : size(p, 2)
		if p(1, ip - 1) >= p(1, ip)
			pp = p(1, [ip ip - 1]); 
		else
			pp = p(1, [ip - 1, ip]); 
		end
		if ~any(x >= pp(1) & x <= pp(2)) & ~any(x >= p(2, ip - 1) & x <= p(2, ip))
			stem = [stem ip];
		else
			stems = [stems; [p(1, stem(1)) p(1, stem(length(stem))) p(2, stem(1)) p(2, stem(length(stem)))]];
			stem = ip;
		end
	end
	stems = [stems; [p(1, stem(1)) p(1, stem(length(stem))) p(2, stem(1)) p(2, stem(length(stem)))]];
	i = find(stems(:, 1) > stems(:, 2));
	stems(i, [1 2]) = stems(i, [2 1]);

% ADD TRAILING UNPAIRED NUCLEOTIDES

	% for istems = 1 : size(stems, 1)
	% 	stems(istems, 1) = stems(istems, 1) - 1;
	% 	while stems(istems, 1) & br(stems(istems, 1)) ~= 'x' & br(stems(istems, 1)) ~= '(' %& br(stems(istems, 1)) ~= '-' 
	% 	   stems(istems, 1) = stems(istems, 1) - 1; 
	% 	end
	% 	stems(istems, 2) = stems(istems, 2) + 1;
	% 	while stems(istems, 2) <= length(br) & br(stems(istems, 2)) ~= 'x' & br(stems(istems, 2)) ~= '(' %& br(stems(istems, 2)) ~= '-' 
	% 	   stems(istems, 2) = stems(istems, 2) + 1; 
	% 	end
	% 	stems(istems, 3) = stems(istems, 3) - 1;
	% 	while stems(istems, 3) & br(stems(istems, 3)) ~= 'x' & br(stems(istems, 3)) ~= ')' %& br(stems(istems, 3)) ~= '-'
	% 	   stems(istems, 3) = stems(istems, 3) - 1; 
	% 	end
	% 	stems(istems, 4) = stems(istems, 4) + 1;
	% 	while stems(istems, 4) <= length(br) & br(stems(istems, 4)) ~= 'x' & br(stems(istems, 4)) ~= ')' %& br(stems(istems, 4)) ~= '-'
	% 	   stems(istems, 4) = stems(istems, 4) + 1; 
	% 	end
	% 	
	% 	stems(istems, :) = stems(istems, :) + [1 -1 1 -1];
	% 	
	% 	br([stems(istems, 1) : stems(istems, 2) stems(istems, 3) : stems(istems, 4)]) = 'x';
	% 
	% 	if ~sum(br=='(' | br==')'), break, end
	% end
	% clear i x p stem pp ip l c br 
end
