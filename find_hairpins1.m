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

function h=find_hairpins1(br)
% function h=find_hairpins1(br), FIND HAIRPINS
i = find(br == '(');
j = find(br == ')');
ij = [i j];
ij_sorted = sort(ij);
% PUTATIVE LOOPS
loops=find(ij_sorted(2:length(ij_sorted)) - ij_sorted(1:length(ij_sorted)-1)>1);
% FROM LOOPS, FIND HAIRPINS
h = [];
for iloops = 1 : length(loops)
% NUCLEOTIDE POSITION OF THE FIRST LEFT PARENTHESIS OF A LOOP 
    li = ij_sorted(loops(iloops));
% NUCLEOTIDE POSITION OF THE FIRST RIGHT PARENTHESIS OF A LOOP
    ri = ij_sorted(loops(iloops)+1);
% IS IT A REAL LOOP...?
    if br(li) == '(' & br(ri) == ')' 
% FIND BOUNDARIES OF AN ENCLOSED ELEMENT (EE), I.E. A STEM LOOP/HAIRPIN FROM ITS LOOP
% THIS PRODUCES EEs WITH THE SAME NUMBER OF NUCLEOTIDE POSITIONS IN STRANDS OF A STEM
        lii = li;
        rii = ri;
        while 1
            if lii > 1 & rii < length(br) 
                if br(lii - 1) ~= ')' & br(rii + 1) ~= '(' 
                    lii = lii - 1; rii = rii + 1; 
                else
                    break
                end
            else
                break
            end
        end

% BALANCE_EEs SO THEY HAVE THE SAME NUMBER OF EXISTING PARENTHESES IN STRANDS OF A STEM. TO USE THIS SET-UP OR NOT...?
        while 1
            il = find(br(lii : li) == '(') + lii - 1;
            ir = find(br(ri : rii) == ')') + ri - 1;
            if length(ir) == length(il), break, end
            if length(ir) > length(il)
% MAKE THE LEFT STRAND OF THE HAIRPIN LONGER
                if lii > 1 & br(lii - 1) == '('
% ADD MORE '('s TO THE LEFT END OF THE HAIRPIN IF THEY ARE AVAILABLE                                        
                    lii = lii - 1;
                else
% REMOVE ')'s FROM THE RIGHT END OF THE HAIRPIN
                    rii = rii - 1; 
                end
            end
            if length(ir) < length(il)
% MAKE THE RIGHT STRAND OF THE HAIRPIN LONGER
                if rii < length(br) & br(rii + 1) == ')'
% ADD MORE ')'s TO THE RIGHT END OF THE HAIRPIN IF THEY ARE AVAILABLE                            
                    rii = rii + 1; 
				else
% REMOVE '('s FROM THE LEFT END OF THE HAIRPIN
					lii = lii + 1;
				end
            end
        end

% JUMP FORWARD OR BACK OVER TRAILING UNPAIRED NUCLEOTIDES AND GAPS AT THE ENDS OF EEs, IF THERE ARE ANY
%         while br(lii) ~= '(' 
%           lii = lii + 1; 
%         end
%         while br(rii) ~= ')'
%           rii = rii - 1;
%         end 
% ANOTHER KIND OF THE ENDS TREATMENT
% 		i = find(br(lii : li) == '(');
% 		j = find(br(lii : li) == '.');
% 		if ~isempty(j) & j(1) < i(1), lii = lii + j(1) - 1; end
% 		i = find(br(ri : rii) == ')');
% 		j = find(br(ri : rii) == '.');
% 		if ~isempty(j) & j(1) > i(1), rii = rii - j(1) + 1; end
		
        h = [h; [lii rii]];
    end
end

% WORK OUT UNAVOIDABLE OVERLAPPING HAIRPINS: THIS SOLUTION SUPPOSES THAT THERE ARE NO PARENTHESES IN THE OVERLAP

	for ih = 1 : size(h, 1) - 1
		if h(ih, 2) >= h(ih + 1, 1)
% FIND POSITIONS OF THE OVERLAP
			d = h(ih, 2) - h(ih + 1, 1);
			h1 = h(ih, 2) - d;
			h2 = h(ih + 1, 1) + d;
% FIND POSITIONS OF 3s IN THE OVERLAP
			i3 = find(br(h1 : h2)=='3') + h1 - 1;
			if ~isempty(i3)
% FIND DISTANCES OF 3s IN THE OVERLAP FROM BEGINNING AND END OF THE OVERLAP 
				d1i3 = abs(i3 - h1);
				d2i3 = abs(i3 - h2);
% WHICH DISTANCES OF 3s ARE BIGGER...?
				ddi3 = d1i3 - d2i3;
% NEGATIVES MEAN CLOSER TO THE BEGINNING OF THE OVERLAP WHEREAS POSITIVE CLOSER TO THE END
				i = find(ddi3 <= 0);
				if ~isempty(i)
% ADD 3s CLOSER TO THE BEGINNING OF THE OVERLAP TO THE END OF THE LEFT HAIRPIN
					h1 = max(i3(i));
% ADD 3s CLOSER TO THE END OF THE OVERLAP TO THE BEGINNING OF THE RIGHT HAIRPIN
					h2 = max(i3(i)) + 1;
				else
% ALL 3s ARE CLOSER TO THE END OF THE OVERLAP
					h2 = min(i3);
					h1 = min(i3) - 1;
				end
			else
				d = (d - 1) / 2;
				h1 = h1 + ceil(d);
				h2 = h2 - floor(d);
			end
			h(ih, 2) = h1;
			h(ih + 1, 1) = h2;
		end
	end

