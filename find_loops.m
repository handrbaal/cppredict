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

function loops=find_loops(br)
% function loops=loops(br), FIND LOOPS, I.E. '(.)'
i = find(br == '(');
j = find(br == ')');
ij = [i j];
ij_sorted = sort(ij);
% PUTATIVE LOOPS
l=find(ij_sorted(2:length(ij_sorted)) - ij_sorted(1:length(ij_sorted)-1)>1);
% GET REAL LOOPS
loops = [];
for il = 1 : length(l)
% IS IT A REAL LOOP...? CHECK THE FIRST LEFT PARENTHESIS OF A PUTATIVE 
% LOOP AND THE FIRST RIGHT PARENTHESIS OF A PUTATIVE LOOP
    if br(ij_sorted(l(il))) == '(' & br(ij_sorted(l(il)+1)) == ')' 
        loops = [loops; [ij_sorted(l(il)) ij_sorted(l(il)+1)]];
    end
end
