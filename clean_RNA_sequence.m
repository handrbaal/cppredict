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

function [i sq] = clean_RNA_sequence(sq)
% function sq = clean_RNA_sequence(sq)
sq = upper(sq);
sq(find(sq == 'N')) = 'G';
i = find(sq ~= 'A' & sq ~= 'C' & sq ~= 'G' & sq ~= 'T' & sq ~= 'U');
% if length(i) > 
% 	sq = 0;
% 	return
% end
for ii = 1 : length(i)
	s=baselookup('code',sq(i(ii)));
	s=cut_string(s,char(9));
	s=cut_string(s{1},'|');
	sq(i(ii)) = s{1};
end
sq(find(sq == 'T')) = 'U';

