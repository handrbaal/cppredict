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

function nci = get_nc_bps(sq, br)
% function nci = get_nc_bps(sq, br)
nci=[];
nt1='GCAUGU';
nt2='CGUAUG';
m = br2m(br);
for i=1:size(m,1)
    j=find(m(i,:)==1); 
    if ~isempty(j)
        nt1i = find((nt1 == upper(sq(i))));
        nt2i = find((nt2 == upper(sq(j))));
		if isempty(nt1i) | isempty(nt2i)
			disp(['get_nc_bps: ' sq(i) ', ' sq(j) '.'])
		end
        if ~sum(nt1i == nt2i)
          nci = [nci; [i j]];
        end
    end   
end
