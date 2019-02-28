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

function ss = cut_string(s, delm)
%put substrings of the string s delimited with delm into a cell array ss

s = [delm s delm]; 
i = find(s == delm); 
ss = {};
for ii = 1 : length(i) - 1
  ss{ii} = s(i(ii)+1 : i(ii + 1)-1);
end
  

