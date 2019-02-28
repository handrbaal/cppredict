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

function vi = fvi(v1, v2)
% function vi = fvi(v1, v2)
% Search for v2 values in the v1. vi contains v1 indexes of v2 values found in v1.

vi = [];
for i = 1 : length(v2)
	vi = [vi find(v1 == v2(i))];
end

% for i = 1 : length(v2)
%  j = find(v1 == v2(i));
%  if ~isempty(j)
%   vi = [vi j(1)];
%  end
% end
