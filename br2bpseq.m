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

function [fn_bpseq,bpseq] = br2bpseq(fn_br)
% function fn_bpseq = br2bpseq(fn_br)
id = fopen(fn_br); br=textscan(id, '%s','delimiter',char(10));fclose(id);
if length(br{1}) ~= 3, disp(['br2bpseq.m: Not enough lines in ' fn_br '.']), pause, end
m=br2m(br{1}{3});
fn_bpseq = [fn_br(1 : strfind(fn_br, '.br')) 'bpseq'];
bpseq=zeros(size(m,1), 2);
for i=1:size(m,1)
    j=find(m(i,:)==1); 
    if ~isempty(j)
        bpseq(i, :) = [i j]; 
        bpseq(j, :) = [j i]; 
    end   
end
bpseq(:, 1) = 1 : size(m,1);
id = fopen(fn_bpseq, 'w');
for i=1:size(m,1)
    fwrite(id, [num2str(bpseq(i, 1)) ' ' br{1}{2}(i) ' ' num2str(bpseq(i, 2)) char(10)]);
end
fclose(id);
