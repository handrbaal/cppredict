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

sqs = fastaread(fn_);
id = fopen([fn_ '.cleaned'],'w');

for isqs = 1 : length(sqs)
% CLEAN SEQUENCE
	sqs(isqs).Sequence = upper(sqs(isqs).Sequence);

	[j sq] = clean_RNA_sequence(sqs(isqs).Sequence);
	if ~isempty(j)
		disp(['Non-sequence chars in ''' sqs(isqs).Header ''' (' sqs(isqs).Sequence(j) ') found and cleaned.'])
		sqs(isqs).Sequence = sq;
	end

	fwrite(id, ['>' sqs(isqs).Header char(10) sqs(isqs).Sequence char(10)]);
end

fclose(id);
clear id i j sq isqs sqs 
