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

function flag = check_fasta(fn_)

sqs = fastaread(fn_);
BAD_CHARS = ':/''()"';
flag = 1;
for isqs = 1 : length(sqs)
% CHECK SEQUENCE
	j = find(upper(sqs(isqs).Sequence) ~= 'A' & upper(sqs(isqs).Sequence) ~= 'C' & upper(sqs(isqs).Sequence) ~= 'G' & upper(sqs(isqs).Sequence)~= 'U');
	if ~isempty(j)
		disp(['check_fasta.m: Non-sequence chars in ''' sqs(isqs).Header ''' (' sqs(isqs).Sequence(j) ') in ' fn_ '.'])
		flag = 0;
		return
	end
% CHECK HEADER
	j = fvi(sqs(isqs).Header, BAD_CHARS);
	if ~isempty(j)
		disp(['check_fasta.m: Bad chars (' sqs(isqs).Header(j) ')in ''' sqs(isqs).Header ''' header in ' fn_ '.'])
		flag = 0;
		return
	end
% CHECK SPECIES
	j = strmatch(sqs(isqs).Header, {sqs([1 : isqs - 1 isqs + 1 : length(sqs)]).Header}, 'exact');
	if ~isempty(j)
		disp(['check_fasta.m: Ambiguous sequence headers for ' sqs(isqs).Header ' in ' fn_ '.'])
		flag = 0;
		return
	end
end


