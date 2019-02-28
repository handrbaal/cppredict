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

% YOU DON'T NEED gaps.m HERE 'COS YOU TAKE ALL GAPS

% INPUT PARS: aln_, sr__q_, sr__br_, sr__bpseq_, OUTPUT PARS: tg_br_, diff_

% eval(['!echo ''>''' aln_(1).Header ' | cat > alm.br && echo ' aln_(1).Sequence ' | cat >> alm.br'])
% REPLACE BAD CHARACTERS IN SOURCE BR
sr__br_(find(sr__br_ == '-')) = '.';

% MARK GAPS IN SOURCE
n2asr = find(aln_(1).Sequence ~= '-');
sr_br_aln = char(ones(1, length(aln_(1).Sequence)) .* '0');
sr_br_aln(n2asr) = sr__br_;
% eval(['!echo ''' sr_br_aln ''' | cat >> alm.br'])

% eval(['!echo ''>''' aln_(2).Header ' | cat >> alm.br && echo ' aln_(2).Sequence ' | cat >> alm.br'])
tg_br_aln_ = sr_br_aln;
% % eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])
% MARK GAPPED POSITIONS IN TARGET 
tg_br_aln_(find(aln_(2).Sequence == '-')) = '-';
% % eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])
% AMONG TARGET GAPS, MARK THOSE THAT ARE UNPAIRED IN SOURCE
% tg_br_aln_(find(sr_br_aln == '.' & tg_br_aln_ == '-')) = '1';
% % eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])

% '0' - GAPS IN SOURCE, I.E. NO CORRESPONDING NUCLEOTIDES IN SOURCE TO THOSE IN TARGET. NOTHING (I.E. NO 
% PAIRS)TO COPY.
% '-' - GAPS IN TARGET. NO CORRESPONDING NUCLEOTIDES IN TARGET TO THOSE IN SOURCE. WHEN COPY THE SOURCE PAIRS, IT 
% CAN BE ONE OF THE FOLLOWING:
% '1' - ONE OF THE SOURCE PAIRING NUCLEOTIDES MAPS INTO A GAP IN TARGET,
% '2' - BOTH SOURCE PAIRING NUCLEOTIDES MAP INTO GAPS IN TARGET.

% FIND GAPS IN THE ALIGNED TARGET SEQUENCE
i = find(tg_br_aln_ == '-');
% FIND THESE GAPS IN THE UNALIGNED SOURCE SEQUENCE
j = fvi(n2asr, i);
% FIND NUCLEOTIDES THAT PAIR WITH THESE GAPS
k = sr__bpseq_(j, 2);
% FIND THESE GAP-PAIRED NUCLEOTIDES IN THE ALIGNED BOTH SOURCE AND TARGET SEQUENCES
l = n2asr(k(find(k)));
% MARK THESE GAP-PAIRED NUCLEOTIDES IN THE ALIGNED TARGET SEQUENCE
tg_br_aln_(l) = '1';
% % eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])
% MARK THOSE OUT OF THESE GAP-PAIRED NUCLEOTIDES THAT ALREADY WERE GAPS, I.E. BOTH SOURCE PAIRING 
% NUCLEOTIDES MAP INTO GAPS IN TARGET
tg_br_aln_(find(aln_(2).Sequence == '-' & tg_br_aln_ == '1')) = '2';
% % eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])

% THROW OUT AND MARK BY '3' NON-CANONICAL BASE PAIRS. THIS WAS MADE ORIGINALLY FOR 6S RNA THAT HASN'T GOT 
% A CONSERVED CORE AND THEREFORE ITS SEQUENCE ALIGNMENT DOESN'T CORRESPOND TO A STRUCTURE ALIGNMENT. THEN
% BASE PAIRS ARE COPIED FROM SOURCE TO WRONG NUCLEOTIDES IN TARGET, PRODUCING NONSENSE NON-CANONICAL 
% BASE PAIRS. CAN BE USED FOR STRUCTURES WITHOUT NON-CANONICAL BASE PAIRS ONLY AS YOU ARE NOT ABLE TO DISTINGUISH
% RIGHT NON-CANONICAL BASE PAIRS FROM NONSENSE, I.E. WRONG NON-CANONICAL BASE PAIRS.
if THROW_OUT_COPIED_NC_BPS & iscanstr(sr__q_, sr__br_)
	if ~isbalancedbr(tg_br_aln_)
		[lun run] = getunbalancedparsCAN(tg_br_aln_);
		tg_br_aln_([lun run]) = '.';
	end
	nci = get_nc_bps(aln_(2).Sequence, tg_br_aln_);
	tg_br_aln_(nci(:)) = '3';
end
% eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])

if PASTE
	paste
% 	eval(paste_proc___);
end

% eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])

% FIND DIFFERENCES BETWEEN SOURCE AND TARGET IN ALIGNED STRUCTURES
tg_br_aln_(find(tg_br_aln_ == '0' | tg_br_aln_ == '1' | tg_br_aln_ == '3')) = '.';
% eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])
diff_ = find(sr_br_aln ~= tg_br_aln_);
j = find(tg_br_aln_ == '2' | tg_br_aln_ == '-');
diff_(fvi(diff_, j)) = [];
n2asr = find(aln_(2).Sequence ~= '-');
diff_ = fvi(n2asr, diff_);

% PRODUCE UNALIGNED TARGET BR
tg_br_ = tg_br_aln_;
tg_br_(j) = [];

clear i j k l n2asr sr_br_aln tg_br_aln nt1 nt2 nt1i nt2i m 
