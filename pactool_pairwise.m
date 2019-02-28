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

if RAND
% 	rng('shuffle'); 	
	rng('default');
	fnrn = ['rn' num2str(rand) '.fa'];
	id=fopen(fnrn,'w'); 
	if exist('TEST_EACH') & TEST_EACH
% TEST WITH RANDOM TARGET SEQUENCES. MEAN LENGTH OF SEQUENCES IN RF00013.fasta IS 179.9865 NUCLEOTIDES.
% RE-READ NAME OF FILE WITH TARGET SQS AFTER RAND COMPUTATION
		sqs = fastaread(fn__);
		for i = 1 : length(sqs)
% 	 		fwrite(id, ['>RAND_' num2str(i) char(10) randseq(length(sqs(i).Sequence),'alphabet','rna') char(10)]); 
			fwrite(id, ['>RAND_' num2str(i) char(10) sqs(i).Sequence(randperm(length(sqs(i).Sequence))) char(10)]);
% 			eval(['!ushuffle -s ' sqs(i).Sequence ' -n ' num2str(RANDSQS) ' -k 2 > ushuffle.fasta'])
% 			s = fopen('ushuffle.fasta'); r = textscan(s, '%s', 'delimiter', char(10)); fclose(s);
% 			for ir = 1 : length(r{1})
% 				fwrite(id, ['>RAND_' num2str(ir) char(10) r{1}{ir} char(10)]);
% 			end
% 			!rm ushuffle.fasta
		end
	end
	if exist('TEST_SOURCE') & TEST_SOURCE
% TEST WITH RANDOM SOURCE SEQUENCE.
		s = fopen(sr__fn_br_); t = textscan(s, '%s', 'delimiter', char(10)); fclose(s);
		
		for isqs = 1 : RANDSQS
			fwrite(id, ['>RAND_' num2str(isqs) char(10) t{1}{2}(randperm(length(t{1}{2}))) char(10)]);
% 			fwrite(id, ['>RAND_' num2str(i) char(10) randseq(length(t{1}{2}),'alphabet','rna') char(10)]); 
		end

% 		eval(['!ushuffle -s ' t{1}{2} ' -n ' num2str(RANDSQS) ' -k 2 > ushuffle.fasta'])
% 		s = fopen('ushuffle.fasta'); r = textscan(s, '%s', 'delimiter', char(10)); fclose(s);
% 		for ir = 1 : length(r{1})
% 			fwrite(id, ['>RAND_' num2str(ir) char(10) r{1}{ir} char(10)]);
% 		end
% 		!rm ushuffle.fasta
	end
	fclose(id);
	fn__ = fnrn;
	clear fnrn
end

sqs_ = fastaread(fn__);

if ~RAND, sqs_original_ = sqs_; end

id=fopen(sr__fn_br_); t = textscan(id, '%s', 'delimiter', char(10)); fclose(id);
sr__h_ = t{1}{1}(2 : length(t{1}{1}));
sr__q_ = t{1}{2};
sr__br_ = t{1}{3};
if getunbalancedparsCAN(sr__br_), disp('Source structure is unbalanced!'), return, end
id=fopen(sr__fn_bpseq_); t = textscan(id, '%u%s%u', 'delimiter', ' '); fclose(id);
sr__bpseq_ = [t{1} t{3}];
if length(sr__br_) ~= size(sr__bpseq_, 1), disp([char(10) 'pactool_pairwise.m: length(sr__br_) ~= size(sr__bpseq_, 1)' char(10)]), pause, return, end

sc_ = ones(1, length(sqs_)) .* -1;
ldiff_ = zeros(1, length(sqs_));
sc_sq_ = sc_;
sq_l_rat_ = ldiff_;

fn_ = ['1.fasta.aln']; 

ind = 1:length(sqs_);
% ind = 173;
if ~RAND, ind_ = ind; what = 'Doing '; else, what = 'Bootstrap ';  end
for iind = 1 : length(ind) 
    isqs_ = ind(iind);
	disp([what num2str(iind) ' / ' num2str(length(ind)) '.'])
	
	id = fopen(sr__fn_br_); t = textscan(id, '%s', 'delimiter', char(10)); fclose(id);
	id = fopen('1.fasta', 'w');	fwrite(id, [t{1}{1} char(10) t{1}{2} char(10) '>' sqs_(isqs_).Header char(10) sqs_(isqs_).Sequence char(10)]); fclose(id); 
% 	eval(['!head -n2 ' sr__fn_br_ ' > 1.fasta']); eval(['!echo -e ''>' sqs_(isqs_).Header '\n' sqs_(isqs_).Sequence ''' | cat >> 1.fasta'])
	
    switch ALM
		case 'clustalw2'
			eval(['!clustalw2 -GAPOPEN=7 -GAPEXT=0.5 -infile=1.fasta -outfile=' fn_ ' -outorder=input -output=fasta 2>&1 > NULL'])
%			eval(['!clustalw2 -infile=1.fasta -outfile=' fn_ ' -outorder=input -output=fasta 2>&1 > NULL'])
		case 'muscle'
			eval(['!~/muscle3.8.31_i86linux64 -in 1.fasta -out ' fn_ ' -quiet'])
		case 'clustalo'
			eval(['!clustalo --infile=1.fasta --outfile=' fn_ ' --output-order=input-order --outfmt=fasta --force 2>&1 > NULL']) 
		otherwise
			disp(['Disallowed value of ALM parameter (' ALM ').'])
    end
    reformat_fasta
    aln_=fastaread(fn_);
%     gaps, copy__pairwise_gaps
    copy__pairwise_hard
	
	if EXTEND_MECHANICALLY_LONELY_PAIRS
		tg_br_ = try_to_extend_mechanically_lonely_pairs(sqs_(isqs_).Sequence, tg_br_);
	end
	
    if ~isbalancedbr(tg_br_)
		disp([sqs_(isqs_).Header '''s not balanced! Unbalanced pairs thrown out.'])
		[lun run] = getunbalancedparsCAN(tg_br_);
		if isempty(lun) & isempty(run), disp('But this is weird!'), return, end
		tg_br_([lun run]) = '.';
	end
    
% SAVE PREDICTED BR INTO A FILE
	
% 	eval(['!echo -e ''>' sqs_(isqs_).Header '\n' sqs_(isqs_).Sequence '\n' tg_br_ ''' | cat > ''' sqs_(isqs_).Header '.br'''])
	id = fopen([sqs_(isqs_).Header '.br'], 'w'); 
	fwrite(id, ['>' sqs_(isqs_).Header char(10) sqs_(isqs_).Sequence char(10) tg_br_ char(10)]); 
	fclose(id);

% SHOW THE PREDICTED STRUCTURE
    if SHOW_PREDICTED_STRUCTURES_INDIVIDUALLY
        s = [];
        if exist('ref__') & ~isempty('ref__')
            [~, str] = unix(['sed -n "3p" < ' ref__]);
            str(find(str == char(10))) = [];
            if length(str) == length(tg_br_)
                diff_ = find(str ~= tg_br_);
                s = ['-basesStyle1 outline=#FF0000 -applyBasesStyle1on "' num2str(diff_(1))];
                for isqs=2:length(diff_), s = [s ', ' num2str(diff_(isqs))]; end
                s = [s '"'];
            else
                disp(['Length of ref and target structure does not agree.'])
            end
        else
            if ~isempty(diff_)
                s = ['-basesStyle1 outline=#FF0000 -applyBasesStyle1on "' num2str(diff_(1))];
                for isqs=2:length(diff_), s = [s ', ' num2str(diff_(isqs))]; end
                s = [s '"'];
            end
        end
        eval(['!java -cp ../VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -algorithm radial -i ' sqs_(isqs_).Header '.br ' s ' &'])
	end
    
% GET SOME INFORMATIVE NUMBERS: IT IS NOT PRECISION EVALUATION! IT JUST SAYS HOW MANY POSITIONS DIFFER in SR AND TG BRs, WHEN THEY ARE ALIGNED.

	sq_l_rat_(isqs_) = length(sqs_(isqs_).Sequence) / length(sr__q_);
    s = [sr__h_ ' -> ' sqs_(isqs_).Header '. Aligned target and source differ in ' num2str(round((100*2*length(diff_)/(length(sr__q_) + length(sqs_(isqs_).Sequence))) * 10)/10) '% of br. (Seq. length ratio=' num2str(sq_l_rat_(isqs_)) ').'];
    ldiff_(isqs_) = abs(length(sqs_(isqs_).Sequence) - length(sr__q_));
	sc_sq_(isqs_) = round((sum((aln_(1).Sequence == aln_(2).Sequence) & aln_(1).Sequence ~= '-' & aln_(2).Sequence ~= '-') / ((length(sr__q_) + length(sqs_(isqs_).Sequence))/2))*100);
    
% MATCH TARGET TO A REFERENCE STRUCTURE, IF THE LATTER EXISTS OR IF IS USED

    if exist('ref__') & ~isempty('ref__')
        id=fopen(ref__); t = textscan(id, '%s', 'delimiter', char(10)); fclose(id);
        if length(t{1}{3}) == length(tg_br_)
            s = [s ' Acc. to ref.=' num2str((round((sum(tg_br_ == t{1}{3}) / length(t{1}{3})) * 1000))/10) '%.'];
            if iscanstr(sqs_(isqs_).Sequence, tg_br_) & iscanstr(t{1}{2}, t{1}{3})
                eval(['!cat ' ref__ ' > 1.str && cat ''' sqs_(isqs_).Header '.br'' >> 1.str'])
                eval(['!RNAdistance < 1.str > 1.str.dist'])
                [~, ss] = unix(['grep ''f:'' 1.str.dist']);
                sc_(isqs_) = str2num(ss(find(ss == ':') + 1 : length(ss)));
                s = [s ' Sim. to ref__: ' num2str(sc_(isqs_)) '.'];
            end
        else
            disp(['Length of ref and target structure does not agree.'])
        end
    end
    
% MATCH TARGET TO A REFERENCE STRUCTURE, WHICH CAN BE A WHATEVER STRUCTURE YOU WANNA MATCH TO THE TARGET
    
    if exist('reference_structure_template_') & ~isempty(reference_structure_template_)
		id = fopen(reference_structure_template_); t = textscan(id, '%s', 'delimiter',char(10)); fclose(id);
        if isbalancedbr(tg_br_) & iscanstr(sqs_(isqs_).Sequence, tg_br_) & iscanstr(t{1}{2}, t{1}{3})
            eval(['!cat ' reference_structure_template_ ' > 1.str && cat ''' sqs_(isqs_).Header '.br'' >> 1.str'])
            eval(['!RNAdistance < 1.str > 1.str.dist'])
            [~, ss] = unix(['grep ''f:'' 1.str.dist']);
			if ~isempty(ss)
				sc_(isqs_) = str2num(ss(find(ss == ':') + 1 : length(ss)));
				s = [s ' Sim. to templ.: ' num2str(sc_(isqs_)) '.'];
			end
        end
    end
    if VERBOSE 
		disp(s)
	end
    if PAUSE_BETWEEN_INDIVIDUAL_STRUCTURES
        pause
        if CLOSE_PREDICTED_STRUCTURES_INDIVIDUALLY
            [~, s]=unix('ps -fu pepik | grep -basesStyle1');
            s=strsplit(s,'\n');
            if ~isempty(strfind(s{1}, 'VARNA'))
                s=strsplit(s{1},' ');
                eval(['!kill ' s{2}])
            end
        end
    end
end
!rm -f alm.br 2.fasta* NULL 1.fasta.aln 1.fasta 1.dnd 1.str 1.str.dist rna.ps 
if RAND, eval(['!rm ' fn__]), end
clear isqs isqs_ id di t fn_ aln_ diff_ i s iind ss str sq sqs ind tg_br_aln_ tg_br_ idout sqs_ what 
