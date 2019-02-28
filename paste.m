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

% SOMETIMES IF YOU PREDICT EVERY SHIT, I.E. HAIRPINS AND STEMS WITH SMALLEST INCONSISTENCY, E.G. WITH ONE INCONSISTENT NUCLEOTIDE ONLY, IT IS CONTRAPRODUCTIVE 'COS YOU ADD 
% MORE BASE PAIRS BY RNAFOLD PREDICTION THAN NEEDED
% HAIRPIN_INCONSISTENCE_TRESHOLD_RNAFOLD = 0;
% STEM_INCONSISTENCE_TRESHOLD_RNAFOLD = 0;

if ~exist('USE_RNAFOLD_WITH_C'), USE_RNAFOLD_WITH_C = 0; end
% if ~exist('COMBINE_COPIED_AND_RNAFOLD_BPS'), COMBINE_COPIED_AND_RNAFOLD_BPS = 0; end
% THE USE OF USE_RNAFOLD_WITH_C=0 & COMBINE_COPIED_AND_RNAFOLD_BPS=1 produces imbalanced targets. 
if ~USE_RNAFOLD_WITH_C, COMBINE_COPIED_AND_RNAFOLD_BPS = 0; end
% RNAfold -C (I.E. USE_RNAFOLD_WITH_C=1) can miss trailing pairs. If COMBINED WITH COMBINE_COPIED_AND_RNAFOLD_BPS = 1, they are put back, if they were copied before 
% RNAfold prediction.
if USE_RNAFOLD_WITH_C, COMBINE_COPIED_AND_RNAFOLD_BPS = 1; end 

hairpins = find_hairpins1(tg_br_aln_); 
for ihairpins = 1 : size(hairpins, 1)
    str = tg_br_aln_(hairpins(ihairpins, 1) : hairpins(ihairpins, 2));

% % THIS IS TO FIX OVERLAPPING HAIRPINS, WHEN THE PREVIOUS ONE WAS ALREADY PREDICTED AND TH EPREDICTED BPS ARE IN THE OVERLAP OF THE NEXT HAIRPIN AS UNBALANCED BPS. THAT'S WHY THE UNBALANCED BPS ARE ALWAYS ON THE LEFT OF THE NEXT HAIRPIN. THEY CAUSE ONBALNCED CONSTRAINT WHEN THEY GO THE RNAFOLD -C. BUT FOR SOME F-ING REAOSN IT WORKS WORSE THAN WHEN IT IS LET BE.
% 	if ~isbalancedbr(str)
% 		[l r]=getunbalancedparsCAN(str);
% 		if ~isempty(r) & isempty(l)
% 			hairpins(ihairpins, 1) = hairpins(ihairpins, 1) + max(r);
% 			str = tg_br_aln_(hairpins(ihairpins, 1) : hairpins(ihairpins, 2));
% 			if ~isbalancedbr(str)
% 				[l r]=getunbalancedparsCAN(str);
% 				str([l r]) = '.';
% 			end
% 		end
% 	end
	
%     sum(str == '0' | str == '1' | str == '2' | str == '3' | str == '-') / sum(str == '.' | str == '(' | str == ')');
	if sum(str == '0' | str == '1' | str == '2' | str == '3' | str == '-') / sum(str == '.' | str == '(' | str == ')') <= HAIRPIN_INCONSISTENCE_TRESHOLD_RNAFOLD
%         disp(['The hairpin # ' num2str(ihairpins) ' too consistent for RNAfold prediction.'])
	else
        sq = aln_(2).Sequence(hairpins(ihairpins, 1) : hairpins(ihairpins, 2));
        sqi = find(sq ~= '-'); 
        sqgapsi = find(sq == '-');
        sq(sqgapsi) = [];
% 		disp(sq)
        str(find(str == '0' | str == '1' | str == '2' | str == '3')) = '.'; 
        str(sqgapsi) = [];
% 		disp(str)
        if ~USE_RNAFOLD_WITH_C
			id = fopen('2.fasta','w'); fwrite(id, sq); fclose(id); %eval(['!echo ' sq ' > 2.fasta'])
			
			if RNA_OR_UNA == 10
				eval(['!RNAfold < 2.fasta > 2.fasta.RNAfold'])
				[~, strRNAfold] = unix(['head 2.fasta.RNAfold']);
				strRNAfold = cut_string(strRNAfold, char(10));
				strRNAfold = strRNAfold{2}(1:length(strRNAfold{1}));
			else
% 				!UNAFold.pl --noisolate -X 1 2.fasta 2>&1 > NULL
% 				[s s]=unix('head -n1 2.fasta.ct');
% 				if s(1)=='0'
% 					strRNAfold = char(ones(1, length(sq)) .* '.');
% 				else
% 					!~/ViennaRNA-2.2.8/src/Utils/ct2b.pl 2.fasta_1.ct > 2.fasta_1.ct.br
% 					[~, strRNAfold] = unix('head 2.fasta_1.ct.br');
% 					strRNAfold = cut_string(strRNAfold, char(10));
% 					strRNAfold = strRNAfold{2}(1:length(strRNAfold{1}));
% 				end
				!hybrid-ss-min --noisolate --suffix DAT --mfold=5,1,1 2.fasta 2>&1 > NULL
				[s s]=unix('head -n1 2.fasta.ct');
				if s(1)=='0'
					strRNAfold = char(ones(1, length(sq)) .* '.');
				else
					!~/ViennaRNA-2.2.8/src/Utils/ct2b.pl 2.fasta.ct > 2.fasta.ct.br
					[~, strRNAfold] = unix('head 2.fasta.ct.br');
					strRNAfold = cut_string(strRNAfold, char(10));
					strRNAfold = strRNAfold{2}(1:length(strRNAfold{1}));
				end
			end
			
            if WORK_OUT_PREDICTED_MULTIPLE_LOOPS_HAIRPINS_AUTOMATICALLY
                strloops = find_loops(strRNAfold);
                if size(strloops, 1) > 1
                    if VERBOSE, disp(['paste.m: more than 1 loop in a RNAfold predicted hairpin structure for ' sqs_(isqs_).Header '. Trying RNAfold -C combined with copied base pairs instead.']), end
                    USE_RNAFOLD_WITH_C__ORIG = USE_RNAFOLD_WITH_C; COMBINE_COPIED_AND_RNAFOLD_BPS__ORIG = COMBINE_COPIED_AND_RNAFOLD_BPS;
                    USE_RNAFOLD_WITH_C = 1; COMBINE_COPIED_AND_RNAFOLD_BPS = 1;
                end
            end
		end
        if USE_RNAFOLD_WITH_C
			id = fopen('2.fasta', 'w'); fwrite(id, [sq char(10) str]); fclose(id); %	eval(['!echo -e ''' sq '\n' str ''' > 2.fasta'])
			eval(['!RNAfold -C < 2.fasta > 2.fasta.RNAfold'])
            [~, strRNAfold] = unix(['head 2.fasta.RNAfold']);
            strRNAfold = cut_string(strRNAfold, char(10));
			if isempty(strRNAfold{1})
				strRNAfold = [];
			else
				strRNAfold = strRNAfold{2}(1:length(strRNAfold{1}));
			end
        end
        if COMBINE_COPIED_AND_RNAFOLD_BPS
% THE STUPID RNAFOLD DOES NOT PREDICT SOME BASE PAIRS IN THE CONSTRAINT. PUT THEM BACK THEN...! 
            if ~isempty(strRNAfold)
				i=find(str ~= strRNAfold);
				ii=find(str~='.');
				strRNAfold(ii)=str(ii);
			end
        end
        if exist('USE_RNAFOLD_WITH_C__ORIG') & exist('COMBINE_COPIED_AND_RNAFOLD_BPS__ORIG')
            USE_RNAFOLD_WITH_C = USE_RNAFOLD_WITH_C__ORIG; COMBINE_COPIED_AND_RNAFOLD_BPS = COMBINE_COPIED_AND_RNAFOLD_BPS__ORIG;
			clear USE_RNAFOLD_WITH_C__ORIG COMBINE_COPIED_AND_RNAFOLD_BPS__ORIG
		end
		if isempty(strRNAfold)
% 			if HELP_HELICES_TO_FOLD
% 				
% 			else
				if VERBOSE, disp('paste.m: the RNAfold-predicted structure is empty. No paste.'), end
% 			end
		else
% CHECK IF THE PREDICTED EE STRUCTURE IS NOT EMPTY
			if all(strRNAfold == '.')
				 if VERBOSE, disp('paste.m: a RNAfold predicted hairpin structure is unstructured. No paste.'), end
			else
% CHECK IF THE PREDICTED EE STRUCTURE IS MATCHING
				if sum(strRNAfold == '(') ~= sum(strRNAfold == ')')
					 if VERBOSE, disp('paste.m: a RNAfold predicted hairpin structure is imbalanced.'), end
				else
% CHECK IF THE PREDICTED EE STRUCTURE HAS ONE LOOP ONLY, AS IT MUST STAY AN EE, I.E. AN ENCLOSED (STRUCTURE) ELEMENT 
					strloops = find_loops(strRNAfold);
					if size(strloops, 1) > 1
						 if VERBOSE, disp('paste.m: more than 1 loop in a RNAfold predicted hairpin structure.'), end
					else
						s = char(ones(1, hairpins(ihairpins, 2) - hairpins(ihairpins, 1) + 1) .* '-');
						s(sqi) = strRNAfold;
% CHECK IF OVERLAPPING HAIRPINS DO NOT REWRITE PARENTHESIS OF EACH OTHER
						if ihairpins > 1 & hairpins(ihairpins, 1) < hairpins(ihairpins - 1, 2) & (any(tg_br_aln_(hairpins(ihairpins, 1) : hairpins(ihairpins - 1, 2)) == '(') | any(tg_br_aln_(hairpins(ihairpins, 1) : hairpins(ihairpins - 1, 2)) == ')'))
% COPY ONLY THAT PART OF THIS HAIRPIN THAT DOES NOT REWRITE PARENTHESIS OF THE PREVIOUS HAIRPIN
							o = hairpins(ihairpins - 1, 2) - hairpins(ihairpins, 1) + 1;
							s = s(o + 1 : length(s));
							[lun run] = getunbalancedparsCAN(s);
							s([lun run]) = '.';
							tg_br_aln_(hairpins(ihairpins, 1) + o : hairpins(ihairpins, 2)) = s;
						else
							tg_br_aln_(hairpins(ihairpins, 1) : hairpins(ihairpins, 2)) = s;
						end
					end
				end
			end
		end
	end
end
% eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])
if PASTE_STEMS
% FIND_STEMS FINDS STEMS THAT MAY INCLUDE HAIRPINS IN ONE OF STRANDS OF STEMS. MOST LIKELY, THEY ARE THE REASON WHY BRS OF RANDOM 18S SEQUENCES ARE UNBALANCED. 
% MAYBE THE OCCASIONAL ERRORS FOR RANDOM 18S SEQUENCES ARE CAUSED BY THIS TOO. SO FAR I WATCHED THE ERRORS FOR RANDOM 18S RRNAS ONLY. INTERESTING IS THAT IT 
% MAKES NO ERRORS FOR REGULAR 18S SEQUENCES.
% 	stems=find_stems(tg_br_aln_, hairpins);
stems=find_stems1(tg_br_aln_, hairpins);
for istems = 1 : size(stems, 1)
    str = [tg_br_aln_(stems(istems, 1) : stems(istems, 2)) tg_br_aln_(stems(istems, 3) : stems(istems, 4))];
	
%     sum(str == '0' | str == '1' | str == '2' | str == '3' | str == '-') / sum(str == '.' | str == '(' | str == ')');
if sum(str == '0' | str == '1' | str == '2' | str == '3' | str == '-') / sum(str == '.' | str == '(' | str == ')') <= STEM_INCONSISTENCE_TRESHOLD_RNAFOLD
%         disp(['The stem # ' num2str(istems) ' too consistent for RNAfold prediction.'])
	else 
        sq = aln_(2).Sequence(stems(istems, 1) : stems(istems, 2));
        sqi = find(sq ~= '-'); 
        sq(find(sq == '-')) = [];

        sq1 = aln_(2).Sequence(stems(istems, 3) : stems(istems, 4));
        sq1i = find(sq1 ~= '-');
        sq1(find(sq1 == '-')) = [];

        id = fopen('2.fasta', 'w'); fwrite(id, [sq char(10) sq1]); fclose(id);  %eval(['!echo -e ''' sq '\n' sq1 ''' > 2.fasta'])
        eval(['!RNAduplex < 2.fasta > 2.fasta.RNAduplex'])
        [~, str] = unix(['head 2.fasta.RNAduplex']);
        s = strsplit(str, ' ');
        p = strsplit(s{2}, ',');
        p1 = strsplit(s{4}, ',');
        s=strsplit(s{1},'&');
    % FORMAT AND WRITE RNADUPLEX OUTPUT
    %     id=fopen([workdir_ '/2.fasta.RNAduplex.txt'],'w');
    %     fwrite(id, [sq char(10) ones(1, str2num(p{1}) - 1) .* ' ' s{1} char(10) sq1 char(10) ones(1, str2num(p1{1}) - 1) .* ' ' s{2}]);
    %     fclose(id);
    % CHECK IF THE PREDICTED EE STRUCTURE IS NOT EMPTY
        if all(s{1} == '.') | all(s{2} == '.')
             if VERBOSE, disp('paste.m: a RNAduplex predicted stem structure is unstructured. No paste.'), end
        else
    % CHECK IF THE PREDICTED EE STRUCTURE IS MATCHING
            if sum(s{1} == '(') ~= sum(s{2} == ')')
                if VERBOSE, disp('paste.m: a RNAduplex predicted stem structure is imbalanced.'), end
            else
                strR = char(ones(1, length(sq)) .* '.');
                strR(str2num(p{1}) : str2num(p{2})) = s{1};
                str = char(ones(1, stems(istems, 2) - stems(istems, 1) + 1) .* '-');
                str(sqi) = strR;
                tg_br_aln_(stems(istems, 1) : stems(istems, 2)) = str;

                strR = char(ones(1, length(sq1)) .* '.');
                strR(str2num(p1{1}) : str2num(p1{2})) = s{2};
                str = char(ones(1, stems(istems, 4) - stems(istems, 3) + 1) .* '-');
                str(sq1i) = strR;
                tg_br_aln_(stems(istems, 3) : stems(istems, 4)) = str;
            end
        end
	end
end
% eval(['!echo ''' tg_br_aln_ ''' | cat >> alm.br'])
end
clear i j i ki ks kis r ir li ri il ir d sq sqi s str loops iloops strloops sqgapsi strRNAfold ii stems
clear istems sq1 sq1i p p1 strR hairpins ihairpins lun run o 
