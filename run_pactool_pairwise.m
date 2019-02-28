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

% disp([char(10) 'sr__fn_br_=' sr__fn_br_ ', fn__=' fn__ char(10) 'BOOTSTRAP=' num2str(BOOTSTRAP) ', RAND=' num2str(RAND) char(10) 'THROW_OUT_COPIED_NC_BPS=' num2str(THROW_OUT_COPIED_NC_BPS) ', PASTE=' num2str(PASTE) char(10) 'ALM=' ALM char(10) 'PASTE_STEMS=' num2str(PASTE_STEMS) ', WORK_OUT_PREDICTED_MULTIPLE_LOOPS_HAIRPINS_AUTOMATICALLY=' num2str(WORK_OUT_PREDICTED_MULTIPLE_LOOPS_HAIRPINS_AUTOMATICALLY) ', USE_RNAFOLD_WITH_C=' num2str(USE_RNAFOLD_WITH_C) ', COMBINE_COPIED_AND_RNAFOLD_BPS=' num2str(COMBINE_COPIED_AND_RNAFOLD_BPS) ', EXTEND_MECHANICALLY_LONELY_PAIRS=' num2str(EXTEND_MECHANICALLY_LONELY_PAIRS) char(10) 'SHOW_PREDICTED_STRUCTURE_INDIVIDUALLY=' num2str(CLOSE_PREDICTED_STRUCTURES_INDIVIDUALLY) ', CLOSE_PREDICTED_STRUCTURES_INDIVIDUALLY=' num2str(CLOSE_PREDICTED_STRUCTURES_INDIVIDUALLY) ', PAUSE_BETWEEN_INDIVIDUAL_STRUCTURES=' num2str(PAUSE_BETWEEN_INDIVIDUAL_STRUCTURES) ', MAKE_FIGURE_HEADERS=' num2str(MAKE_FIGURE_HEADERS) ', MAKE_FIGURES=' num2str(MAKE_FIGURES) ', MERGE_FIGURES_INTO_PDF=' num2str(MERGE_FIGURES_INTO_PDF) ', SHOW_FIGURES_MERGED=' num2str(SHOW_FIGURES_MERGED) ', SHOW_TEST_FIGURE=' num2str(SHOW_TEST_FIGURE) ', VERBOSE=' num2str(VERBOSE) '.'])
% if ~RAND, s = input('Are the pars OK? [y]/n', 's'); if s == 'n', return, end, end

if ~exist(fn__) | ~exist(sr__fn_br_), disp([fn__ ' and/or ' sr__fn_br_ ' do not exist.']), return, end

id=fopen(sr__fn_br_); f=fread(id); fclose(id); if f(length(f)) ~= 10, disp([sr__fn_br_ ' does not end by CF (char(10)). Stop.']), pause, end

if 0
	fn_ = fn__;
	clean_headers
	fn_=[fn_ '.cleaned']
	clean_fasta_RNA
	fn__ = [fn_ '.cleaned']
else
	if ~check_fasta(fn__), disp(['run_pactool_pairwise.m: Stop ''cos check_fasta.m in ' fn__ '.']), return, end
	if ~check_br(sr__fn_br_), disp(['run_pactool_pairwise.m: Stop ''cos check_br.m in ' sr__fn_br_ '.']), return, end
end

sr__fn_bpseq_ = [sr__fn_br_(1 : strfind(sr__fn_br_, '.br') - 1) '.bpseq'];
br2bpseq(sr__fn_br_); 

if BOOTSTRAP == 1
	RAND = 0; SHOW_PREDICTED_STRUCTURES_INDIVIDUALLY = 0; 
else
	SHOW_TEST_FIGURE = 0;
end

pactool_pairwise

if BOOTSTRAP == 1
	testalg='test_JJ_live';
% 	testalg='test_JJ_live_FE';
% 	testalg='test_live';
% disp(testalg)
	eval(testalg);
	eval(['!rm -f RAND_*.br'])
else
	testalg='';
end

% MAKE THE OUTPUT LIST

s = ['[sequence name: sequence similarity to template (%), seq. length query/template ratio, seq. length query/template difference (# of nucleotides), RNAdistance score to template[, z-score (if BOOTSTRAP=1)]]'];
disp([char(10) char(10) 'These been done: ' char(10) s char(10)])
for iind = 1 : length(ind_)
	isqs = ind_(iind);
% 	s = [sqs_original_(isqs).Header ': ' char(9) num2str(sc_sq_(isqs)) '%, ' num2str(round(sq_l_rat_(isqs) * 100) / 100 ) ' (' num2str(ldiff_(isqs)) ' nts), ' num2str(sc_(isqs))];
	s = [sqs_original_(isqs).Header char(9) num2str(sc_sq_(isqs)) char(9) num2str(round(sq_l_rat_(isqs) * 100) / 100 ) char(9) num2str(ldiff_(isqs)) char(9) num2str(sc_(isqs))];	
	if exist('aval_')
% 		s = [s ', ' num2str(aval_(isqs))];
		s = [s char(9) num2str(aval_(isqs))];
	end
% 	s =  [s '.'];
	disp(s)
end

% if BOOTSTRAP & SHOW_TEST_FIGURE
% 	labels = {sqs_original_.Header};
% figure('Visible','off'); hold on, grid, box
% 	ylabel(['Structure similarity to reference structure (''-''), ' char(10) 'biological relevance (''*'', ''+'', ''x''),' char(10) 'Sq. l. ratio with reference sq. (''--'').'])
% 	set(gca,'TickLength',[0 0],'xtick',1:1:length(labels),'xticklabel',labels,'TickLabelInterpreter','none','XTickLabelRotation',90,'FontSize',7)
% 	if exist('aval_')
% 		i=find(aval_>=2); plot(i,aval_(i),'k*'), i = find(aval_<2 & aval_>=1); plot(i, aval_(i),'k+'), i = find(aval_<1); plot(i, aval_(i),'kx')
% 	end
% 	plot((max(sc_)-sc_)./10,'k')
% 	plot(sq_l_rat_ .* 100, 'k--')
% 	% 	i = find(pval_<=.05); plot(i, pval_(i),'r*'), i = find(pval_>.05); plot(i, pval_(i),'rx')
% 	a = axis; axis([a(1) length(labels) a(3:4)])
% % 	ss = cut_string(sr__h_, '_'); s = ss{1}; for iss = 2 : length(ss), s = [s '\_' ss{iss}]; end, title([s ', ' ALM ', ' testalg]) % ', ''-'' str. sim. ''x'' (non-relevant)\in<-\infty, \sigma), ''+'' (marginally relevant)\in<\sigma, 2\sigma), ''*'' (relevant)\in<2\sigma, \infty>.'],'interpreter','tex','FontSize',10)
% 	title([sr__h_ ', ' ALM ', ' testalg], 'interpreter','none') 
% 	saveas(gcf, [workdir_ '/' fn__ '.eps'], 'eps');
% 	eval(['!evince ' workdir_ '/' fn__ '.eps &'])
% end

if MAKE_FIGURES
pslist = []; 
for iind = 1 : length(ind_)
		isqs = ind_(iind);
		eval(['!RNAplot -t0 < ' '' sqs_original_(isqs).Header '.br' ''])
% RECONSTRUCT A .PS FILE NAME IF RNAPLOT 1) MADE PS FILE, 2) CUT A HEADER FROM A BR FILE IN THE PS FILE NAME
		if length(sqs_original_(isqs).Header) >= 43
			psname_ = [sqs_original_(isqs).Header(1 : 42) '_ss.ps'];
		else
			psname_ = [sqs_original_(isqs).Header '_ss.ps'];
		end
% SOMETIMES RNAPLOT DOES NOT MAKE A .PS FILE FOR UNKNOWN REASON
		if exist(psname_)
			eval(['!mv ' psname_ ' ' sqs_original_(isqs).Header '.ps']);
			psname_ = [sqs_original_(isqs).Header '.ps'];
% 			info_ = [sqs_original_(isqs).Header ', ' num2str(sc_(isqs))];
			info_ = sqs_original_(isqs).Header;
			if exist('aval_')
				info_ = [info_ ', ' num2str(aval_(isqs))];
			end
			info_ = [info_ ', ' num2str(sc_(isqs))];
			put_info_into_ps_pic
			pslist = [pslist psname_ ' '];
		end
	end
	if MERGE_FIGURES_INTO_PDF
		eval(['!cat ' pslist ' > merged.ps'])
		eval(['!ps2pdf merged.ps ' sr__h_ '__' fn__ '__' ALM '__' testalg '__pactool_pairwise.m.pdf'])
		!rm merged.ps
		if SHOW_FIGURES_MERGED
			eval(['!evince ' sr__h_ '__' fn__ '__' ALM '__' testalg '__pactool_pairwise.m.pdf &'])
		end
	end
end

% THIS IS FOR WHEN IT IS CALLED BY subopt_run_pactool_pairwise.m
if exist('it___')
	brlist = [];
	for iind = 1 : length(ind_)
		brlist = [brlist sqs_original_(ind_(iind)).Header '.br '];
	end
	eval(['!cat ' brlist ' > all_' num2str(floor((it___ - 1) / 3) + 1) '.br'])
	if MAKE_FIGURES	
		eval(['!mv ' sr__h_ '__' fn__ '__' ALM '__' testalg '__pactool_pairwise.m.pdf all_' num2str(floor((it___ - 1) / 3) + 1) '.pdf'])
	end
	eval(['!rm ' brlist])
end

if MAKE_FIGURES
	eval(['!rm ' pslist])
end

clear isqs id di t fn_ aln_ i s ss str sq labels d info_ psname_ ss iss s a name a sqs_original_ sc_rnd_ idall is ldiff_ 
clear lun nci run sc_sq_ ind_ iind pars_ sq_l_rat_ sr__q_ sr__h_ sr__br_ sr__bpseq_ avals head spsmax sps pslist brlist 
