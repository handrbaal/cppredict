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

sc__ = sc_;
sc_sq___ = sc_sq_;
ldiff___ = ldiff_;
sq_l_rat___ = sq_l_rat_;

RAND = 1; 
pactool_pairwise

if TEST_EACH
% PACTOOL_PAIRWISE.m WILL GENERATE EXACTLY # OF QUERIES * RANDSQS OF RANDOM STRS
	aval_ = ones(1, length(sc__)) .* 100000;
	j = 1;
	for i = 1 : RANDSQS : length(sc_)
		aval_(j) = round((mean(sc_(i : i + RANDSQS - 1)) - sc__(j)) ./ std(sc_(i : i + RANDSQS - 1)) .* 10) / 10;	
% 		[h aval_(j)] = ttest2(sc__(j), sc_(i : i + RANDSQS - 1));
		j = j + 1;
	end
end

if TEST_SOURCE
% TTEST2
% pval_ = ones(1, length(sc__)) .* realmax;
% for isc = 1 : length(sc__)
% 	[h pval_(isc)] = ttest2(sc__(isc), sc_);
% end

% MY TEST

aval_ = round((mean(sc_) - sc__) ./ std(sc_) .* 10) / 10;
% i = find(aval_ <= 1);
% j = find(aval_ > 1 & aval_ <= 2);

% PVALS
% 
% pd = fitdist(sc_','normal');
% Y = cdf(pd, 0 : max([sc_ sc__]));
% pval_ = zeros(1, length(sc__));
% j = find(sc__ > 0);
% pval_(j) = Y(sc__(j));
% % i = j(find(Y(sc__(j)) > .05));
end

sc_ = sc__;
sc_sq_ = sc_sq___;
ldiff_ = ldiff___;
sq_l_rat_ = sq_l_rat___;

clear Y i j pd sc__ sc_sq___ ldiff___ sq_l_rat___ h 
