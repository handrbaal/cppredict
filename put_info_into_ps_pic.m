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

%put info_ into picture and title of *.ps files
%cd ~/rna/temp

id=fopen(psname_, 'r'); t = char(fread(id))'; fclose(id);
id=fopen(psname_, 'w');

i=strfind(t, '/sequence (\') + length('/sequence (\') + 1;
j=strfind(t(i : length(t)), '\')+i-2;
seqlen = length(t(i:j));
if seqlen < 50
	fntsz=7; 
else
		if seqlen < 100
			fntsz=15;
		else
			if seqlen < 150
				fntsz = 25;
			else
				if seqlen < 200
					fntsz = 30;
				else
					fntsz = 40;
				end
			end
		end
	end

i=strfind(t, '% show it');
if ~isempty(i)
	t(i : length(t)) = []; 
	t = [t ['% show it' char(10) '/Times-Roman findfont' char(10) num2str(fntsz) ' scalefont' char(10) 'setfont' char(10) 'newpath' char(10) 'xmin ymax moveto' char(10) '(' info_ ') show' char(10) 'showpage' char(10) 'end' char(10) '%%EOF' char(10)]];
end
fwrite(id, t);
fclose(id);

%cd ~/rna
clear id t i j fntsz 
