function res=write_emstats(filein,root, Asum,g);
filein='/Users/alex/Desktop/andy/ematpathlab/FORNIX/NOS/140228_1_NOS_f/2650/140228_1\ NOS_001.TIF/';

root='FORNIX/NOS';
label=regexp(root, '/', 'split');
label=char(label{1});%'FORNIX';


[mystrings,pos1,pos2]=regexp(filein,root, 'split');
mystrings2=regexp(mystrings{2},'/', 'split');
myname=mystrings2{2}
pathout='/Users/alex/Desktop/andy/ematpathlab/stats_em/';
[g1,numg]=hist(g,10);
numaxons=floor(Asum);

fileout=[pathout label '_stats.csv' ];
mystats=[ numaxons mean(g) std(g) g1 numg];
numreps=numel([ numaxons mean(g) std(g) g1 numg]);
hdr=['%s ' repmat('%4.4f ', 1,numreps) '\n'];
    
fid=fopen(fileout,'a');
fprintf(fid, hdr,myname, mystats);
fclose(fid)
%dlmwrite(fileout,[ numaxons mean(g) std(g) g1 numg], 'delimiter', '\t', '-append', 'roffset', 1, 'coffset',1, 'precision', '%6.4f%');;
    