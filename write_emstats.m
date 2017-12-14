function res=write_emstats(filein,root, imindex,Asum,g, axdiameter);

%write num axons, diameters and g ratio stats to a file
%one row per image, aquired at 7100 mag (neccessary for estimating g ratio)
%one file per white matter tract and genotype, including all images aquired
%sample inputs

% filein='/Users/alex/Desktop/andy/ematpathlab/FORNIX/NOS/140228_1_NOS_f/2650/140228_1\ NOS_001.TIF/';
% root='FORNIX/NOS';
%res=write_emstats(filein,root, Asum,mygratios, mydiameters); 

g=g(~isnan(g));
pathout='/Users/alex/Desktop/andy/ematpathlab/stats_em/';

label=regexp(root, '/', 'split');
label=char(strcat(label{1}, '_', label{2})); %FORNIX_NOS


[mystrings,pos1,pos2]=regexp(filein,root, 'split');
mystrings2=regexp(mystrings{2},'/', 'split');
myname=[mystrings2{2} num2str(imindex)];
[g1,numg]=hist(g,10);
[fi1,numfi]=hist(axdiameter,10);
numaxons=floor(Asum);

fileout=[pathout label '_stats.csv' ];
mystats=[ numaxons mean(axdiameter) std(axdiameter)  mean(g) std(g) g1 -1000 numg -1000 fi1 -2000 numfi -2000];
numreps=numel(mystats);

hdr=['%s ' repmat('%4.4f ', 1, numreps), '\n'];
    
fid=fopen(fileout,'a');
fprintf(fid, hdr,myname, mystats);
fclose(fid)

res=1
%dlmwrite(fileout,[ numaxons mean(g) std(g) g1 numg], 'delimiter', '\t', '-append', 'roffset', 1, 'coffset',1, 'precision', '%6.4f%');;
    