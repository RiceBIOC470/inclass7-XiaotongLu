%Inclass assignment 7. 
% 1. The gene Cdx2 is a crucial transcription factor involved in number of
% developmental stages. Use the UCSC genome browser to answer the following questions
% about it:

% A. What human chromosome is it located on?
%XiaotongLu:
It is located on Chro 13.
% B. How many exons does it have?
The total exon counts 3. Their position in coding region is hg38 chr13:27,963,115-27,969,006 and there total size is 5,892;
% C. What is the precise position of its stop codon in the genome?T
The stop codon is in the position 27,969,006.

% D. Identify at least one difference in sequence between human and mouse
% CDX2.
They have different sizes in coding region as the size of Cdx2 in human is 5892 but 5141 in mouse.

% E. In which human tissues is it expressed most abundantly?
The tissue is colon-transverse.

%2. A. Use the unigene database to find the accession number for a genbank
% entry containing the complete coding sequence of Cdx2. 

%XiaotongLu
NM_001265

% B. Use MATLAB to read the genbank information corresponding to that
% accession number.
gene=getgenbank('NM_001265')

% C. Use the information read in to find the position of the start and stop
% codon within the sequence. What are the parts of the sequence before the start codon 
% and after the stop codon?

%Xiaotong Lu
Before the start codon, there should be TATA box which is also a part for promoter. After the stop
codon, there will be polyA tail to stop the translation of DNA.

first solution:
gene=getgenbank('NM_001265');
geneinfo=gene.CDS;
Startnstoppostion=geneinfo.location


Startnstoppostion =

    '363..1304'

%XiaotongLu
second solution:
 J=gene.Sequence % J=fastaread('sequence.fasta')
 Dnaseq=upper(J);
startcodon=strfind(Dnaseq,'ATG');
stopcodon=[strfind(Dnaseq,'TAA'),strfind(Dnaseq,'TAG'),strfind(Dnaseq,'TGA')];
firststopcodon=zeros(1,length(startcodon));
for ii=1:length(startcodon)%do the circulation to every startcodon and start the next criculation until find the nearest stopcodon for the 1st startcodon. 
    ORFlen=stopcodon-startcodon(ii);%find all ORF and store in one array %use all stopcodon to minus every startcodon one by one
    good_length=1e8;
    good_index=0;
    for jj=1:length(ORFlen)
        if ORFlen(jj)>0&&mod(ORFlen(jj),3)==0&&ORFlen(jj)<good_length
            good_length=ORFlen(jj);%find all proper ORF
            good_index=jj;%number for every startcodon
        end
    end
    if good_index>0
        firststopcodon(ii)=stopcodon(good_index);%length(ORFlen)=length(stopcodon) %not place the ORF in order so the ORF will be a stopcodon corresponding to the startcodon
    else
        firststopcodon(ii)=startcodon(ii);
end
end
stopposition=firststopcodon
startposition=startcodon


% D. Use the protein_id to read the information on the protein. Use the
% information you read to determine where the homeobox domain of the protein is.
% Hint: see the field "Features". 
proteindat=getgenpept(gene.CDS.protein_id);
proteinInfo=proteindat.Features;
Region          190..242                                                  '
    '                /region_name="Homeobox"                                   ';
    '                /note="Homeobox domain; pfam00046"                        ';