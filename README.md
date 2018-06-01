# Homolog_search

1. Get target genomes from NCBI.

- Download assemblies from https://www.ncbi.nlm.nih.gov/assembly/?term= . All files prot,nucl,gbf,gbk,gff, etc...

    - extract info from assemblies. Strain, Species, Completeness, size, etc.
        
    - Extract all compressed files in genome folder
        
            gunzip */GCF_*.faa.gz
            gunzip */GCF_*[0-9]_genomic.fna.gz
    - Copy protein file and assembly file on front folder
           
            cp */*.faa .
            cp */*[0-9]_genomic.fna .
            
    - Get strains info and change name of assembly (.fna) or (.faa) files to Gender_Species_Strains_Source.faa

            for i in $(ls */*assembly_stats.txt); do echo $i  | cat $i | grep "Assembly level"  | cut -d ':' -f2 | sed 's/Genome//g' ; done

            for i in $(ls */*assembly_stats.txt); do echo $i  | cat $i | grep "total-length"  | grep -v 'Primary' | grep -v 'sequence' | cut -d'-' -f2 | sed 's/length//'; done

            for i in $(ls */*assembly_stats.txt); do echo $i  | cat $i | grep "total-gap-length" | grep -v 'of' | cut -d'-' -f3 | sed 's/length//' ; done
            for i in $(ls */*assembly_report.txt); do echo $i  | cat $i | grep "strain=" | cut -d':' -f2 | sed 's/strain=//' | sed 's/\//_/g' | sed 's/,/_/g' | tr '\(' ' ' | tr '\)' ' ' | sed 's/ //g' ; done
            
            for i in $(ls */*assembly_report.txt); do echo $i  | cat $i | grep "Assembly name"  | cut -d':' -f2; done

    - Make summary info strains in excel or text editor.

- Alternative: Generation new annotation with Prokka.

        for i in $(ls *.fna); do echo $i;  ~/software/prokka-1.12/prokka/bin/prokka --outdir Annotation_$(echo $i | cut -d'_' -f3) --genus $(echo $i | cut -d'_' -f1) --species $(echo $i | cut -d'_' -f2) --strain $(echo $i | cut -d'_' -f3) --locustag Ab_$(echo $i | cut -d'_' -f3) --prefix $(echo $i | cut -d'_' -f1,2,3)_Prokka --rfam --usegenus $i  ;  done


2. Blast query sequence to target strains.

- Make protein databases of genomes

        a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1,2,3) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_prot_genomes/$(echo $i | cut -d'_' -f1,2,3)_db ; done

- Blast proteins

        a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1,2,3) ;blastp -db db_prot_genomes/$(echo $i | cut -d'_' -f1,2,3)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_prot_genomes/blastProt_DNA-uptake_$(echo $i | cut -d'_' -f1,2,3).xml -query ~/Documents/EPFL/Competence_in_NosoPath/Reference_Competence_Proteins/List_competence_Acinetobacter_baylyi_prot_V2.faa ; done

- Extract sequences of proteins only complete ones

        for gene in $(echo ACIAD0360_pilD ACIAD0361_pilC ACIAD0362_pilB ACIAD0558_pilF ACIAD0695_fimT ACIAD0911_pilU ACIAD0912_pilT ACIAD3314_comF ACIAD3315_comE ACIAD3316_comC ACIAD3318_comB ACIAD3319_pilV ACIAD3321_fimU ACIAD3338_comP ACIAD3355_comQ ACIAD3356_comL ACIAD3357_comO ACIAD3359_comN ACIAD3360_pilM ACIAD2639_comA ACIAD3064_comEA ACIAD0209_dprA ACIAD1385_recA ACIAD3449_ssb ACIAD0242_comM) ; do echo $gene;
        for i in $(ls blastProt*.xml); do echo $i ; star="$(cat $i | cut -f2 |sed 's/^/blastdbcmd /g' | sed 's/ / -entry /g' | awk '$3="\x27"$3"\x27"')"; 
        # reorder coordinates; range="$(cat $i | cut -f9,10 | while read line; do echo $line | sed 's/ /\n/g' | sort | gawk '{line=line " " $0} END {print line}' ; done | sed 's/^ /-range/g' | sed 's/ /-/' | sed 's/range/range /')" ; 
        # extract name and add db part; # nb of lines in xml; lines="$(cat $i | wc -l)"; 
        export lines # put variable in open environment ; db="$(echo $i | cut -d'_' -f3,4,5,6,7 | sed 's/\.xml/_db/g' | sed 's/^/-db /' | perl -ne 'print $_ x $ENV{lines}')"; 
        # to save in file;
        nam1="$(cat $i | cut -f1  | cut -f1 -d'/')";
        nam2="$(echo $i | cut -d'_' -f3,4,5,6,7 | sed 's/\.xml//g' | perl -ne 'print $_ x $ENV{lines}')";end="$(paste <(echo "$nam1") <(echo "$nam2") --delimiters '_' | sed 's/^/> Seq_/'| sed 's/$/.fa/' )";
        # aseemble everything ;paste <(echo "$star") <(echo "$db") <(echo "$range") <(echo "$end") --delimiters ' ' | grep $gene | head -1; done ; done

- Change name of sequences

        for i in $(ls Seq*.fa); do cat $i | sed "s/>.*/>${i} /" | sed 's/\.fa//'| sed 's/,//' | sed 's/whole genome shotgun sequence//' > $(echo $i)2; done


3. Alternative to point 2. 
HMM search for finding homologs

- Assemble in one file protein sequence in different strains 

        cat Seq_comA_Streptococcus_pneumoniae_*faa2 > ALL_comA.faa

- Align
Mafft https://mafft.cbrc.jp/alignment/server/index.html
clustalW 
clustalO

- Create HMM profile of alignment

        for i in $(ls ALL*.aln); do echo $i; hmmbuild $(echo $i | sed 's/.aln/.hmm/g') $i ; done

- Find homolog on other strain of other species using hmm profile
 
    -on a single gene
 
        hmmsearch ~/Documents/EPFL/Competence_in_NosoPath/GramNEG_competencegenes/Acinetobacter_baumannii/db_prot_genomes/ALL_comM.hmm  ~/Documents/EPFL/Competence_in_NosoPath/GramNEG_competencegenes/Campylobacter_sps/Campylobacter_jejuni_ATCC35925_ASM202830v1_Complete.faa > ~/Documents/EPFL/Competence_in_NosoPath/GramNEG_competencegenes/Campylobacter_sps/Competence_genes_hmm/Campylobacter_jejuni_ATCC35925.out

    - on multiple genes
    
           for i in $(ls *.hmm); do echo $i ; for strain in $(ls ~/Documents/EPFL/Competence_in_NosoPath/GramNEG_competencegenes/Campylobacter_sps/*.faa); do echo $strain ; hmmsearch $i $strain > $(echo $strain | cut -d'.' -f1).out ; done; done

- Look results and Select top hit

        for i in $(ls ACIAD*comM*); do echo $i ; head -30 $i ; done

- Select top hit and extract protein in new file (i.e. Gene_Genus_species_strain.fa)


4. Compare results of HMM profile to blast search

blast perform well when gene is conserved within strains. 
case of looking homolog in close related species or within the same species

HMM profile perform better when query is not well conserved. 
case of looking homolog in far related species. Need to get sequences of the genes that are not conserved to create the HMM profile.



