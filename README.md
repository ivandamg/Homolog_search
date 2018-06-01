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


2. Get info
