Notebooks:

gff_analyser -> Parser für GTF Files. object_list = gffBuilder.build_gff3_class(file=GTF_PATH) zum parsen der Datei. Die Liste object_list enthält ein Object zu jeder input File.
		Über die generate_*_gtf() Methoden können die entsprechenden GTF Files erzeugt werden. Als Input wird hier die erzeugte object_list benötigt.


get_gtf.ipynb -> Erstellen der GTF Files. Enhancer und Blacklisted Regions Files müssen noch seperat runtergeladen werden. Für Promotoren, TSS und gene bodies muss zunachst eine seperate File.
                 über die Methode .generate_feature_gtf() erzeugt werden. Die erstellten Files werden in einem "/out" Ordner abgelegt. 

calc_overlap_pct.py -> Angepasste Wrapperfunktionen der sctoolbox. pct_fragments_in_promoters() zum verrechnen einzelner Feature GTF Files und pct_fragments_in_features() zum verrechnen aller
		       GTF-Files in "/out" Ordner.

scATAC_frame.ipynb -> gff_analyser angebunden an Yousef's Framework zur scATAC analyse.
