#/bin/bash

#Folders with list of DEGs
#./list_genes_FC2 -> lista completa y solo la columna de genes para cada comparación (q_val<0.05, FC>2)
	#./inf	-> lista completa de valores infinitos para cada comparación
			#List DEGs in which the log2FC were inf or -inf because 
			in some condiions the RPKM values were 0.
	#./up_regulated -> lista completa de up-regulated genes in each regime/condition without inf values
		#./inf ->  lista completa de valores inf para los genes up-regulated

#./list_genes -> lista copleta de genes para cada comparación (q_val<0.05)
	#./up_regulated -> lista completa de up-regulated genes in each regime/condition with inf values
