# Mebarcoding-en-peces-con-Qiime2-para-el-locus-COX1

_Martin Holguin Osorio_\
_junio de 2021_\
_Version 3_ 

## Creacion de archivo de metadatos y preparacion de espacio de trabajo
El archivo "metadata.tsv" tiene que seguir un orden  y una organizacion especifica para funcionar corectamente, mas informacion en [tutorial metadata en Qiime2](https://docs.qiime2.org/2020.11/tutorials/metadata/) , se recomienda usar el add-on keemei de google sheets para crear el archivo de "metadata.tsv" de la manera mas rapida y sencilla, mas [info aqui](https://keemei.qiime2.org).

* Tras descargar metadata.tsv usando google sheets:

```
#abro qiime2 en una terminal
conda activate qiime2-2021.8

#navego y ubico el archivo de metadata en la direccion de la carpeta donde voy a trabajar(en mi caso /home/martin/eDNA/1 )
cd /home/martin/eDNA/6/COI

#tomo las lecturas de los datos crudos y los pongo dentro de esta nueva carpeta "datos"
mkdir datos

#creo una carpeta donde se ubicaran las salidas de cada proceso (resultados)
mkdir salidas
```

## Importacion de datos
Luego de ver la naturaleza de estos datos (demultiplexados, con secuenciacion pareada y con los barcodes en las secuencias), defino que tengo que usar otro comando para importar los datos, asi que, dejo los nombres de los datos crudos (dejandolos como estan, sin alterar nada) y los ubico en la carpeta "datos"

```
#importo los datos a qiime2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path datos/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path salidas/secs_multiplexadas.qza
```

## Demultiplexacion
Separo la informacion de cada muestra, contenida en los archivos fastq, usando los barcodes de referencia que hay para cada una, en el archivo "metadata.tsv" 

```
#demultiplexo
qiime cutadapt demux-paired \
  --i-seqs salidas/multiplexed-seqs.qza \
  --m-forward-barcodes-file metadata.tsv \
  --m-forward-barcodes-column BarcodeSequence \
  --p-error-rate 0 \
  --o-per-sample-sequences salidas/secs_demultiplexadas.qza \
  --o-untrimmed-sequences salidas/untrimmed.qza \
  --verbose
  
#creo resumen de demultiplezacion y  creo un visualizador de este
qiime demux summarize \
  --i-data salidas/secs_demultiplexadas.qza \
  --o-visualization salidas/secs_demultiplexadas.qzv

#abro visualizador en el navegador
qiime tools view salidas/secs_demultiplexadas.qzv
#en el archivo "untrimmed.qza" quedan todas las secuencias sin asginar
```

## Eliminación de ruido (denoising con DADA2 )
Este Plugin corrige múltiples artefactos producto de la amplificación y secuenciación, une las lecturas pareadas y corta los extremos de las lecturas con baja calidad, para al final interpretar y agrupar la diversidad de secuencias como ASVs (variantes de secuencias de amplicones) en lugar de OTUs (unidades taxonómicas operativas). La eliminación de ruido detecta secuencias erróneas y las fusiona con la secuencia "madre" correcta, interpretándola como ASV, mientras que su alternativa (el agrupamiento ó clustering) intenta combinar un conjunto de secuencias (sin importar si contienen o no errores) en entidades biológicas significativas, las cuales se interpretan como OTUs. Ambos métodos de tratamiento de las lecturas (secuencias) buscan asignar secuencias al nivel taxonómico de especie y han sido presentados como alternativas mutuamente excluyentes (Antich et al., 2021)

```
#hago denoising escojiendo las posiciones de trimero a ojo
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs salidas/secs_demultiplexadas.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 3 \
  --p-trunc-len-f 222 \
  --p-trunc-len-r 201 \
  --o-table salidas/tabla.qza \
  --o-representative-sequences salidas/secs_representativas.qza \
  --o-denoising-stats salidas/resumen_denoising.qza
  
#genero visualizacion del archivo "table.qzv", el cual brinda informacion cuantitativa de las secuencias significativas
qiime feature-table summarize \
--i-table salidas/tabla.qza \
--o-visualization salidas/tabla.qzv \
--m-sample-metadata-file metadata.tsv

#resumen de denoising
qiime metadata tabulate \
--m-input-file salidas/resumen_denoising.qza \
--o-visualization salidas/resumen_denoising.qzv

#secuencias representativas
qiime feature-table tabulate-seqs \
--i-data salidas/secs_representativas.qza \
--o-visualization salidas/secs_representativas.qzv

#visualizo
qiime tools view salidas/tabla.qzv
```

## asignacion taxonomica 

#etapa experimental; entre las lineas esta la zona de experimentos, se hizo de esta manera ya que se tenia que escoger la mejor ruta de asignacion para mejorar y corregir varios detalles de la version anterior del pipeline
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************
### asignacion taxonomica con sklearn
qiime feature-classifier classify-sklearn \
--i-classifier ../BdD/Locales/COI/classifier_SuppMaterial2_FISHF2.qza \
--i-reads salidas/secs_representativas.qza \
--p-confidence 0.90 \
--p-read-orientation same \
--o-classification salidas/taxonomy_SuppMaterial2.qza
#Plugin error from feature-classifier:
#
#  pop from empty list
#
#Debug info has been saved to /tmp/qiime2-q2cli-err-65qmtwhk.log



#***********************************************************************************************************************
### asignacion taxonomica con vsearch
qiime feature-classifier classify-consensus-vsearch \
  --i-query salidas/secs_representativas.qza \
  --i-reference-reads ../BdD/Locales/COI/SuppMaterial2.qza \
  --i-reference-taxonomy ../BdD/Locales/COI/taxo_SuppMaterial2.qza \
  --p-perc-identity 0.90 \
  --p-min-consensus 0.51 \
  --p-maxaccepts 10 \
  --p-threads 12 \
  --o-classification salidas/vsearch-taxonomy_SuppMaterial2.qza

#creo visualizador para taxonomia
qiime metadata tabulate \
--m-input-file salidas/vsearch-taxonomy_SuppMaterial2.qza \
--o-visualization salidas/vsearch-taxonomy_SuppMaterial2.qzv

#creo un barplot de taxonomia
qiime taxa barplot \
--i-table salidas/tabla.qza \
--i-taxonomy salidas/vsearch-taxonomy_SuppMaterial2.qza \
--m-metadata-file metadata.tsv \
--o-visualization salidas/vsearch_taxa_bar_plots.qzv
#visualizo
qiime tools view salidas/vsearch_taxa_bar_plots.qzv



#***********************************************************************************************************************
### asignacion taxonomica hibrida con vsearch y sklearn
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query salidas/secs_representativas4.qza \
  --i-reference-reads ../BdD/Locales/COI/SuppMaterial2.qza \
  --i-reference-taxonomy ../BdD/Locales/COI/taxo_SuppMaterial2.qza \
  --i-classifier ../BdD/Locales/COI/classifier_SuppMaterial2.qza \
  --o-classification salidas/0.97-hybrid-vsearch-sklearn-taxonomy_SuppMaterial2.qza
#Plugin error from feature-classifier:
#
#  Command '['vsearch', '--fastx_subsample', '/tmp/qiime2-archive-w29c6cac/320ba9cd-ec4d-4aeb-a771-6050de6d7490/data/dna-sequences.fasta', '--sample_size', '1000', '--randseed', '0', '--fastaout', '/tmp/tmpljeec60b']' returned non-zero exit status 1.
#
#Debug info has been saved to /tmp/qiime2-q2cli-err-cbadm2i7.log




#***********************************************************************************************************************
### asignacion taxonomica con blast
qiime feature-classifier classify-consensus-blast \
--i-query salidas/secs_representativas.qza \
--i-reference-reads ../BdD/Locales/COI/SuppMaterial2.qza \
--i-reference-taxonomy ../BdD/Locales/COI/taxo_SuppMaterial2.qza \
--p-perc-identity 0.90 \
--o-classification salidas/blast-classified-table_SuppMaterial2.qza 

#creo visualizador de las secuencias agrupadas
qiime metadata tabulate \
  --m-input-file salidas/blast-classified-table_SuppMaterial2.qza \
  --o-visualization salidas/blast-classified-table_SuppMaterial2.qzv

#creo un barplot de taxonomia
qiime taxa barplot \
--i-table salidas/tabla.qza \
--i-taxonomy salidas/blast-classified-table_SuppMaterial2.qza  \
--m-metadata-file metadata.tsv \
--o-visualization salidas/blast-taxa_bar_plots_SuppMaterial2.qzv
#visualizo
qiime tools view salidas/blast-taxa_bar_plots_SuppMaterial2.qzv
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Escojo vsearch al 90% ya que este marcador (COX1) no suele relajarse tanto en su porcentaje de identidad segun he visto en la literatura, sin embargo si se relaja al 70 muestra mayor numero de
#spp y una mayor proporcion de eDNA asignada en cada muestra, PERO, detecta varias especies ajenas al area de estudio
#Blast al 90% tambien da resultados similares, en cuanto a las asignaciones



#########################################################################################################################################################################################################################
###################################################### Filogenia ambiental ###################################################

#hago filogenia con la base de datos
#el comando phylogeny align-to-tree-mafft-fasttree hace todo el pipeline para generar la filogenia
    #qiime alignment mafft ...
    #qiime alignment mask ...
    #qiime phylogeny fasttree ...
    #qiime phylogeny midpoint-root ...
#tecnicamente es un pipeline dentro de este pipeline (hace lo mismo en menos lineas, dentro de este pipeline)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences salidas/secs_representativas.qza \
  --output-dir salidas/filogenia_ambiental

#creo visualizador del arbol filogenetico con las spp detectadas en la asignacion taxonomica con vsearch al 90%
qiime empress community-plot \
    --i-tree salidas/filogenia_ambiental/rooted_tree.qza \
    --i-feature-table salidas/tabla.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file salidas/vsearch-taxonomy_SuppMaterial2.qza \
    --o-visualization salidas/empress-tree.qzv
    
    
#######################################################################################################################################################################################################################################
################################### filtracion #################################
#como hay rasgos(ASVs) de otras cosas ajenas a peces filtro y dejo solo la informacion de peces

#filtro la tabla 
qiime taxa filter-table \
  --i-table salidas/tabla.qza \
  --i-taxonomy salidas/vsearch-taxonomy_SuppMaterial2.qza  \
  --p-include Actinopterygii \
  --o-filtered-table salidas/tabla_filtrada.qza
#creo visualizador de la tabla
qiime feature-table summarize \
--i-table salidas/tabla_filtrada.qza \
--o-visualization salidas/tabla_filtrada.qzv \
--m-sample-metadata-file metadata.tsv
#visualizo
qiime tools view salidas/tabla_filtrada.qzv  


#filtro las secuencias
qiime taxa filter-seqs \
  --i-sequences salidas/secs_representativas.qza \
  --i-taxonomy salidas/vsearch-taxonomy_SuppMaterial2.qza  \
  --p-include Actinopterygii \
  --o-filtered-sequences salidas/secs_representativas_filtradas.qza \
  

#**************************************************************************************************************************
### vuelvo y hago asignacion taxonomica con vsearch para verificar que solo quedan secuencias de peces
qiime feature-classifier classify-consensus-vsearch \
  --i-query salidas/secs_representativas_filtradas.qza \
  --i-reference-reads ../BdD/Locales/COI/SuppMaterial2.qza \
  --i-reference-taxonomy ../BdD/Locales/COI/taxo_SuppMaterial2.qza \
  --p-perc-identity 0.90 \
  --p-min-consensus 0.51 \
  --p-maxaccepts 10 \
  --p-threads 12 \
  --o-classification salidas/vsearch-taxonomy_SuppMaterial2_filtrada.qza

#creo visualizador para taxonomia
qiime metadata tabulate \
--m-input-file salidas/vsearch-taxonomy_SuppMaterial2_filtrada.qza \
--o-visualization salidas/vsearch-taxonomy_SuppMaterial2_filtrada.qzv

#creo un barplot de taxonomia para verificar visualmente que solo quedan secuencias de peces
qiime taxa barplot \
--i-table salidas/tabla_filtrada.qza \
--i-taxonomy salidas/vsearch-taxonomy_SuppMaterial2_filtrada.qza \
--m-metadata-file metadata.tsv \
--o-visualization salidas/vsearch_taxa_bar_plots_filtrada.qzv
#visualizo
qiime tools view salidas/vsearch_taxa_bar_plots_filtrada.qzv


#******************************************************************************************
### vuelvo y hago la Filogenia ambiental ###################################################

#hago filogenia con la base de datos
#el comando phylogeny align-to-tree-mafft-fasttree hace todo el pipeline para generar la filogenia
    #qiime alignment mafft ...
    #qiime alignment mask ...
    #qiime phylogeny fasttree ...
    #qiime phylogeny midpoint-root ...
#tecnicamente es un pipeline dentro de este pipeline (hace lo mismo en menos lineas, dentro de este pipeline)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences salidas/secs_representativas_filtradas.qza \
  --output-dir salidas/filogenia_ambiental_filtrada

#creo visualizador del arbol filogenetico con las spp detectadas en la asignacion taxonomica con vsearch al 90%
qiime empress community-plot \
    --i-tree salidas/filogenia_ambiental_filtrada/rooted_tree.qza \
    --i-feature-table salidas/tabla_filtrada.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file salidas/vsearch-taxonomy_SuppMaterial2_filtrada.qza \
    --o-visualization salidas/empress-tree_filtrada.qzv


#######################################################################################################################################################################################################################################
################################### analisis de rarefacción y diversidades alfa_beta #################################
#la conversión de recuentos a proporciones y la rarefacción de los recuentos son dos formas de tener en cuenta y corregir la profundidad de muestreo desigual de las muestras. Aunque las
#proporciones son más intuitivas y fáciles de entender, la rarefacción es mas popular a la hora de calcular los índices de diversidad, especialmente para algunos índices de 
#diversidad que tienen en cuenta el número de especies utilizando la presencia/ausencia en lugar de la abundancia, como los índices de Sørensen y Jaccard. 
#Cuando los recuentos de las lecturas se enrarecen, algunos recuentos bajos llegan a cero, reduciendo efectivamente el número de especies corrigiendo la profundidad de muestreo desigual en estos índices. 

#**************************************************************************************************************************
####IMPORTANTE: esto solo se hace en este primer ensayo piloto
### En este primer ensayo piloto es necesario hacer una rerefaccion con todas las ASVs para definir cual fue el kit de extraccion que mas diversiadd (ASVs rescato)
#con el plugin diversity core-metrics, hago la rarefaccion y ademas se calculan todas las medidas de diversidad estandar
qiime diversity core-metrics-phylogenetic \
--i-phylogeny salidas/filogenia_ambiental/rooted_tree.qza \
--i-table salidas/tabla.qza \
--p-sampling-depth 10000 \
--m-metadata-file metadata.tsv \
--output-dir salidas/core-metrics-results_KitsE
##El parámetro "--p-sampling-depth" es el número de secuencias para enrarecer. Si alguna de las muestras 
#tiene un número menor de secuencias que las enrarecidas, será excluida del análisis. 

#creo visualizador de la tabla relativizada
qiime feature-table summarize \
--i-table salidas/core-metrics-results_KitsE/rarefied_table.qza \
--o-visualization salidas/core-metrics-results_KitsE/rarefied_table.qzv \
--m-sample-metadata-file metadata.tsv
#visualizo
qiime tools view salidas/core-metrics-results_KitsE/rarefied_table.qzv 


### hago la curva de rerefaccion con todas las ASVs sobre los datos que ahora se encuentran relativizados, para definir cual fue el kit de extraccion que mas diversiadd (ASVs rescato)
qiime diversity alpha-rarefaction \
--i-table salidas/core-metrics-results_KitsE/rarefied_table.qza \
--p-max-depth 10000 \
--m-metadata-file metadata.tsv \
--p-steps 25 \
--o-visualization salidas/alpha_rarefaction_KitsE.qzv
#visualizo
--o-visualization salidas/alpha_rarefaction_KitsE.qzv

#visualizo
qiime tools view salidas/alpha_rarefaction_KitsE.qzv
#Aunque el kit de suelos (Powersoil) haya rescatado mayor numero de secuencias y las asignaciones de peces se hallan dado en una muestra que fue extraida con dicho kit,
#Se puede ver claramente que el kit de aguas (KM) rescata mayor diversidad en general, una vez que se corrige las disparidades dadas por la profundidad de muestreo desigual de las muestras
####IMPORTANTE: estas metricas de diversidad estan sobreestimadas ya que tienen en cuenta toda la diversidad rescatada por los primers, lo cual se aleja de nuestro objetivo
#que es el calcular las metricas de la diverdad de peces, por eso de aqui en adelante se calculan con los datos filtrados
#**************************************************************************************************************************



#corro (ejecuto) las métricas de diversidad predeterminadas  : 
#se debe comparar un numero de secuencias que considere el mismo numero para cada una de las muestras, por ende el valor maximo que puede tener
#el argumento --p-sampling-depth sera el de la muestra con menos secuencias (features)
qiime diversity core-metrics-phylogenetic \
--i-phylogeny salidas/filogenia_ambiental_filtrada/rooted_tree.qza \
--i-table salidas/tabla_filtrada.qza \
--p-sampling-depth 204 \
--m-metadata-file metadata.tsv \
--output-dir salidas/core-metrics-results_filtrada
#Plugin error from diversity:
#
#  Ordinations with less than two dimensions are not supported.
#
#Debug info has been saved to /tmp/qiime2-q2cli-err-a4i0hxdp.log


#genero la curva de rarefaccion, usando la informacion proveeida por la tabla para escoger el valor de "--p-max-depth"
qiime diversity alpha-rarefaction \
--i-table salidas/tabla_filtrada.qza \
--p-max-depth 204 \
--m-metadata-file metadata.tsv \
--p-steps 25 \
--o-visualization salidas/alpha_rarefaction.qzv

#visualizo
qiime tools view salidas/alpha_rarefaction.qzv


#ya que solo hay una muestra tengo que sacar las metricas que pueda por aparte
#creo una carpeta donde estaran las metricas
mkdir salidas/core-metrics-results

#*******************diversidad filogenetica (faith)
#creo vizualizador para ver los valores del indice de diversidad alfa escojido, en cada muestra
qiime diversity alpha-phylogenetic \
--i-table salidas/tabla_filtrada.qza \
--i-phylogeny salidas/filogenia_ambiental_filtrada/rooted_tree.qza \
--p-metric "faith_pd" \
--o-alpha-diversity salidas/core-metrics-results/faith_pd_vector.qza
#creo vizualizador para ver los valores del indice de diversidad alfa escojido
qiime metadata tabulate \
--m-input-file salidas/core-metrics-results/faith_pd_vector.qza \
--o-visualization salidas/core-metrics-results/faith_pd_vector.qzv
#vizualizo
qiime tools view salidas/core-metrics-results/faith_pd_vector.qzv


#*******************diversidad de especies (shannon)
qiime diversity alpha \
--i-table salidas/tabla_filtrada.qza \
--p-metric "shannon" \
--o-alpha-diversity salidas/core-metrics-results/shannon_vector.qza
#creo vizualizador para ver los valores del indice de diversidad alfa escojido, en cada muestra
qiime metadata tabulate \
--m-input-file salidas/core-metrics-results/shannon_vector.qza \
--o-visualization salidas/core-metrics-results/shannon_vector.qzv
#vizualizo
qiime tools view salidas/core-metrics-results/shannon_vector.qzv



#*******************equidad de Pielou (evenness) (pielou)
qiime diversity alpha \
--i-table salidas/tabla_filtrada.qza \
--p-metric "pielou_e" \
--o-alpha-diversity salidas/core-metrics-results/evenness_vector.qza
#creo vizualizador para ver los valores del indice de diversidad alfa escojido, en cada muestra
qiime metadata tabulate \
--m-input-file salidas/core-metrics-results/evenness_vector.qza \
--o-visualization salidas/core-metrics-results/evenness_vector.qzv
#vizualizo
qiime tools view salidas/core-metrics-results/evenness_vector.qzv



#*******************diversidad de simpson (simpson)
qiime diversity alpha \
--i-table salidas/tabla_filtrada.qza \
--p-metric "simpson" \
--o-alpha-diversity salidas/core-metrics-results/simpson_vector.qza
#creo vizualizador para ver los valores del indice de diversidad alfa escojido, en cada muestra
qiime metadata tabulate \
--m-input-file salidas/core-metrics-results/simpson_vector.qza \
--o-visualization salidas/core-metrics-results/simpson_vector.qzv
#vizualizo
qiime tools view salidas/core-metrics-results/simpson_vector.qzv


#*******************diversidad de Chao (chao1)
qiime diversity alpha \
--i-table salidas/tabla_filtrada.qza \
--p-metric "chao1" \
--o-alpha-diversity salidas/core-metrics-results/chao1_vector.qza
#creo vizualizador para ver los valores del indice de diversidad alfa escojido, en cada muestra
qiime metadata tabulate \
--m-input-file salidas/core-metrics-results/chao1_vector.qza \
--o-visualization salidas/core-metrics-results/chao1_vector.qzv
#vizualizo
qiime tools view salidas/core-metrics-results/chao1_vector.qzv



#**************************************************************DIVERSIDAD BETA


#*******************distancia Unifrac ponderada (weighted_unifrac)
qiime diversity beta \
--i-table salidas/tabla_filtrada.qza \
--p-metric "braycurtis" \
--o-distance-matrix salidas/core-metrics-results/bray_curtis_pcoa_results.qza

qiime diversity beta \
--i-table salidas/tabla_filtrada.qza \
--p-metric "jaccard" \
--o-distance-matrix salidas/core-metrics-results/jaccard_pcoa_results.qza


####################################################################################################################################################################################################################################
######################################### exportacion datos #########################################
#exporto la tabla de frecuencias final en formato BIOM
qiime tools export \
  --input-path salidas/tabla.qza \
  --output-path tabla_biom_exportada
  
#convierto la tabla en formato .tsv
biom convert \
  -i tabla_biom_exportada/feature-table.biom \
  -o tabla_biom_exportada/feature-table.tsv \
  --to-tsv
  
#importo la taxonomia asignada  
qiime tools export \
  --input-path salidas/vsearch-taxonomy_SuppMaterial2.qza \
  --output-path tabla_biom_exportada/  

#añado la taxonomia a la tabla
biom add-metadata \
  -i tabla_biom_exportada/feature-table.biom \
  -o tabla_biom_exportada/table-with-taxonomy.biom \
  --observation-metadata-fp tabla_biom_exportada/taxonomy.tsv \
  --sc-separated taxonomy

#convierto la tabla en formato .tsv
biom convert \
  -i tabla_biom_exportada/table-with-taxonomy.biom \
  -o tabla_biom_exportada/table-with-taxonomy.tsv \
  --to-tsv
