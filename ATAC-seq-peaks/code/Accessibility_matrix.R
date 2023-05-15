#Figure 3E in paper

#setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2")
#system("gunzip 2E_Accessibility_score.tsv.gz")
#setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject")
library(readr)
Accessibility_scores=read_tsv(file=paste0(scratch,"CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Accessibility_score.tsv"))

coordinates=Accessibility_scores$...1
cell_lines=names(Accessibility_scores)

mat=as.matrix(Accessibility_scores[,-1])
rownames(mat)=coordinates

library("ComplexHeatmap")
library("viridis")
pdf(output_hs, width=8.0, height=7)
hm <- Heatmap( mat,
               col=viridis(99),
               cluster_rows=T, cluster_columns=T,
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "ward.D",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "ward.D",
               show_row_names=F,
               show_column_names=T, column_names_side="top",
               column_names_gp = gpar(fontsize=5),
               heatmap_legend_param = list(
                 title = "mean(log(1+acc.))",
                 title_gp = gpar(fontsize = 6),
                 title_position = "leftcenter-rot",
                 labels_gp = gpar(fontsize = 5)
               )
)
draw(hm)




# tableFl <- lookupConfig "accessibility" "/projects/ps-renlab/kai/project/Atlas/output/SCATACSeq/Feature/Peak/Cluster/relative_accessibility_scores.tsv"
# dir <- lookupConfig "output_dir" "output/"
# orderFl <- lookupConfig "cluster_order" "../Cluster/output/phylo.txt"
# annoFl <- lookupConfig "cluster_annotation" "/projects/ren-transposon/home/kai/Atlas/annotation.tsv" 
# let output = dir <> "/restricted_peaks.pdf"
# liftIO $ do
# orders <- reverse . T.lines <$> T.readFile orderFl
# anno <- readAnno annoFl
# peaks <- forM input $ \(_, peakFl) ->
# map (T.pack . B.unpack . showBed) <$> (readBed peakFl :: IO [BED3])
# print $ length $ nubSort $ concat peaks
# peaks' <- create >>= sampling' 5000 peaks
# df <- DF.map (logBase 2 . (+1)) . (`DF.csub` orders) <$> DF.readTable tableFl
# plotRestrictedPeaks output $ changeName anno $
# DF.rbind $ diagonize' average $ 
#   map (DF.map (min 6) . (df `DF.rsub`)) peaks'
#                               {-
#                                   let tables = V.fromList $ parMap rdeepseq (getSignal df) peaks
#                                   (names, vecs) = unzip $ flatten $ hclust Ward tables (euclidean `on` snd)
#                                   df' = DF.mkDataFrame names (DF.colNames df) $ map V.toList vecs
#                 num_peak = map (\(x,y) -> (x, length y)) peaks
# 
#             plotRestrictedPeaks (dir <> "/restricted_peaks.pdf") num_peak $ df' `DF.csub` names
#                                   
#                                   DF.writeTable output (T.pack . show) $ df' `DF.csub` names
#             return (output, num_peak)
#             -}
#         |] $ return ()

hclust Ward tables (euclidean `on` 
                    
                    plotRestrictedPeaks output df = R.runRegion $ do
                    mat <- toRMatrix df
                    
                    library("ComplexHeatmap")
                    library("viridis")
                    pdf(output_hs, width=8.0, height=7)
                    hm <- Heatmap( mat_hs,
                                   col=viridis(99),
                                   cluster_rows=F, cluster_columns=F,
                                   show_row_names=F,
                                   show_column_names=T, column_names_side="top",
                                   column_names_gp = gpar(fontsize=5),
                                   heatmap_legend_param = list(
                                     title = "mean(log(1+acc.))",
                                     title_gp = gpar(fontsize = 6),
                                     title_position = "leftcenter-rot",
                                     labels_gp = gpar(fontsize = 5)
                                   )
                    )
                    draw(hm)
                    dev.off()                    