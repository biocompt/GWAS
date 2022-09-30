library(ramwas)

setwd('/hdd/Estudios_finalizados/gwas/')
ruta <- getwd()
gwas <- list.files(ruta, pattern = "*.txt")

for (i in 1:length(gwas)){
  resultados <- read.table(gwas[i], header = T, sep = "\t")
  head(resultados)
  resultados <- na.omit(data.frame(resultados))
  resultados <- resultados[order(resultados$P),]
  head(resultados)
  resultados$CHR <- as.factor(resultados$CHR)
  
  # ---- Hacemos el manhattan plot ----
  max_pvalor <- -log10(min(resultados$P)/100)
  if (max_pvalor<=-log10(1e-10)){
    max_pvalor <- -log10(1e-10)
  }
  
  tiff(paste("plots/", sub(".txt", ".jpeg", gwas[i]), sep =""),
       width = 1200, height = 600, units = 'px')
  manPlotFast(manPlotPrepare(resultados$P, resultados$CHR, resultados$BP), ylim=c(0,max_pvalor), cex=0.5)
  abline(h=c(-log10(1e-8), -log10(1e-5)), col=c('blue','red'), lty=c(2,2))
  dev.off()
  
  input <- resultados
  stat_type = "PVAL"
  color <- "#BB0000"
  xmax <- as.numeric(10)
  output <- paste("plots/QQplot_", sub(".txt", ".png", gwas[i]), sep ="")
  source("/home/carlos/scripts/qqplot_function.R")
  rm(z)
  dev.off()
}

