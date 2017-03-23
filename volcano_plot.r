library("ggplot2")
library("dplyr")
library("ggrepel")
library("reshape2")

volcano.plot<-function(df, log.fold.change = TRUE) {
  print(head(df))
  if (log.fold.change == TRUE) {
  p<- ggplot(data = df, (aes(log2(Fold.change), -log(padj) )))
  p+
  geom_point(data = df, aes(colour = 'grey'), size = rel(4))+
  geom_hline(yintercept = -log(0.05), linetype = 2) + 
  geom_vline(xintercept = 1, linetype =2) +
  geom_vline(xintercept = -1, linetype =2) +
  labs(x = "log2(Fold Change)", y = "-log(FRD Adjusted p-value)") +
  theme(axis.text =     element_text(size = rel(1)),
	axis.title.x =  element_text(size = rel(1.25)), 
	axis.title.y =  element_text(size = rel(1.25)),
	legend.text =   element_text(size = rel(1)), 
	legend.title =  element_text(size = rel(1.25)),
	legend.position = c(0.15,0.85)
	) +
  geom_text_repel( data = filter(df, 
		  (padj < 0.05 | log2(Fold.change) > 1 | log2(Fold.change) < -1)),
		  aes(label = Metabolite), 
	size = rel(4),
	box.padding = unit(0.35, "lines"),
	point.padding = unit(0.3,"lines"))
  } else {
  
  p<- ggplot(data = df, (aes((Fold.change), -log(padj) )))
  p+
  geom_point(data = df, aes(colour = 'grey'), size = rel(4))+
  geom_hline(yintercept = -log(0.05), linetype = 2) + 
  geom_vline(xintercept = 1, linetype =2) +
  geom_vline(xintercept = -1, linetype =2) +
  labs(x = "Fold Change", y = "-log(FRD Adjusted p-value)") +
  theme(axis.text =     element_text(size = rel(1)),
	axis.title.x =  element_text(size = rel(1.25)), 
	axis.title.y =  element_text(size = rel(1.25)),
	legend.text =   element_text(size = rel(1)), 
	legend.title =  element_text(size = rel(1.25)),
	legend.position = c(0.15,0.85)
	) +
  geom_text_repel( data = filter(df, 
		  (padj < 0.05 | (Fold.change) > 1.5 | (Fold.change) < -1.5)),
		  aes(label = Metabolite), 
	size = rel(4),
	box.padding = unit(0.35, "lines"),
	point.padding = unit(0.3,"lines"))
          
  }
          
}

volcano.plot<-function(df, log.fold.change = TRUE) {
  print(head(df))
  p<- ggplot(data = df, (aes(log2(Fold.change), -log(padj) )))
  p+
  geom_point(data = df, aes(colour = 'grey'), size = rel(4))+
  geom_hline(yintercept = -log(0.05), linetype = 2) + 
  geom_vline(xintercept = 1, linetype =2) +
  geom_vline(xintercept = -1, linetype =2) +
  labs(x = "log2(Fold Change)", y = "-log(FRD Adjusted p-value)") +
  theme(axis.text =     element_text(size = rel(1)),
	axis.title.x =  element_text(size = rel(1.25)), 
	axis.title.y =  element_text(size = rel(1.25)),
	legend.text =   element_text(size = rel(1)), 
	legend.title =  element_text(size = rel(1.25)),
	legend.position = c(0.15,0.85)
	) +
  geom_text_repel( data = filter(df, 
		  (padj < 0.05 | log2(Fold.change) > 1 | log2(Fold.change) < -1)),
		  aes(label = Metabolite), 
	size = rel(4),
	box.padding = unit(0.35, "lines"),
	point.padding = unit(0.3,"lines"))
          
}
