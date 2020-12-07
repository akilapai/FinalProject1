#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)

### James Function Source
source("James_Func.R")

###Gabrielle data
load("Part2ResultsAllLegaspi.rda")

###Girish and Aislyn data
load("FinalData.rda")
source("binscatter.R")

### Zelong data
load("Part2ResultsALLZelong.rda")
graphdata<-ResultALL

###Akila data
Graph1Data<-subset(FinalData, 
                   ProteinName=="TNFRSF9" | 
                     ProteinName=="NT-3" |
                     ProteinName=="IL-22RA1" |
                     ProteinName=="CCL4")

# Aislyn's Setup --------------------------------------------------------

DrugA_aislyn<-subset(FinalData, 
                     Treatment=="DrugA")
DrugB_aislyn<-subset(FinalData, 
                     Treatment=="DrugB")

aislyn_sds <- function(x) {
  mean <- mean(x)
  sd1 <- mean+sd(x)
  sd2 <- mean-sd(x)
  return(c(y=mean, ymin=sd2, ymax=sd1))
}

## JANE'S DATA
JaneData<-subset(FinalData, 
                 ProteinName=="TNFRSF9" | 
                   ProteinName=="NT-3" |
                   ProteinName=="IL-22RA1" |
                   ProteinName=="CCL4"|
                   ProteinName== "VEGF-A"|
                   ProteinName == "BDNF" |
                   ProteinName=="MCP-3" |
                   ProteinName== "hGDNF"|
                   ProteinName== "IL-8"|
                   ProteinName== "CDCP1")

JaneDataLess <- subset(JaneData, Age<=40)
JaneDataOver <-subset(JaneData, Age>40 )
ExpressionLess <- JaneDataLess$Expression
ExpressionOver <- JaneDataOver$Expression

# Define server logic required to draw a histogram
shinyServer(function(input, output){
  
  ###Gabrielle's section###
  gabi <- ResultsAll
  output$gabrielle <-renderPlot({
    if(input$predictor == "Time"){hist(gabi[[3]], main = "Histogram of P-values for Protein Expression by Wk0 vs Wk4", 
                                       xlab = "P-values", ylim = c(0,50), cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.6, 
                                       col = "skyblue", breaks=seq(0,1,0.05), xaxp=c(0,1,20), labels = TRUE)
      abline(v = 0.05, col="red", lwd=1, lty=2)
      
      a <- NULL
      for(i in which(gabi[[3]]<=0.05)){
        a <- rbind(a,rownames(gabi)[i])}
      colnames(a) <- "p <= 0.05"
      output$names <- renderTable({a})}
    
    if(input$predictor == "Treatment"){hist(gabi[[6]], main = "Histogram of P-values for Protein Expression by DrugA vs DrugB", 
                                            xlab = "P-values", ylim = c(0,30), cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.6, 
                                            col = "lightgreen", breaks=seq(0,1,0.05), xaxp=c(0,1,20), labels = TRUE)
      abline(v = 0.05, col="red", lwd=1, lty=2)
      
      b <- NULL
      for(i in which(gabi[[6]]<=0.05)){
        b <- rbind(b,rownames(gabi)[i])}
      colnames(b) <- "p <= 0.05"
      output$names <- renderTable({b})}
    
    
    if(input$predictor == "Gender"){hist(gabi[[9]], main = "Histogram of P-values for Protein Expression by M vs F", 
                                         xlab = "P-values", ylim = c(0,40), cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.6, 
                                         col = "pink", breaks=seq(0,1,0.05), xaxp=c(0,1,20), labels = TRUE)
      abline(v = 0.05, col="red", lwd=1, lty=2)
      
      c <- NULL
      for(i in which(gabi[[9]]<=0.05)){
        c <- rbind(c,rownames(gabi)[i])
      }
      colnames(c) <- "p <= 0.05"
      output$names <- renderTable({c})}})
  ###END of Gabrielle's section###
  
  ####Girish's Section####
  
  sliderValues <- reactive({
    
    data.frame(
      Name = c("Integer"),
      Value = as.character(c(input$integer)),
      stringsAsFactors = FALSE)
    
  })
  
  
  output$girish <-renderPlot({
    
    FinalData_violin_gv <- FinalData %>% filter(ProteinName %in% c("MMP-10", "FGF-19", "CCL28",  "CST5",   "SCF",    "BDNF",  "TRANCE", "MCP-3",  "IL-6" ,  "FGF-21"))
    expr_range_gv<-range(FinalData_violin_gv$bmi)
    numbers_of_bins_gv <- input$integer
    finaldata_gv<-FinalData%>%mutate(MyQuantileBins = cut(bmi, 
                                                          breaks = seq(expr_range_gv[1],expr_range_gv[2],length.out = numbers_of_bins_gv), 
                                                          include.lowest=TRUE,labels =FALSE))#))#
    
    
    a_gv<-seq(expr_range_gv[1],expr_range_gv[2],length.out = numbers_of_bins_gv)
    b_gv<-a_gv[-length(a_gv)] + diff(a_gv)/2
    
    #ints <- findInterval(test$x, brks, all.inside = T)
    #(brks[ints] + brks[ints + 1]) / 2 
    
    finaldata_gv_bins<-finaldata_gv %>%
      group_by(MyQuantileBins,ProteinName) %>%
      mutate(meanexpr = mean(Expression))
    
    finaldata_gv_bins_unique<- finaldata_gv_bins %>% 
      filter( !duplicated(paste0(ProteinName,MyQuantileBins)))
    finaldata_gv_bins_unique2 <-finaldata_gv_bins_unique %>% mutate(bmi2 = b_gv[MyQuantileBins])
    
    finaldata_gv_bins_unique2_3 <- finaldata_gv_bins_unique2 %>% filter(ProteinName %in% c("MMP-10", "FGF-19", "CCL28",  "CST5",   "SCF",    "BDNF",  "TRANCE", "MCP-3",  "IL-6" ,  "FGF-21"))
    
    
    
    ggplot(finaldata_gv_bins_unique2_3, aes(y=meanexpr, x=bmi2, color=ProteinName)) + geom_point()+geom_smooth(se=FALSE) +labs(x="BMI",y="Expression",title="Expression vs. BMI for the 5 highest and lowest expression~BMI coefficients", subtitle = "Mean expression by BMI bin with loess smoothed lines")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
  })
  
  #output$girish <-renderPlot(
  #  ggplot(FinalData_violin, aes(y=Expression, x=bmi, color=ProteinName)) + geom_point() )
  
  #+ scatter.thinning(bmi,Expression, resolution=100,max.per.cell=100)
  
  #binscatter(formula="Expression ~ bmi", key_var = "bmi",
  #   data=FinalData_violin, bins=1000, partial=FALSE)
  
  
  ####END Girish's Section####
  
  
  ###Zelong's Section###
  
  time1<-vector()
  time1[1:92]<-"week0"
  time1[93:184]<-"week4"
  graphdata$time<-time1
  output$distPlot<-renderPlot({ggplot(graphdata[graphdata$time==input$Time,],aes(x=ProteinName, y=Expression_change))+geom_point(aes(size=3, color = Significance)) +
      scale_color_manual(labels=c("insignificant","significant"),values = c("black","red"))+labs(x="Protein Name", y=" Expression Change(DrugA-DrugB)", title = "Scatterplot of Protein Expression Change of different drug treatments")+
      theme(plot.title = element_text(size=24,color="orange"),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))},
      height = 600, width = 800)
  
  ###End Zelong's Section###
  
  
  ###Akila's Section###
  output$mygraph <-renderPlot(ggplot(Graph1Data[Graph1Data$Gender == input$Gender,], aes(x=Time, y=Expression, fill=Time))+
                                geom_violin() + geom_boxplot(width=0.1, fill="white") +
                                labs(title="Protein Expression vs. Treatment and Time") +
                                facet_wrap(~ProteinName+Treatment, nrow = 2, scales = "free_y") + 
                                theme_classic() + scale_fill_brewer(palette="PuRd"))
  ###END Akila's Section###
  
  ###Parker's Section###
  plot3 <- subset(FinalData, Treatment=="DrugA")
  plot3 <- subset(plot3, ProteinName=="VEGF-A")
  plot3 <- sample_n(plot3, 200)
  plot6 <- subset(FinalData, Treatment=="DrugB")
  plot6 <- subset(plot6, ProteinName=="VEGF-A")
  plot6 <- sample_n(plot6, 200)
  output$plotgraph1 <- renderPlot({
    x    <- plot3$bmi
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#75AADB", border = "white",
         xlab = "BMI's for Treatment A",
         main = "Histogram for BMI")})
  
  output$plotgraph2 <- renderPlot({
    y   <- plot6$bmi
    bins <- seq(min(y), max(y), length.out = input$bins + 1)
    
    hist(y, breaks = bins, col = "#75AADB", border = "white",
         xlab = "BMI's for Treatment B",
         main = "Histograms for BMI")})
  
  ###END Parker's Section###
  
  # Aislyn's Section --------------------------------------------------------
  
  output$aislyn_plot <-renderPlot(
    ggplot(FinalData[FinalData$ProteinName == input$ProteinName,], aes(x=Time, y=Expression))+
      geom_violin(aes(fill=Time)) + 
      ggtitle(unique(input$ProteinName)) +
      scale_fill_brewer(palette=2, labels=c("Week 0", "Week 4")) +
      facet_wrap(~Treatment) +
      labs(fill="Timepoint", x="Time", y="Expression") +
      stat_summary(fun.data=aislyn_sds, color="black", geom="pointrange", shape=5, size=.5) +
      theme_bw(base_size=18) + theme(legend.title = element_text(face="bold"), legend.text = element_text(size=18), axis.text.x = element_blank(), legend.position = "right",  axis.title.y=element_text(face="bold"))
  )
  
  output$aislyn_title <- renderText({
    aislyn_1 = "Comparing Expression at Week 0 vs Week 4: "
    paste(aislyn_1)
    
  })
  output$aislyn_drugA <- renderText({
    aislyn_2 = "  Drug A (p-value):  " 
    paste0(aislyn_2, 
           paste(
             format(
               t.test(DrugA_aislyn[DrugA_aislyn$ProteinName == input$ProteinName,]$Expression ~ DrugA_aislyn[DrugA_aislyn$ProteinName == input$ProteinName,]$Time, paired=TRUE)$p.value, 
               digits=4)), 
           sep = " ")
    
  })
  output$aislyn_drugB <- renderText({
    aislyn_3 = "  Drug B (p-value): "
    paste0(aislyn_3, 
           paste(
             format(
               t.test(DrugB_aislyn[DrugB_aislyn$ProteinName == input$ProteinName,]$Expression ~ DrugB_aislyn[DrugB_aislyn$ProteinName == input$ProteinName,]$Time, paired=TRUE)$p.value, 
               digits=4)), 
           sep = " ")
    
  })
  
  
  
  
  
  ### START JANE'S SECTION###
  output$jane<-renderPlot ({
    # if (input$Age == "Less than or Equal to 40")
    # {ggplot(JaneDataLess, aes(x=ProteinName, y=ExpressionLess))+
    #     geom_point() +
    #     labs(title="Protein vs. Expression", xlab="Protein", ylab="Expression")}
    # 
    
    
    if (input$Age == "Over 40")
    {ggplot(JaneDataOver, aes(x=ProteinName, y=ExpressionOver))+
        geom_point() +
        labs(title="Protein vs. Expression", xlab="Protein", ylab="Expression")}
    
    else if (input$Age == "Less than or Equal to 40")
    {ggplot(JaneDataLess, aes(x=ProteinName, y=ExpressionLess))+
        geom_point() +
        labs(title="Protein vs. Expression", xlab="Protein", ylab="Expression")}
    
    
    
    })
  
  
  ## End Jane's Section 

  
  ### START Wendy'S SECTION###
  output$wendy <- renderPlot({
    if(input$gt == "Males on Drug A"){
      final_wendy<- FinalData %>%
        distinct(PatientID,.keep_all = TRUE) %>%
        filter(Gender =="M") %>%
        filter(Treatment =="DrugA")
      hist(final_wendy$Age, main="Age Distribution among Males on Drug A", ylab= "Frequency",xlab="Age", col="lightblue")}
    
    if(input$gt == "Males on Drug B"){
      final_wendy<- FinalData %>%
        distinct(PatientID,.keep_all = TRUE) %>%
        filter(Gender =="M") %>%
        filter(Treatment =="DrugB")
      hist(final_wendy$Age, main="Age Distribution among Males on Drug B", ylab= "Frequency",xlab="Age", col="lightblue")}
    
    if(input$gt == "Females on Drug A"){
      final_wendy<- FinalData %>%
        distinct(PatientID,.keep_all = TRUE) %>%
        filter(Gender =="F") %>%
        filter(Treatment =="DrugA")
      hist(final_wendy$Age, main="Age Distribution among Females on Drug A", ylab= "Frequency",xlab="Age", col="lightblue")}
    
    if(input$gt == "Females on Drug B"){
      final_wendy<- FinalData %>%
        distinct(PatientID,.keep_all = TRUE) %>%
        filter(Gender =="F") %>%
        filter(Treatment =="DrugB")
      hist(final_wendy$Age, main="Age Distribution among Females on Drug B", ylab= "Frequency",xlab="Age", col="lightblue")}
  })
  
  ### Start James
  output$AllProt.JL <- renderPlot({
    in_SIZE = input$SIZE.JL
    in_MOD = input$MOD.JL
    data.JL <- James_Func(DATA = FinalData, MOD = in_MOD, SIZE = in_SIZE)[,2]
    bins <- seq(from=min(data.JL), to=max(data.JL), length.out=100)
    hist(data.JL, breaks = bins, col = rgb(.2,1,1), border = 'white', 
         main = paste0('Bootstrapped Regression Slopes of Expression by ', in_MOD),
         xlab = 'Regression slopes',
         ylab = 'Density')})
  
  output$SingleProt.JL <- renderPlot({
    in_SIZE = input$SIZE.JL
    in_MOD = input$MOD.JL
    in_prot = input$protein_number.JL
    data_sub.JL = subset(FinalData,ProteinNumber == in_prot)
    data.JL <- James_Func(DATA = data_sub.JL,
                          MOD = in_MOD, SIZE = in_SIZE)[,2]
    bins <- seq(from=min(data.JL), to=max(data.JL), length.out=20)
    hist(data.JL, breaks = bins, col = 'red', border = 'white', 
         main = paste0('Bootstrapped Regression Slopes of ',
                       in_prot,
                       ' Expression by ',
                       in_MOD),
         xlab = 'Regression slopes',
         ylab = 'Density',
         ylim = c(0,5))})
  
  } ) 





