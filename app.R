library(shiny)
library(shinythemes)
library(htmltools)
library(Biostrings)
library(DT)
library(data.table)
library(stringr)
library(shinyjs)
library(shiny)
library(ggmsa)
library(pwalign)
library(ggplot2)
library(dplyr)
library(Rsamtools)

options(shiny.maxRequestSize=100*1024^2)
data_table1 <- read.table('www/table1.txt',header = T,sep = '\t')

Orthogroups.tbl=data.table::fread('www/strainseq/Orthogroups.txt',header = T,sep = '\t')


RSAG1.stat=read.table('www/strainseq/RSAG1.stat',header = T,sep = '\t')


geneinfo=read.table('www/strainseq/geneinfo.txt',header = T,sep = '\t',quote = '')




source("source.R")


# Define UI
ui <- fluidPage(

  theme = shinytheme('flatly'),
  
  
  
  
  h1(strong("立枯丝核菌AG1-IA基因组数据平台"), align = "center"),
  
  
  navbarPage(
    'V1.0',
    
 ########################## UI 1 ######################################################   
    tabPanel(
      'Home',
            fluidRow( 
               div(style = "border: 2px solid black; padding: 10px; border-radius: 5px;",
               div(style = "width: 100%; text-align: left;", 
                   p( tags$span(style = "font-size: 20px; color: black;","Rice sheath blight, caused by Rhizoctonia solani anastomosis group 1-intraspecific group IA (AG1-IA), is a devastating disease threatening global rice production. While several R. solani AG1-IA genomes are available, limitations in existing annotations hinder our understanding of its pathogenicity mechanisms. Here, we present a high-quality genome annotation of R. solani AG1-IA strain YZ19, a virulent rice sheath blight pathogen, generated using full-length transcriptome data and rigorous manual curation. Our analysis identified 37,417 transcripts at 12,464 gene loci, including 10,815 protein-coding genes with 25,834 unique open reading frames (ORFs). Notably, we discovered 1,840 genes containing micro-exons, a feature previously overlooked in R. solani annotations. Time-course transcriptome analysis revealed that approximately 50% of pathogenicity-related genes, including those encoding secretory proteins, putative CAZymes, and effectors, harbor multiple ORFs and micro-exons. Furthermore, we identified 2,983 lncRNAs, with 325 showing co-localization and transcriptional correlation with protein-coding genes during rice infection, suggesting a role in regulating virulence. A web-based genome database has been developed based on this enhanced annotation to facilitate gene discovery and comparative genomic analyses. Our findings provide a valuable resource for elucidating the molecular basis of R. solani AG1-IA pathogenicity and developing effective disease control strategies.
"))
               ),

               )
        ),
      hr(),
      fluidRow(
        tags$style(HTML("
        #my_table td, #my_table th {
            white-space: nowrap; /* 防止换行 */
            overflow: hidden;    /* 超出部分隐藏 */
            text-overflow: ellipsis; /* 使用省略号表示截断 */
            text-align: right; /* 表格内容右对齐 */
        }
    ")),
        column(6,
               div(style = "padding: 5px;", 
                   div(style = "width: 100%; display: flex; justify-content: center; ", 
                       img(src = "fig1.png", style = "max-width: 100%; height: auto;")),
        )),

        column(6,
               div(style = "padding: 5px;", 
               div(style = "width:100%; display: flex; justify-content: center;",  
                   div(style = "overflow-x: auto; width: auto;",  
                       br(),
                       p(tags$span(style = "font-size: 20px; color: black;",
                                   "Genome statistics of the 12 rice-infecting Rhizoctonia solani AG1-IA")),
                       tableOutput("my_table")
                   )
               )
        ))
        
      ),

      
      tags$div(
        hr(), 
        h2('Reference'),
        
        p( tags$span(style = "font-size: 20px; color: black;","HG81")),
        a(  tags$span(style = "font-size: 16px; color: red;", "Yang, Q., Yang, L., Wang, Y., Chen, Y., Hu, K., Yang, W., Zuo, S., Xu, J., Kang, Z., Xiao, X., & Li, G. (2022). 
                     A High-Quality Genome of Rhizoctonia solani, a Devastating Fungal Pathogen with a Wide Host Range. Molecular Plant-Microbe Interactions®, 35(10), 954–958."),
            href = "https://doi.org/10.1094/MPMI-06-22-0126-A",
            target = "_blank"),
        br(),
        br(),
        
        p( tags$span(style = "font-size: 20px; color: black;","BRS1")),
        a(  tags$span(style = "font-size: 16px; color: red;", "Francis, A., Ghosh, S., Tyagi, K., Prakasam, V., Rani, M., Singh, N. P., Pradhan, A., Sundaram, R. M., Priyanka, C., Laha, G. S., Kannan, C., Prasad, M. S., Chattopadhyay, D., & Jha, G. (2023). 
                      Evolution of pathogenicity-associated genes in Rhizoctonia solani AG1-IA by genome duplication and transposon-mediated gene function alterations. BMC Biology, 21(1)."),
            href = "https://doi.org/10.1186/s12915-023-01526-0",
            target = "_blank"),
        br(),
        br(),   
        
        p( tags$span(style = "font-size: 20px; color: black;","XN")),
        a(  tags$span(style =  "font-size: 16px; color: red;", "Li, C., Guo, Z., Zhou, S., Han, Q., Zhang, M., Peng, Y., Hsiang, T., & Chen, X. (2021). 
                      Evolutionary and genomic comparisons of hybrid uninucleate and nonhybrid Rhizoctonia fungi. Communications Biology, 4(1)."),
            href = "https://doi.org/10.1038/s42003-021-01724-y",
            target = "_blank"),
        br(),
        br(),
        
        p( tags$span(style = "font-size: 20px; color: black;","BM1")),
        a(  tags$span(style = "font-size: 16px; color: red;", "Kaushik, A., Roberts, D. P., Ramaprasad, A., Mfarrej, S., Nair, M., Lakshman, D. K., & Pain, A. (2022). 
                      Pangenome Analysis of the Soilborne Fungal Phytopathogen Rhizoctonia solani and Development of a Comprehensive Web Resource: RsolaniDB. Frontiers in Microbiology, 13."),
            href = "https://doi.org/10.3389/fmicb.2022.839524",
            target = "_blank"),
        br(),
        br(),  
        
        p( tags$span(style = "font-size: 20px; color: black;","B2 ADB WGL YN-7")),
        a(  tags$span(style = "font-size: 16px; color: red;", "Lee, D.-Y., Jeon, J., Kim, K.-T., Cheong, K., Song, H., Choi, G., Ko, J., Opiyo, S. O., Correll, J. C., Zuo, S., Madhav, S., Wang, G.-L., & Lee, Y.-H. (2021). 
                      Comparative genome analyses of four rice-infecting Rhizoctonia solani isolates reveal extensive enrichment of homogalacturonan modification genes. BMC Genomics, 22(1)."),
            href = "https://doi.org/10.1186/s12864-021-07549-7",
            target = "_blank"),
        br(),
        br(),
        
        p( tags$span(style = "font-size: 20px; color: black;","1802/KB")),
        a(  tags$span(style = "font-size: 16px; color: red;", "Nadarajah, K., Mat Razali, N., Cheah, B. H., Sahruna, N. S., Ismail, I., Tathode, M., & Bankar, K. (2017). 
                      Draft Genome Sequence of Rhizoctonia solani Anastomosis Group 1 Subgroup 1A Strain 1802/KB Isolated from Rice. Genome Announcements, 5(43)."),
            href = "https://doi.org/10.1128/genomea.01188-17",
            target = "_blank"),
        br(),
        br(),
        
        p( tags$span(style = "font-size: 20px; color: black;","GD118")),
        a(  tags$span(style = "font-size: 16px; color: red;", "Zheng, A., Lin, R., Zhang, D., Qin, P., Xu, L., Ai, P., Ding, L., Wang, Y., Chen, Y., Liu, Y., Sun, Z., Feng, H., Liang, X., Fu, R., Tang, C., Li, Q., Zhang, J., Xie, Z., Deng, Q., … Li, P. (2013). 
                      The evolution and pathogenic mechanisms of the rice sheath blight pathogen. Nature Communications, 4(1). "),
            href = "https://doi.org/10.1038/ncomms2427",
            target = "_blank"),
        br(),
        br()
      )
    ),
 ########################## UI 2 ######################################################   
 
    tabPanel('Blast', 
             HTML('<iframe src="http://124.222.22.158:4567" height="2000" width="100%"></iframe>')
    ),
    
    tabPanel('JBrowse',  
             p('Use JBrowse to view the genome sequences of Rhizoctonia solani AG-1A strain YZ19 and 
                the full-length transcript annotation. You can also view the coverage of Illumina strand-specific cDNA reads (BPM normalized) and 
               ONT cDNA reads (CPM normalized) at different hours post-infection, generated by poly-A selection.'),
             HTML('<iframe src="http://124.222.22.158:80/yn19web" height="2000" width="100%"></iframe>')
    ),
    
 ########################## UI 3 ######################################################   

    tabPanel('Orthologous genes & Sequence extraction',
             
      
             h3("Gene search"),
             
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 textAreaInput("search_id1", "Gene ID, transcript ID or Key words：",resize = 'vertical',
                               value = '',height = 200,
                               placeholder = 'examples:\nchr016g08480\nchr016g08490.mRNA1\nTPR-like\n...'),
                 actionButton("search1", "Search",icon = icon("magnifying-glass")),
                 actionButton("example_button1", "Example",icon = icon('pen')),  
                 h5('Note: One gene or keyword per line.')
               ),
               
               mainPanel(
                 width = 10,
                 verbatimTextOutput("error_message1"),  # 用于显示错误信息
                 DT::DTOutput("search_tbl")
               )
             ),
             
             hr(),
             
             ###
             h3("Orthologous genes"),
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 textInput("gene_id1", "Gene ID or transcript ID："),
                 actionButton("submit1", "Search",icon = icon("magnifying-glass")),
                 h5('Note: Multiple IDs separated by commas.')
               ),
               
               mainPanel(
                 width = 10,
                 verbatimTextOutput("error_message2"),  # 用于显示错误信息
                 DT::DTOutput("geneTable")
               )
             ),
             
             ###
             hr(),
             h3("Sequence extraction"),
             sidebarLayout(
               
               
               
               
               sidebarPanel(
                 width = 2,
                 selectInput("strain2", "Strain：",
                             choices = c("YZ19" = "YZ19", "XN" = "XN", "HG81" = "HG81","BRS1" = "BRS1",
                                         "B2" = "B2", "ADB" = "ADB", "WGL" = "WGL","YN-7" = "YN7",
                                         "GD118" = "GD118", "BM1" = "BM1")),
                 textInput("gene_id2", "Gene ID or Transcript ID：", ""),

                 sliderInput("upstream", "Upstream:",
                             min = 50, max = 5000,
                             value = 2000, step = 100,
                             post = "bp", sep = ",",
                             animate = F),
                 sliderInput("downstream", "Downstream:",
                             min = 50, max = 5000,
                             value = 2000, step = 100,
                             post = "bp", sep = ",",
                             animate = F),  
                 actionButton("submit2", "Submit",icon = icon("dna")),
                 h5('Note: The sequence located on the negative strand will be reverse complementary.')
               ),
               
               mainPanel(

                 
                 width = 10,
                 
                 fluidRow(
                   column(6, 
                          h4("Gene Sequence Output"),
                          verbatimTextOutput("sequence_out1")),
                   column(6,
                          h4("Transcript Sequence Output"),
                          verbatimTextOutput("sequence_out2"),
                          hr(),
                          h4("Peptide Sequence Output"),
                          verbatimTextOutput("sequence_out3"))
                 )
                 
               )
             )
    ),
 ########################## UI 4 ######################################################   
 
    tabPanel('Gene/Transcript Expression Query',
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                  textAreaInput(inputId = 'searchid2',label = 'Gene/Transcript ID',resize = 'vertical',
                                value = '',height = 500,placeholder = 'examples:\nchr01g00790\nchr01g00790.mRNA1\nchr01g00790.mRNA2...'),
                 actionButton("submit3", "Submit"),
                 actionButton("example_button2", "Example",icon = icon('pen')),  
                 p('Note: All the 3rd transcriptome data were normalized using CPM (Counts Per Million).')
               ),
               
               mainPanel(
                 width = 10,
                 verbatimTextOutput("error_message3"),
                 DT::DTOutput("gene_exp_table"),
                 hr(),
                 DT::DTOutput("transcript_exp_table")
                 
               )
             )
    ),
 ########################## UI 5 ######################################################   
 
    tabPanel('Download',
             p(tags$span(style = "font-size: 24px; color: black;","Download the genome sequences and associated transcript annotations of Rhizoctonia solani AG-1A")
               ),
             
             fluidRow(
               column(4, strong("Strain")),
               column(2, strong("Genome sequences (FASTA.gz)")),
               column(2, strong("Annotation (gff3.gz)")),
               column(2, strong("mRNA sequences (FASTA.gz)")),
               column(2, strong("Peptide sequences (FASTA.gz)"))
             ),
             hr(),
             fluidRow(
               column(4, "YZ19"),
               column(2, downloadButton("download1.1", "Download")),
               column(2, downloadButton("download1.2", "Download")),
               column(2, downloadButton("download1.3", "Download")),
               column(2, downloadButton("download1.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "XN"),
               column(2, downloadButton("download2.1", "Download")),
               column(2, downloadButton("download2.2", "Download")),
               column(2, downloadButton("download2.3", "Download")),
               column(2, downloadButton("download2.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "HG81"),
               column(2, downloadButton("download3.1", "Download")),
               column(2, downloadButton("download3.2", "Download")),
               column(2, downloadButton("download3.3", "Download")),
               column(2, downloadButton("download3.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "BRS1"),
               column(2, downloadButton("download4.1", "Download")),
               column(2, downloadButton("download4.2", "Download")),
               column(2, downloadButton("download4.3", "Download")),
               column(2, downloadButton("download4.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "B2"),
               column(2, downloadButton("download5.1", "Download")),
               column(2, downloadButton("download5.2", "Download")),
               column(2, downloadButton("download5.3", "Download")),
               column(2, downloadButton("download5.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "ADB"),
               column(2, downloadButton("download6.1", "Download")),
               column(2, downloadButton("download6.2", "Download")),
               column(2, downloadButton("download6.3", "Download")),
               column(2, downloadButton("download6.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "WGL"),
               column(2, downloadButton("download7.1", "Download")),
               column(2, downloadButton("download7.2", "Download")),
               column(2, downloadButton("download7.3", "Download")),
               column(2, downloadButton("download7.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "YN-7"),
               column(2, downloadButton("download8.1", "Download")),
               column(2, downloadButton("download8.2", "Download")),
               column(2, downloadButton("download8.3", "Download")),
               column(2, downloadButton("download8.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "GD118"),
               column(2, downloadButton("download9.1", "Download")),
               column(2, downloadButton("download9.2", "Download")),
               column(2, downloadButton("download9.3", "Download")),
               column(2, downloadButton("download9.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "BM1"),
               column(2, downloadButton("download10.1", "Download")),
               column(2, downloadButton("download10.2", "Download")),
               column(2, downloadButton("download10.3", "Download")),
               column(2, downloadButton("download10.4", "Download"))
             ),
             hr(),
             fluidRow(
               column(4, "1082/KB"),
               column(2, downloadButton("download11.1", "Download"))
             )
             ),
 ########################## UI 6 ######################################################   
 navbarMenu('Tools',
            tabPanel(
              "rename gff",
              #gff rename
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  fileInput('upload_gff', label = 'upload gff3 file(Max size 50MB, .gz Recommended)', accept = 'gff'),
                  verbatimTextOutput("fileContents"),
                  verbatimTextOutput("old_gff_info")
                ),
                mainPanel(DT::DTOutput("headContents"))
                
              ),
              hr(),
              sidebarLayout(
                sidebarPanel(width = 3,
                             uiOutput("new_gff_id_digits"),
                             uiOutput("new_gff_id_step"),
                             verbatimTextOutput('new_id_sample'),
                             uiOutput("overlap"),
                             uiOutput("after_upload_gff_actionbutton")
                             ),
                mainPanel(width = 9, fluidRow(
                  column(
                    3,
                    uiOutput("after_upload_gff_p2"),
                    uiOutput("after_upload_gff_newchr")
                  ),
                  column(
                    3,
                    uiOutput("after_upload_gff_p1"),
                    verbatimTextOutput("orgin_chr")

                  ), 
                  column(
                    6,
                    uiOutput("after_fixgff_title"),
                    DT::DTOutput("newgfftbl"),
                    hr(),
                    uiOutput('fixgff_download')
                  )
                ))
              )
              
            ),
            
            tabPanel("Pairwise Alignment",#Pairwise Alignment
                     sidebarLayout(
                       sidebarPanel(
                         textAreaInput("align_seq1", label = "sequence 1", 
                                       rows = 5, placeholder = "Nucleic acid or amino acid sequence."),
                         textAreaInput("align_seq2", label = "sequence 2", 
                                       rows = 5, placeholder = "Nucleic acid or amino acid sequence."),
                         
                         radioButtons("alignmentType", label = "Algorithm:",
                                     choices = list("global" = "global", "local" = "local")),
                         selectInput("substitutionMatrix", label = "substitutionMatrix(only work for aa sequence alignment):",
                                      choices = list("BLOSUM62" = "BLOSUM62", 
                                                     "BLOSUM50" = "BLOSUM50",
                                                     "BLOSUM80" = "BLOSUM80",
                                                     "BLOSUM100" = "BLOSUM100",
                                                     "PAM70" = "PAM30",
                                                     "PAM70" = "PAM70",
                                                     "PAM120" = "PAM120",
                                                     "PAM250" = "PAM250")),
                         numericInput('gap_open', label = 'gapOpening:', value = 2),
                         numericInput('gap_extension', label = 'gapExtension:', value = 5),
                         sliderInput('blockwidth','block width of alignment output',
                                     min = 50,max=300,step = 50,value = 100),
                         actionButton("align", "Align"),
                         hr(),
                         hr(),
                         uiOutput("after_align_range"),
                         fluidRow(
                           column(6,  uiOutput("after_align_button")),
                           column(6,  uiOutput("after_align_plot"))
                         )
                       ),
                       
                       mainPanel(
                        br(),
                           verbatimTextOutput("alignmentResult"),
                           plotOutput('alignmentPlot')
                         
                       )
                     )
                     ),
            
            tabPanel("ORF Prediction",
                     sidebarPanel(width = 3,
                                  textAreaInput("seq_need_ORF", label = "Input sequence:", 
                                                rows = 5, placeholder = "Nucleic acid sequence."),
                                  radioButtons("seq_need_ORF_strand", label = "Sequence strand:",
                                               choices = list("Forward(Orgin seq)" = F,
                                                              "Reverse(ReverseComplement seq)" = T)),
                                  sliderInput("orf_len", "Min Peptide length:",
                                              min = 10, max = 100,
                                              value = 30, step = 10,
                                              post = "aa", sep = ",",
                                              animate = F),
                                  actionButton("getORF", "get ORF")
                     ),
                     mainPanel(
                                  verbatimTextOutput("ORF_seq"),
                                  verbatimTextOutput("ORF_pep"))
                     ),
            tabPanel("Micro-exon identification",
                     sidebarPanel(width = 3,
                                  fileInput('upload_BAM_file',
                                            label = 'upload sorted BAM file', accept = 'bam'),
                                  fileInput('upload_fa_file',
                                            label = 'upload reference fasta file'),
                                  actionButton("getme", "Search for reads containing MEs")
                     ),
                     mainPanel(
                       DT::DTOutput("microexon_out")
                    )
            )
  )
 )
)
# Server logic
server <- function(input, output,session) {
  
  
  output$my_table <- renderTable({
    data_table1  # 返回数据框以显示在表格中
  }, rownames = TRUE,bordered = T)  # 可选择显示行名

  
  
#############1.gene search########### 
  observeEvent(input$example_button1, {
    # 当按钮被点击时，填充内容
    updateTextAreaInput(session, "search_id1", value = "chr016g08480\nchr016g08490.mRNA1\nTPR-like")
  })
  
  observeEvent(input$search1, {
    
    id<- trimws(input$search_id1)
    
    if (nchar(id) == 0) {
      output$error_message1 <- renderText("Please enter a gene ID or description.")
      return()  # 提前返回
    }
    
    # 读取用户输入并转为小写
    input_text <- tolower(id)
    input_text=stringr::str_replace_all(input_text,pattern = '\\|',replacement = '\\\\|')
    # 将用户输入分割为行，并去除前后空格
    search_terms <- str_trim(unlist(str_split(input_text, "\\n")))
    
    # 初始化结果数据框
    results <- data.frame()
    
    # 遍历每个搜索条目并查找匹配
    for (term in search_terms) {
      # 确保搜索条目不为空
      if (nzchar(term)) {
        # 使用 str_detect 进行模糊匹配
        matched_rows <- geneinfo[
          str_detect(tolower(geneinfo$GeneID), regex(term, ignore_case = TRUE)) | 
            str_detect(tolower(geneinfo$mRNA.ID), regex(term, ignore_case = TRUE)) |
            str_detect(tolower(geneinfo$description), regex(term, ignore_case = TRUE)) , ] 
          
        
        # 将匹配的行添加到结果
        results <- rbind(results, matched_rows)
      }
    }
    
    # 输出结果
    output$search_tbl <- renderDT({
      datatable(
        unique(results),
        extensions = 'Buttons',
        filter = list(
          position = 'top',
          clear = TRUE,
          plain = FALSE
        ),
        options = list(
          dom = 'Bfrtip',
          searchHighlight = TRUE,
          buttons = c('copy', 'csv', 'excel', 'print', 'pdf')
        )
      )  # 确保结果不重复
    })  
  })
  

  
#############2.Id mapping########### 
  observeEvent(input$submit1, {
    id<- trimws(input$gene_id1)
    if (nchar(id) == 0) {
      output$error_message2 <- renderText("Please enter a gene ID.")
      return()  # 提前返回
    } else if (!(id %in% RSAG1.stat$mRNA.ID | id %in% RSAG1.stat$GeneID)) {
      output$error_message2 <- renderText("No matching Gene ID found.")
      return()  # 提前返回
    }

    # 确定 geneid 和 mRNAid
    if (id %in% RSAG1.stat$GeneID) {
      mRNAid <- RSAG1.stat[RSAG1.stat$GeneID == id, 'mRNA.ID']

    } else {
      geneid <- RSAG1.stat[RSAG1.stat$mRNA.ID == id, 'GeneID']
      mRNAid <- RSAG1.stat[RSAG1.stat$GeneID == geneid, 'mRNA.ID']

    }

    
    strain_id <-  RSAG1.stat[RSAG1.stat$mRNA.ID%in%mRNAid, 'strain']
    
    strain_id <- strain_id[1]

    mRNAid=stringr::str_replace_all(mRNAid,pattern = '\\|',replacement = '\\\\|')
    pattern <- paste(mRNAid, collapse = "|")
    
    Orthogroups.out =  Orthogroups.tbl|>dplyr::filter(stringr::str_detect(.data[[strain_id]],pattern))

    output$geneTable <- renderDT({

      genetable=Orthogroups.out[1,]|>t()|>as.data.frame()
      genetable=data.frame(Strain=rownames(genetable)[-1],Orthogroup.transcript=genetable$V1[-1])
      datatable(genetable)
    })
    
  })

  
  
  
#############3.提取序列########### 
   ###载入所有菌株基因号/转录本号
 
  observeEvent(input$submit2, {
    
    # 在执行比对前显示模态对话框
    showModal(modalDialog(
      title = "Message",
      "Extracting sequence, please wait...",
      easyClose = FALSE,  # 不允许用户关闭对话框
      footer = NULL
    ))
    
    strain_id <- trimws(input$strain2)
    id <- trimws(input$gene_id2)
    shinyjs::disable("submit2")
    # 载入对应菌株的序列和注释

    tryCatch({
      strain.fa <- Biostrings::readDNAStringSet(paste0("www/strainseq/genome/", 
                                                       strain_id, ".genomic.fa"))
    }, error = function(e) {
      output$sequence_out1 <- renderText(paste("Error reading fasta file:", e$message))
    })
    

    strain.mRNA.fa <- Biostrings::readDNAStringSet(paste0('www/strainseq/mRNA/', strain_id, '.mRNA.fa'))
    strain.pep.fa <- Biostrings::readAAStringSet(paste0('www/strainseq/peptide/', strain_id, '.pep.fa'))
    
    
    strain.gff <- tryCatch(
      {
        read.table(paste0('www/strainseq/gff3/', strain_id, '.gff3'), sep = '\t', quote = '')
      },
      error = function(e) {
        output$sequence_out1 <- renderText("Error reading GFF file.")
        return(NULL)
      }
    )
    
    strain.stat <- read.table(paste0('www/strainseq/stat/', strain_id, '.stat'), sep = '\t', header = TRUE, quote = '')
    
    strain.strand <- dplyr::filter(strain.gff, V3 == 'mRNA')
    strain.strand <- dplyr::transmute(strain.strand, strand = V7, mRNA.ID = stringr::str_extract(V9, '(?<=ID=)[^;]+'))
    strain.stat <- dplyr::left_join(strain.stat, strain.strand, by = 'mRNA.ID')
    
    # 检查是否输入为空字符串
    if (nchar(id) == 0) {
      removeModal()
      output$sequence_out1 <- renderText("Please select a Strain and enter a gene ID.")
      return()  # 提前返回
    } else if (!(id %in% strain.stat$mRNA.ID | id %in% strain.stat$GeneID)) {
      removeModal()
      output$sequence_out1 <- renderText("Please check the Strain or Gene ID.")
      output$sequence_out2 <- renderText(" ")
      output$sequence_out3 <- renderText(" ")
      return()  # 提前返回
    }
    
    # 确定 geneid
    if (id %in% strain.stat$GeneID) {
      geneid <- id
    } else {
      geneid <- strain.stat[strain.stat$mRNA.ID == id, 'GeneID']
    }
    
    # 检查 geneid 是否有效
    if (length(geneid) == 0) {
      removeModal()
      output$sequence_out1 <- renderText("No Gene ID found for the given mRNA ID.")
      return()
    }
    
    ### 提取 gene 序列 ###
    gene_chr <- strain.stat[strain.stat$GeneID == geneid, 'Chr.ID'][1]
    gene_range <- strain.stat[strain.stat$GeneID == geneid, 'Gene.Range'][1]
    gene_strand <- strain.stat[strain.stat$GeneID == geneid, 'strand'][1]
    #染色体编号
    chr_num <-  which(stringr::str_detect(names(strain.fa), gene_chr))
    # 使用自定义提取函数
    gene_sequence <- extract_and_merge_sequences(fa = strain.fa, chr = chr_num, range = gene_range, strand = gene_strand)
    
    gene_sequence <- wrap_text(gene_sequence, width = 120)

    ### 提取 gene 上下游序列 ###
    chrlen=length(strain.fa[[chr_num]]) 
    # 使用自定义提取坐标

    upstream_range = get_upstream_coordinates(gene_range, input$upstream, gene_strand, chrlen)

    # 提取序列
    gene_up_sequence <- extract_and_merge_sequences(fa = strain.fa, chr = chr_num, range = upstream_range, strand = gene_strand)
    gene_up_sequence <- wrap_text(gene_up_sequence, width = 120)
    
    ### 下游    
    downstream_range = get_downstream_coordinates(gene_range,input$downstream, gene_strand, chrlen)
    
    # 提取序列
    gene_down_sequence <- extract_and_merge_sequences(fa = strain.fa, chr = chr_num, range = downstream_range, strand = gene_strand)
    
    gene_down_sequence <- wrap_text(gene_down_sequence, width = 120)
    
    sequence_out1 <- paste0('>', geneid, ' gene\n', gene_sequence,
                            '\n',
                           '>', geneid, ' up_stream\n', gene_up_sequence,
                           '\n',
                           '>', geneid, ' down_stream\n', gene_down_sequence)
      
    ### 提取 mRNA 序列 ###
    
    mRNAid=strain.stat[strain.stat$GeneID == geneid, 'mRNA.ID']
    
    mRNAid_num=extract_matching_strings(names(strain.mRNA.fa), mRNAid)

    mRNA_sequence=strain.mRNA.fa[mRNAid_num]
    
    mRNA_sequence=as.character(mRNA_sequence) 
    
    mRNA_sequence <- convert_to_fasta(mRNA_sequence)
    
    sequence_out2 <- mRNA_sequence
      


    ### 提取 CDS 序列 ###
    
    mRNAid=strain.stat[strain.stat$GeneID == geneid, 'mRNA.ID']
    
    mRNAid_num=extract_matching_strings(names(strain.pep.fa), mRNAid)
    
    pep_sequence=strain.pep.fa[mRNAid_num]
    
    pep_sequence=as.character(pep_sequence) 
    
    pep_sequence <- convert_to_fasta(pep_sequence)
    
    sequence_out3 <- pep_sequence
    

    
    # 检查提取的基因序列是否为空
    if (nchar(sequence_out1) == 0) {
      output$sequence_out1 <- renderText("No sequence found for the given Gene ID.")
    } else {
      output$sequence_out1 <- renderText(sequence_out1)
    }
    
    if (nchar(sequence_out2) == 0) {
      output$sequence_out2 <- renderText("No sequence found for the given Gene ID.")
    } else {
      output$sequence_out2 <- renderText(sequence_out2)
    }
    
    if (nchar(sequence_out3) == 0) {
      output$sequence_out3 <- renderText("No sequence found for the given Gene ID.")
    } else {
      output$sequence_out3 <- renderText(sequence_out3)
    }
    
    # 删除大对象
    rm(strain.fa,strain.gff,strain.stat)
    
    # 运行垃圾回收
    gc()
    shinyjs::enable("submit2")
    removeModal()

  })
  
  
#############4.表达查询###########
  observeEvent(input$example_button2, {
    # 当按钮被点击时，填充内容
    updateTextAreaInput(session, "searchid2", value = "chr01g00790\nchr01g00790.mRNA1\nchr01g00790.mRNA2")
  })
  observeEvent(input$submit3, {
    
    #载入表达量数据
    
    gene_exp=read.table('www/combined_gene_tpm.tsv',header = T,sep = '\t')
    mRNA_exp=read.table('www/combined_transcript_tpm.tsv',header = T,sep = '\t')
    
    id<- trimws(input$searchid2)
    
    if (nchar(id) == 0) {
      output$error_message3 <- renderText("Please enter a gene ID or description.")
      return()  # 提前返回
    }

    # 将用户输入分割为行，并去除前后空格
    search_terms <- str_trim(unlist(str_split(id, "\\n")))
    
    # 含有 "mRNA" 的字符串
    mRNA_id <- search_terms[str_detect(search_terms,'mRNA')]
    
    # 不含 "mRNA" 的字符串
    gene_id <- search_terms[!str_detect(search_terms,'mRNA')]
    

    out_gene_tbl=gene_exp|>dplyr::filter(geneid%in%gene_id)
   
    
    if(nrow(out_gene_tbl)==0){
      output$gene_exp_table <- renderDT({
        
        datatable(data.frame(warning='No gene ID found'),options = list(dom = 't'),caption = htmltools::tags$caption(
          htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;', 'Table of Gene Expression Levels'))
        )
        
      })
    }else{
          output$gene_exp_table <- renderDT({
            
            datatable(out_gene_tbl,caption = htmltools::tags$caption(
              
              htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;', 'Table of Gene Expression Levels'),
              htmltools::p(style = 'font-size:16px; text-align: left; color: black;', '"r" represents the mycelial sample of strain YZ19, 
                       while "rN" represents the infected sample taken from the leaf sheath of the rice variety Nipponbare after inoculation.')
              
            )
            )
      
    })
    }
    
    
    
    
    out_mRNA_tbl=mRNA_exp|>dplyr::filter(mRNAid%in%mRNA_id)
    
    
    if(nrow(out_mRNA_tbl)==0){
      output$transcript_exp_table <- renderDT({
        
        datatable(data.frame(warning='No mRNA ID found'),caption = htmltools::tags$caption(
          htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;', 'Table of Transcript Expression Levels'))
        )
      })
    }else{
      output$transcript_exp_table <- renderDT({
        
        datatable(out_mRNA_tbl,caption = htmltools::tags$caption(

          htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;', 'Table of Transcript Expression Levels'),
          htmltools::p(style = 'font-size:16px; text-align: left; color: black;', '"r" represents the mycelial sample of strain YZ19, 
                       while "rN" represents the infected sample taken from the leaf sheath of the rice variety Nipponbare after inoculation.')
          
          )
        )
      })
    }
    
    output$error_message3 <- NULL
    
  })
#############5.下载########
  
  # 下载文件的处理
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  output$download1.1 <- downloadHandler(
    filename = function() {
      "YZ19.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/YZ19.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/YZ19.genomic.fa.gz', file)
    }
  )
  
  output$download1.2 <- downloadHandler(
    filename = function() {
      "YZ19.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/YZ19.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/YZ19.gff3.gz', file)
    }
  )
  
  output$download1.3 <- downloadHandler(
    filename = function() {
      "YZ19.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/YZ19.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/YZ19.mRNA.fa.gz', file)
    }
  )
  
  output$download1.4 <- downloadHandler(
    filename = function() {
      "YZ19.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/YZ19.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/YZ19.pep.fa.gz', file)
    }
  )
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  output$download2.1 <- downloadHandler(
    filename = function() {
      "XN.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/XN.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/XN.genomic.fa.gz', file)
    }
  )
  
  output$download2.2 <- downloadHandler(
    filename = function() {
      "XN.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/XN.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/XN.gff3.gz', file)
    }
  )
  
  output$download2.3 <- downloadHandler(
    filename = function() {
      "XN.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/XN.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/XN.mRNA.fa.gz', file)
    }
  )
  
  output$download2.4 <- downloadHandler(
    filename = function() {
      "XN.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/XN.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/XN.pep.fa.gz', file)
    }
  )
 
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download3.1 <- downloadHandler(
    filename = function() {
      "HG81.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/HG81.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/HG81.genomic.fa.gz', file)
    }
  )
  
  output$download3.2 <- downloadHandler(
    filename = function() {
      "HG81.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/HG81.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/HG81.gff3.gz', file)
    }
  )
  
  output$download3.3 <- downloadHandler(
    filename = function() {
      "HG81.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/HG81.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/HG81.mRNA.fa.gz', file)
    }
  )
  
  output$download3.4 <- downloadHandler(
    filename = function() {
      "HG81.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/HG81.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/HG81.pep.fa.gz', file)
    }
  )

  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download4.1 <- downloadHandler(
    filename = function() {
      "BRS1.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/BRS1.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/BRS1.genomic.fa.gz', file)
    }
  )
  
  output$download4.2 <- downloadHandler(
    filename = function() {
      "BRS1.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/BRS1.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/BRS1.gff3.gz', file)
    }
  )
  
  output$download4.3 <- downloadHandler(
    filename = function() {
      "BRS1.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/BRS1.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/BRS1.mRNA.fa.gz', file)
    }
  )
  
  output$download4.4 <- downloadHandler(
    filename = function() {
      "BRS1.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/BRS1.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/BRS1.pep.fa.gz', file)
    }
  )
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download5.1 <- downloadHandler(
    filename = function() {
      "B2.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/B2.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/B2.genomic.fa.gz', file)
    }
  )
  
  output$download5.2 <- downloadHandler(
    filename = function() {
      "B2.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/B2.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/B2.gff3.gz', file)
    }
  )
  
  output$download5.3 <- downloadHandler(
    filename = function() {
      "B2.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/B2.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/B2.mRNA.fa.gz', file)
    }
  )
  
  output$download5.4 <- downloadHandler(
    filename = function() {
      "B2.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/B2.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/B2.pep.fa.gz', file)
    }
  )
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download6.1 <- downloadHandler(
    filename = function() {
      "ADB.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/ADB.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/ADB.genomic.fa.gz', file)
    }
  )
  
  output$download6.2 <- downloadHandler(
    filename = function() {
      "ADB.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/ADB.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/ADB.gff3.gz', file)
    }
  )
  
  output$download6.3 <- downloadHandler(
    filename = function() {
      "ADB.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/ADB.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/ADB.mRNA.fa.gz', file)
    }
  )
  
  output$download6.4 <- downloadHandler(
    filename = function() {
      "ADB.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/ADB.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/ADB.pep.fa.gz', file)
    }
  )
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download7.1 <- downloadHandler(
    filename = function() {
      "WGL.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/WGL.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/WGL.genomic.fa.gz', file)
    }
  )
  
  output$download7.2 <- downloadHandler(
    filename = function() {
      "WGL.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/WGL.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/WGL.gff3.gz', file)
    }
  )
  
  output$download7.3 <- downloadHandler(
    filename = function() {
      "WGL.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/WGL.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/WGL.mRNA.fa.gz', file)
    }
  )
  
  output$download7.4 <- downloadHandler(
    filename = function() {
      "WGL.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/WGL.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/WGL.pep.fa.gz', file)
    }
  )
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download8.1 <- downloadHandler(
    filename = function() {
      "YN7.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/YN7.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/YN7.genomic.fa.gz', file)
    }
  )
  
  output$download8.2 <- downloadHandler(
    filename = function() {
      "YN7.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/YN7.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/YN7.gff3.gz', file)
    }
  )
  
  output$download8.3 <- downloadHandler(
    filename = function() {
      "YN7.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/YN7.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/YN7.mRNA.fa.gz', file)
    }
  )
  
  output$download8.4 <- downloadHandler(
    filename = function() {
      "YN7.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/YN7.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/YN7.pep.fa.gz', file)
    }
  )
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download9.1 <- downloadHandler(
    filename = function() {
      "GD118.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/GD118.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/GD118.genomic.fa.gz', file)
    }
  )
  
  output$download9.2 <- downloadHandler(
    filename = function() {
      "GD118.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/GD118.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/GD118.gff3.gz', file)
    }
  )
  
  output$download9.3 <- downloadHandler(
    filename = function() {
      "GD118.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/GD118.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/GD118.mRNA.fa.gz', file)
    }
  )
  
  output$download9.4 <- downloadHandler(
    filename = function() {
      "GD118.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/GD118.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/GD118.pep.fa.gz', file)
    }
  )
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download10.1 <- downloadHandler(
    filename = function() {
      "BM1.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/BM1.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/BM1.genomic.fa.gz', file)
    }
  )
  
  output$download10.2 <- downloadHandler(
    filename = function() {
      "BM1.gff3.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/gff3/BM1.gff3.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/gff3/BM1.gff3.gz', file)
    }
  )
  
  output$download10.3 <- downloadHandler(
    filename = function() {
      "BM1.mRNA.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/mRNA/BM1.mRNA.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/mRNA/BM1.mRNA.fa.gz', file)
    }
  )
  
  output$download10.4 <- downloadHandler(
    filename = function() {
      "BM1.pep.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/peptide/BM1.pep.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/peptide/BM1.pep.fa.gz', file)
    }
  )
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  output$download11.1 <- downloadHandler(
    filename = function() {
      "1082KB.genome.fa.gz"  # 客户端下载时的文件名
    },
    content = function(file) {
      # 确保路径正确，并检查文件是否存在
      if (!file.exists('www/strainseq/genome/1082KB.genomic.fa.gz')) {
        stop("File does not exist!")
      }
      # 将服务器上的文件拷贝到下载文件
      file.copy('www/strainseq/genome/1082KB.genomic.fa.gz', file)
    }
  )
  
#############6.序列比对##########
  # 初始化反应值
  reactiveData_align <- reactiveVal(NULL)
  observeEvent(input$align, {
    # 在执行比对前显示模态对话框
    showModal(modalDialog(
      title = "提示",
      "正在比对，请稍候...",
      easyClose = FALSE,  # 不允许用户关闭对话框
      footer = NULL
    ))
    
    seq1 = check_sequence_info(input$align_seq1, sequence_name = 'seq1')
    seq2 = check_sequence_info(input$align_seq2, sequence_name = 'seq2')
    
    alignmentType=input$alignmentType
    gap_open=input$gap_open
    gap_extension=input$gap_open
    substitutionMatrix=input$substitutionMatrix
    
    if (seq1[1] == 0 | seq2[1] == 0) {
      output$alignmentResult <- renderText("Incorrect input!")
      removeModal()  # 关闭模态对话框
      return()  # 提前返回
    }
    
    if (seq1[1] != seq2[1]) {
      output$alignmentResult <- renderText("The sequence types are inconsistent.")
      removeModal()  # 关闭模态对话框
      return()  # 提前返回
    }
    
    # 将输入序列转换为Biostrings对象
    if (seq1[1] == 'nt') {
      p1 = DNAStringSet(seq1[3])
      names(p1) = seq1[2]
      p2 = DNAStringSet(seq2[3])
      names(p2) = seq2[2]
    } else {
      p1 = AAStringSet(seq1[3])
      names(p1) = seq1[2]
      p2 = AAStringSet(seq2[3])
      names(p2) = seq2[2]
    }
    
    # 执行全局比对或局部比对
    if (seq1[1] == 'nt') {
      alignment <- pwalign::pairwiseAlignment(p1, p2, type = alignmentType,
                                              gapOpening = gap_open,
                                              gapExtension = gap_extension)
    } else {
      alignment <- pwalign::pairwiseAlignment(p1, p2, type = alignmentType,
                                              gapOpening = gap_open,
                                              gapExtension = gap_extension,
                                              substitutionMatrix=substitutionMatrix)
    }
    
    # 输出比对结果
    output$alignmentResult <- renderPrint({
      pwalign::writePairwiseAlignments(alignment, block.width = input$blockwidth)
    })
    
      # 比对结果存储在反应值中
      alignment_type=seq1[1]
 
      p1=as.character(alignment@pattern)
      names(p1)=alignment@pattern@unaligned@ranges@NAMES
      
      p2=as.character(alignment@subject)
      names(p2)=alignment@subject@unaligned@ranges@NAMES
      alignment_max=alignment@subject@range@width
      alignment_fa_out=c(p1,p2)   
      
      alignment_plot_data=list(alignment_type,alignment_fa_out)
      
      reactiveData_align(alignment_plot_data)    
      
      #生成绘图范围
      output$after_align_range <- renderUI({
        sliderInput('align_plot_range','Range of alignment visualization',
                    min = 1, max = alignment_max,step = 1,
                    value = c(1, ceiling(alignment_max/10) ))
      })
      
      #生成绘图按钮
      output$after_align_button <- renderUI({
        actionButton('align_plot_button','Visualization',
                     icon = icon('compass-drafting'))
      })
      

    removeModal()  # 比对完成后关闭模态对话框
  })
  
  
  
  
  
  observeEvent(input$align_plot_button, {

     plotData =  reactiveData_align()

         if(plotData[[1]] == 'nt'){
       output$alignmentPlot <- renderPlot({
       ggmsa( DNAMultipleAlignment(plotData[[2]]), 
              start = input$align_plot_range[1],
              end = input$align_plot_range[2],
              color = "Shapely_NT",
                 char_width = 0.7, seq_name = TRUE)
                                         })
     }else{output$alignmentPlot <- renderPlot({
         ggmsa( AAMultipleAlignment(plotData[[2]]), 
                start = input$align_plot_range[1],
                end = input$align_plot_range[2],
                color = "Clustal",
                char_width = 0.7, seq_name = TRUE)
                                              })
     } 
          
     output$after_align_plot <- renderUI({
       req(reactiveData_align())
       downloadButton("downloadPlot", "Print",icon = icon('print'))
     })
     
          
  })
  
  output$downloadPlot <- downloadHandler(

    filename = function() {
      paste("plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      plotData =  reactiveData_align()

      if(plotData[[1]] == 'nt'){
         ggfigure = ggmsa( DNAMultipleAlignment(plotData[[2]]),
                 start = input$align_plot_range[1],
                 end = input$align_plot_range[2],
                 color = "Shapely_NT",
                 char_width = 0.7, seq_name = TRUE)
      }else{
        ggfigure= ggmsa( AAMultipleAlignment(plotData[[2]]),
               start = input$align_plot_range[1],
               end = input$align_plot_range[2],
               color = "Clustal",
               char_width = 0.7, seq_name = TRUE)
      }
      ggsave(file, plot = ggfigure, device = "png",dpi = 300)
    }

  )
      reactiveData_align <- reactiveVal(NULL)
#############7.orf预测########  
  
  observeEvent(input$getORF, {
    req(input$seq_need_ORF)
    seq=input$seq_need_ORF
    min_pep_len=input$orf_len
    include_reverse = input$seq_need_ORF_strand
    if (grepl("^>", seq)) {
      sequence_lines <- unlist(strsplit(seq, "\n"))
      sequence_name <- gsub("^>", "", sequence_lines[1])
      actual_sequence <- paste(sequence_lines[-1], collapse = "")
    } else {
      actual_sequence <- seq
    }
    
    cleaned_seq <- remove_special_characters(actual_sequence)
    sequence_type <- NA
    if (grepl("^[ACGT]+$", cleaned_seq, ignore.case = TRUE)) {
      sequence_type <- "nt"
    } else {
      output$ORF_seq <- renderText("incorrect input sequence")
      return()
    }
    
    cleaned_seq=DNAString(cleaned_seq)
    
    
    out_orf=find_longest_orf(cleaned_seq,include_reverse = include_reverse)
    
    out_orf_seq=out_orf[[1]]
    out_orf_start=out_orf[[2]]
    out_orf_end=out_orf[[3]]

    
    if(nchar(out_orf_seq)<min_pep_len*3){
          
    output$ORF_seq <- renderText('none ORF find!')
    output$ORF_pep <-renderText(" none peptide find!")
    
    return()  
      
    }else{
      
    orf_seq= as.character(out_orf_seq)
    final_orf_out=paste0('>ORF start: ',out_orf_start,' end: ', 
                     out_orf_end+3,' length: ',nchar(out_orf_seq),'\n',
                     wrap_text(orf_seq,120))
    output$ORF_seq <- renderText( final_orf_out)
    out_pep=Biostrings::translate(out_orf_seq)
    pep_seq=as.character(out_pep)
    final_pep_out=paste0(">aa\n",wrap_text(pep_seq,120))
    output$ORF_pep <- renderText(final_pep_out)
    }

  })
  
  
########### 8.gffrename############
  # Define reactive expression to get the uploaded file path
  file_data <- reactive({
    req(input$upload_gff)  # 确保文件已经上传
    infile <- input$upload_gff$datapath   # 获取上传文件的路径
    return(infile)
  })
  reactiveData_chr <- reactiveVal(NULL)
  output$fileContents <- renderText({
    req(input$upload_gff)  # 确保文件已经上传
    paste("file name:", input$upload_gff$name, 
          "\nfile size:", input$upload_gff$size, "\n")
  })
  
  observeEvent(input$upload_gff, {
    gff_checkout <- check_gff3_format(file_data())
    if (gff_checkout[1] != 0) {
      error_message <- paste0('There is an error on line ',
                              gff_checkout[1],
                              ' , looks like incorrect ',
                              gff_checkout[2],
                              ' , please check.')
      showModal(modalDialog(
        title = "Oops!",
        error_message,
        easyClose = TRUE,  
        footer = NULL
      ))
      output$after_upload_gff <- renderUI({ NULL })
      return(NULL)  # 停止执行后续代码
    } else {
      old_gffdata <- read.table(file_data(), header = FALSE, sep = '\t', quote = '', stringsAsFactors = FALSE)
      
      output$headContents <- renderDT({
        colnames(old_gffdata) <- c('seqid','source','type','start','end','score','strand', 'phase','attributes')
        datatable(old_gffdata[1:100,], caption = htmltools::tags$caption(
          htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;',
                       'top 100 rows of input gff3 file'))
        )
      })
      
      output$orgin_chr <- renderText({
        chr = unique(old_gffdata$V1)
        paste(chr, collapse = '\n')
      })
      
      output$old_gff_info <- renderText({
        info = extract_info_from_gff(old_gffdata)
        paste0('Counts \n',
              'chr nums: ',info[1],'\n',
              'gene nums: ',info[2],'\n',
              'transcript nums: ',info[3])
      })
        output$after_upload_gff_newchr <- renderUI({
        chr_num=length(unique(old_gffdata$V1)) 
        digits = nchar(chr_num)
        defaultchr=str_pad(1:chr_num,width = digits,side = 'left',pad = 0)
        reactiveData_chr(
          paste0('chr',defaultchr)
          )         
        defaultchr=paste0('chr',defaultchr,collapse = '\n')
        

        placeholder=paste0( "Enter new ID to replace default, if you wish.\n",
                            defaultchr)
        textAreaInput("newchr", label = "", 
                      rows = 50, placeholder =placeholder)
      })
      output$after_upload_gff_actionbutton <- renderUI({
        actionButton("process_file", "Rename gff3 file")
      })
      output$after_upload_gff_p1 <- renderUI({
        h4("Old seqid of input annotation: ")
      })
      output$after_upload_gff_p2 <- renderUI({
        h4("new seqid(Match old IDs one-to-one.): ")
      })
      
      output$new_gff_id_digits <- renderUI({
       sliderInput('id_digits','Digit count of numerical ID string.',
                    min = 5, max = 10, value = 6,step = 1)
      })

      output$new_gff_id_step <- renderUI({
        radioButtons('id_step','Increment step of numerical ID string.',
                     choices = list("1" = 1,
                                    "5" = 5,
                                    "10" = 10,
                                    "50" = 50), 
                     selected = 10)
      })
      
      output$new_id_sample <- renderText({
        id_digits=input$id_digits
        id_step=as.numeric(input$id_step) 
        id=seq(from = 0, by = id_step, length.out = 5)
        id=str_pad(id,width = id_digits,side = 'left',pad = 0)
        id= paste0('chr01g',id,collapse = '\n')
        return(paste0('examples: \n',id))
      })
      
        output$overlap <- renderUI({
        sliderInput('mrna_overlap','overlap relative to the total length of the transcript: ',
                    min = 0.1, max = 1, value = 0.5,step = 0.1)
      })

    }
  })
  ########
  observeEvent(input$process_file, {
    
    showModal(modalDialog(
      title = "Message",
      "Processing, please wait...",
      easyClose = FALSE,  # 不允许用户关闭对话框
      footer = NULL
    ))  
    
    gff_file = file_data()
    
    ## 如果没有输入新id，则使用默认id
    if (nchar(input$newchr) == 0) {
      newchr = reactiveData_chr()
    } else {
      input_text = input$newchr
      newchr = unlist(strsplit(input_text, "\\n"))
      newchr = newchr[nzchar(newchr)]
      # 检查条件
      if (length(newchr) != length(reactiveData_chr()) || 
          length(unique(newchr)) != length(newchr) ||  # 检查是否有重复
          !all(grepl("^[0-9_]+$", newchr))) {  # 检查是否只含有数字和下划线
        
        showModal(modalDialog(
          title = "seqid error!",
          "Please ensure that the input meets the following conditions:\n
           Contains matched entries\n
           All entries are unique\n
           Contains only numeric characters and underscores",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
    } 
    
      new_digits = as.numeric(input$id_digits) 
      new_step =as.numeric(input$id_step) 
      new_overlap =as.numeric(input$mrna_overlap)

      result <- tryCatch({
        
      outgff =  fix_gff_annotation(gff_file,
                           digits = new_digits,
                           step = new_step,
                           chr = newchr,
                           overlap=new_overlap)
      

      output$after_fixgff_title <- renderUI({
        h4("gff3 rename sucessed!: ")
      })

      output$newgfftbl <- renderDT({
        colnames(outgff) <- c('seqid','source','type','start','end','score','strand', 'phase','attributes')
        datatable(outgff[1:100,], extensions = 'Buttons',
                  filter = list(
                    position = 'top',
                    clear = TRUE,
                    plain = FALSE
                  ),
                  options = list(
                    dom = 'Bfrtip',
                    searchHighlight = TRUE,
                    buttons = list('copy', 'print', 
                                   list(
                                        extend = 'collection',
                                        buttons = c('csv', 'excel'),
                                        text = 'Download')
                                   )
                  ),
                  caption = htmltools::tags$caption(
          htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;',
                       'top 100 rows of output gff3 file'))
        )
      })
        
      
      
      # 获取当前时间并格式化为字符串
      timestamp <- format(Sys.time(), "%Y%m%d")  # 例如：20231009_153045
      
      # 指定压缩文件名，添加时间戳
      file_name <- paste0("output/output.", timestamp, ".gff3.gz")
      
      # 保存数据框并压缩为 gz 格式
      write.table(outgff, gzfile(file_name), 
                  row.names = F,quote = F,sep = '\t',col.names = F)
      
      
      
      output$fixgff_download <- renderUI({
        downloadButton('downloadgff','Download full gff3 file',icon = icon('download'))
      })
      
      
      }, error = function(e) {
        return(FALSE)  # 捕获错误并返回自定义信息
      })
      
      if (is.logical(result)) {   
        showModal(modalDialog(
          title = "Message",
          "error, please check your gff3 file!",
          easyClose = TRUE,  # 不允许用户关闭对话框
          footer = NULL
        ))  
      }else{
        showModal(modalDialog(
          title = "Message",
          "ok",
          easyClose = TRUE,  # 不允许用户关闭对话框
          footer = NULL
        ))  
      }
    })
  
  output$downloadgff <- downloadHandler(
    filename = function() {
    paste0("Renamed.gff3.gz")
    },
    content = function(file) {
      file.copy( paste0("output/output.", 
                           format(Sys.time(), "%Y%m%d"),
                           ".gff3.gz"), file)
    }
  )
 
##############microexon######
  observeEvent(input$getme, {
    req(input$upload_BAM_file)
    req(input$upload_fa_file)
    infile_bam <- input$upload_BAM_file$datapath  
    infile_fa <- input$upload_fa_file$datapath  
    param <- ScanBamParam(what=scanBamWhat())
    result <- tryCatch({

    bam <- scanBam(infile_bam, param=param)
    ref.fa = readDNAStringSet(infile_fa)
    
    

    metbl=identity_micro_exon(bam=bam,fa=ref.fa)

        output$microexon_out <- renderDT({
        colnames(metbl) <- c('ReadsID ','SeqID','Strand','start','qwidth',
                             'cigar','Reads_seq', 'end','ref_seq')
        datatable(metbl, extensions = 'Buttons',
                  filter = list(
                    position = 'top',
                    clear = TRUE
                  ),
                  options = list(
                    dom = 'Bfrtipl',
                    searchHighlight = TRUE,
                    buttons = list('copy', 'print',
                                   list(
                                     extend = 'collection',
                                     buttons = c('csv', 'excel'),
                                     text = 'Download')
                    )
                  ),
                  caption = htmltools::tags$caption(
                    htmltools::p(style = 'font-size:20px; text-align: left; color: black;font-weight: bold;',
                                 'It is recommended to use gsaman for gene structure correction.'))
        )
      })

     },error = function(e) {
     return(FALSE)  # 捕获错误并返回自定义信息
     })

    if (is.logical(result)) {
      showModal(modalDialog(
        title = "Message",
        "error, please check your BAM or fasta file!",
        easyClose = TRUE,  # 不允许用户关闭对话框
        footer = NULL
      ))
    }
    
    
    
  })
  
}
  


shinyApp(ui, server)





