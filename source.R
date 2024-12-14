wrap_text <- function(input_string, width = 40) {
  # 获取字符串的总长度
  input_length <- nchar(input_string)
  
  # 创建一个空的结果字符串
  result <- ""
  
  # 逐段提取每一段并添加换行符
  for (start in seq(1, input_length, width)) {
    end <- min(start + width - 1, input_length)  # 计算每段的结束位置
    result <- paste0(result, substring(input_string, start, end), "\n")
  }
  
  return(result)
}


# 定义函数提取并合并序列
extract_and_merge_sequences <- function(fa, chr, range, strand = '+') {
  # 解析区间字符串
  intervals <- strsplit(range, ",")[[1]]
  intervals <- lapply(intervals, function(x)
    as.numeric(unlist(strsplit(x, ":"))))
  

  # 提取片段
  extracted_sequences <- sapply(intervals, function(interval) {
    start <- interval[1]
    end <- interval[2]
    as.character(Biostrings::subseq(fa[[chr]], start, end))
  })
  
  # 拼接提取的片段
  
  new_sequence <- paste(extracted_sequences, collapse = "")
  
  if (strand == '-') {
    new_sequence <- DNAString(new_sequence)
    new_sequence <- Biostrings::reverseComplement(new_sequence)
    new_sequence <- as.character(new_sequence)
    
  }
  
  return(new_sequence)
}

# 定义函数返回基因上下游坐标
get_upstream_coordinates <- function(gene_coords, promoter_length, strand, chromosome_length) {
  # 将基因坐标字符串分割为起始和结束位置
  coords <- unlist(strsplit(gene_coords, ":"))
  start <- as.numeric(coords[1])  # 起始位置
  end <- as.numeric(coords[2])      # 结束位置
  
  if (strand == "+") {
    # 正链：启动子在基因的左侧
    promoter_start <- start - promoter_length
    promoter_end <- start - 1
    
    # 限制在染色体边界内
    promoter_start <- max(promoter_start, 1)  # 启动子最小坐标为 1
    promoter_end <- max(promoter_end, promoter_start)  # 确保结束坐标不小于起始坐标
    
  } else if (strand == "-") {
    # 负链：启动子在基因的右侧
    promoter_start <- end + 1
    promoter_end <- end + promoter_length
    
    # 限制在染色体边界内
    promoter_end <- min(promoter_end, chromosome_length)  # 启动子最大坐标不超过染色体长度
    promoter_start <- min(promoter_start, promoter_end)  # 确保起始坐标不大于结束坐标
    
  } else {
    stop("Invalid strand: must be '+' or '-'")
  }
  
  # 返回启动子的坐标
  return(paste(promoter_start, ":", promoter_end, sep = ""))
}


get_downstream_coordinates <- function(gene_coords, promoter_length, strand, chromosome_length) {
  # 将基因坐标字符串分割为起始和结束位置
  coords <- unlist(strsplit(gene_coords, ":"))
  start <- as.numeric(coords[1])  # 起始位置
  end <- as.numeric(coords[2])      # 结束位置
  
  if (strand == "+") {
    # 正链：启动子在基因的下游（结束位置之后）
    promoter_start <- end + 1
    promoter_end <- end + promoter_length
    
    # 限制在染色体边界内
    promoter_start <- max(promoter_start, 1)  # 启动子最小坐标为 1
    promoter_end <- min(promoter_end, chromosome_length)  # 启动子最大坐标不超过染色体长度
    
  } else if (strand == "-") {
    # 负链：启动子在基因的上游（起始位置之前）
    promoter_start <- start - promoter_length
    promoter_end <- start - 1
    
    # 限制在染色体边界内
    promoter_start <- max(promoter_start, 1)  # 启动子最小坐标为 1
    promoter_end <- max(promoter_end, promoter_start)  # 确保结束坐标不小于起始坐标
    
  } else {
    stop("Invalid strand: must be '+' or '-'")
  }
  
  # 返回启动子的坐标
  return(paste(promoter_start, ":", promoter_end, sep = ""))
}

extract_matching_strings <- function(string_vector, target_vector) {
  # 使用 paste0 将目标向量中的元素组合成正则表达式模式
  # 元素名内存在‘|’时，先将其转义
  target_vector=stringr::str_replace_all(
    target_vector, 
    '\\|','\\\\|')
  pattern <- paste(target_vector, collapse = "|")
  
  # 使用 grepl 查找匹配的字符串
  matching_strings_num <- which(
    stringr::str_detect(
      string_vector, 
      pattern)) 
  
  return(matching_strings_num)
}


convert_to_fasta <- function(dna_vector) {
  # 检查输入
  if (is.null(names(dna_vector))) {
    sequence_names <- paste("seq", seq_along(dna_vector), sep = "_")  
  }else{
    sequence_names <- names(dna_vector)# 默认序列名
  }
  
  if (length(dna_vector) != length(sequence_names)) {
    stop("Length of dna_vector and sequence_names must be the same.")
  }
  
  # 创建 FASTA 格式的字符串
  fasta_strings <- sapply(seq_along(dna_vector), function(i) {
    seq=wrap_text(dna_vector[i],120)
    paste0(">", sequence_names[i], "\n", seq)
  }, USE.NAMES = FALSE)
  
  # 合并为单个字符串
  fasta_output <- paste(fasta_strings, collapse = "\n")
  
  return(fasta_output)
}

# 从序列中清理特殊字符
remove_special_characters <- function(str) {
  cleaned_str <- gsub("[\n\t\r ]", "", str)  # 去除换行符、制表符、回车符和空格
  cleaned_str <- gsub("[^A-Za-z0-9]", "", cleaned_str)  # 去除非字母和数字的字符
  return(cleaned_str)
}

# 检查序列是nt还是aa，返回序列类型，序列名称，序列本体
check_sequence_info <- function(sequence_str, sequence_name = 'seq1') {
  if (grepl("^>", sequence_str)) {
    sequence_lines <- unlist(strsplit(sequence_str, "\n"))
    sequence_name <- gsub("^>", "", sequence_lines[1])
    actual_sequence <- paste(sequence_lines[-1], collapse = "")
  } else {
    actual_sequence <- sequence_str
  }
  
  cleaned_seq <- remove_special_characters(actual_sequence)
  sequence_type <- NA
  if (grepl("^[ACGT]+$", cleaned_seq, ignore.case = TRUE)) {
    sequence_type <- "nt"
  } else if (grepl("^[ACDEFGHIKLMNPQRSTUVWY]+$", cleaned_seq, ignore.case = TRUE)) {
    sequence_type <- "aa"
  } else {
    return(0)
  }
  
  return(c(sequence_type, sequence_name, cleaned_seq))
}

find_longest_orf <- function(dna_sequence, include_reverse = TRUE) {# 找到最长的 ORF
  
  n <- nchar(dna_sequence)  # 获取 DNA 序列的长度
  longest_orf <- ""
  
  find_orf <- function(sequence,seqlen) {  
    local_longest_orf <- ""
    start_codon <- "ATG"
    stop_codons <- c("TAA", "TAG", "TGA")
    orf_start=0
    orf_end=0
    start_codon=Biostrings::matchPattern(start_codon,sequence)

    start_loc=start_codon@ranges@start
    for (i in start_loc) {      # 找到起始密码子，查找终止密码子
      for (j in seq(i+3, seqlen - 2, by = 3)) {
        stop_codon <- substr(sequence, j, j + 2)
        stop_codon=as.character(stop_codon)
        if (stop_codon %in% stop_codons) {
          current_orf <- substr(sequence, i, j + 2)
          if (nchar(current_orf) > nchar(local_longest_orf)) {
            local_longest_orf <- current_orf  # 更新最长 ORF
            orf_start=i
            orf_end=j
          }
          break  # 找到终止密码子后退出内层循环
        }
      }
    }
    return(list(local_longest_orf,orf_start,orf_end))
  }
  
  longest_orf <- find_orf(dna_sequence,n)  # 查找正向 ORF

  if (include_reverse) {  # 如果选择包括反向互补链
    reverse_seq <- Biostrings::reverseComplement(dna_sequence)
    longest_orf <- find_orf(reverse_seq,n)
  }
  
  return(longest_orf)
}



check_gff3_format <- function(file_path) {
  # 读取文件
  con <- file(file_path, "r")
  lines <- readLines(con)
  close(con)
  
  # 初始化检查结果

  error_line <- 0
  error_type='none' 
  
  # 定义允许的类型
  allowed_types <- c("gene", "mRNA", "exon", "CDS",'region','tRNA','rRNA',
                     'intron','five_prime_UTR','three_prime_UTR','pseudogene')

  # 检查每一行
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # 跳过注释行
    if (grepl("^#", line)) {
      next
    }
    
    # 分割字段
    fields <- strsplit(line, "\t")[[1]]
    
    # 检查列数
    if (length(fields) != 9) {
      is_valid <- FALSE
      error_line <- i
      error_type='obs'
      break
    }
    
    # 检查第三列的类型
    type <- fields[3]
    if (!(type %in% allowed_types)) {
      is_valid <- FALSE
      error_line <- i
      error_type='type'
      break
    }
  }
  
  # 返回检查结果

    return(c(error_line,error_type))

  
}


judge_mrna2gene <- function(mrnaA, mrnaB,overlap=0.5) {#接受mrna坐标向量，返回布尔值
  if (max(mrnaA[1], mrnaB[1]) <= min(mrnaA[2], mrnaB[2]) &
      (abs(mrnaA[2] - mrnaB[1]) > mrnaA[3]*overlap | abs(mrnaA[2] - mrnaB[1]) > mrnaB[3]*overlap)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


fix_gff_annotation <- function(gff_file,digits=7,step=10,chr,overlap=0.5) {
  # 读取gff文件为数据框
  gff <-  read.table(gff_file, header = FALSE, sep = '\t', quote = '', stringsAsFactors = FALSE)
  id_tbl=data.frame(V1=unique(gff$V1),seqid=chr)
  gff=gff|>dplyr::left_join(id_tbl,by='V1')
  gff=gff[,-1]
  gff=gff|>dplyr::rename(V1=seqid)
  gff=gff|>dplyr::select(9,1:8)
  #替换seqid#
  
  # 提取旧id到新的一列mRNAid
  gff <- gff |>
    dplyr::mutate(mRNAid = case_when(
      V3 == 'mRNA' ~ str_extract(V9, '(?<=ID=)[^;]+'),
      V3 == 'exon' ~ str_extract(V9, '(?<=Parent=)[^;]+'),
      V3 == 'CDS' ~ str_extract(V9, '(?<=Parent=)[^;]+')
    ))
  
  ##### 按照正反链拆分转录本,生成转录本长度mrnalen #####
  
  # 提取mRNA
  gff.mrna <- gff |>
    dplyr::filter(V3 == 'mRNA')
  
  # 正链转录本处理
  gff.mrna.plus <- gff.mrna |>
    dplyr::filter(V7 == '+') |>
    dplyr::arrange(V1, V4) |>
    dplyr::mutate(
      gene = 1,
      gene_start = 1,
      gene_end = 1,
      mrnalen = V5 - V4 + 1
    )
  # 初始基因编号
  gff.mrna.plus[1, 11] <- 1
  # 初始基因位置
  gff.mrna.plus[1, 12] <- gff.mrna.plus[1, 4]
  gff.mrna.plus[1, 13] <- gff.mrna.plus[1, 5]
  
  # 如果判定同一基因， 在上一个转录本的基础上拓宽基因范围，否则刷新为新的转录本的范围
  gff.mrna.plus <- as.data.frame(gff.mrna.plus)
  for (i in 2:nrow(gff.mrna.plus)) {
    if (judge_mrna2gene(mrnaA=gff.mrna.plus[i-1, c(12, 13, 14)],
                        mrnaB=gff.mrna.plus[i, c(4, 5, 14)],
                        overlap=overlap
                        )) {
      gff.mrna.plus[i, 11] <- gff.mrna.plus[i-1, 11]
      gff.mrna.plus[i, 12] <- min(gff.mrna.plus[i-1, 12], gff.mrna.plus[i, 4])
      gff.mrna.plus[i, 13] <- max(gff.mrna.plus[i-1, 13], gff.mrna.plus[i, 5])
    } else {
      gff.mrna.plus[i, 11] <- gff.mrna.plus[i-1, 11] + 1
      gff.mrna.plus[i, 12] <- gff.mrna.plus[i, 4]
      gff.mrna.plus[i, 13] <- gff.mrna.plus[i, 5]
    }
  }
  
  # 反链转录本处理
  gff.mrna.minus <- gff.mrna |>
    dplyr::filter(V7 == '-') |>
    dplyr::arrange(V1, V4) |>
    dplyr::mutate(
      gene = 1,
      gene_start = 1,
      gene_end = 1,
      mrnalen = V5 - V4 + 1
    )
  # 初始基因编号
  gff.mrna.minus[1, 11] <- 1000001
  # 初始基因位置
  gff.mrna.minus[1, 12] <- gff.mrna.minus[1, 4]
  gff.mrna.minus[1, 13] <- gff.mrna.minus[1, 5]
  
  # 如果判定同一基因， 在上一个转录本的基础上拓宽基因范围，否则刷新为新的转录本的范围
  gff.mrna.minus <- as.data.frame(gff.mrna.minus)
  for (i in 2:nrow(gff.mrna.minus)) {
    if (judge_mrna2gene(mrnaA=gff.mrna.minus[i-1, c(12, 13, 14)], 
                        mrnaB=gff.mrna.minus[i, c(4, 5, 14)],
                        overlap=overlap
                        )) {
      gff.mrna.minus[i, 11] <- gff.mrna.minus[i-1, 11]
      gff.mrna.minus[i, 12] <- min(gff.mrna.minus[i-1, 12], gff.mrna.minus[i, 4])
      gff.mrna.minus[i, 13] <- max(gff.mrna.minus[i-1, 13], gff.mrna.minus[i, 5])
    } else {
      gff.mrna.minus[i, 11] <- gff.mrna.minus[i-1, 11] + 1
      gff.mrna.minus[i, 12] <- gff.mrna.minus[i, 4]
      gff.mrna.minus[i, 13] <- gff.mrna.minus[i, 5]
    }
  }
  
  # 重新合并正反链转录本，整合基因信息
  gff.mrna.merge <- rbind(gff.mrna.plus, gff.mrna.minus) |>
    as.data.frame()
  
  # 按基因位置排序
  gff.mrna.merge <- gff.mrna.merge |>
    dplyr::arrange(V1, gene_start)
  
  # 正反基因在一起重新编号
  gff.mrna.merge <- gff.mrna.merge |>
    dplyr::ungroup() |>
    dplyr::group_by(V1) |>
    dplyr::mutate(diff = c(0, diff(gene)), geneid = cumsum(diff != 0))
  


  gff.mrna.merge <- gff.mrna.merge |>
    dplyr::mutate(geneid = paste0(V1,
                           'g',
                           str_pad(geneid + step
                                   , digits, side = 'left', '0')
                           ))
  
  # 提取gene成表
  gff.gene <- gff.mrna.merge |>
    dplyr::ungroup() |>
    dplyr::group_by(geneid) |>
    dplyr::reframe(
      V1 = V1,
      V2 = V2,
      V3 = 'gene',
      V4 = min(gene_start),
      V5 = max(gene_end),
      V6 = V6,
      V7 = V7,
      V8 = V8,
      V9 = paste0('ID=', geneid)
    ) |>
    unique() |>
    dplyr::select(-geneid) |>
    dplyr::ungroup()
  
  # 转录本编号按基因内位置顺序依次给号
  gff.mrna.merge <- gff.mrna.merge |>
    dplyr::arrange(V1, V3, V4) |>
    dplyr::group_by(geneid) |>
    dplyr::mutate(num = 1:n()) |>
    dplyr::mutate(m_id = paste0(geneid, '.t', num))
  
  # 提取mRNA成表,注意保留旧id用于更新exon表
  gff.mrna.rename <- gff.mrna.merge |>
    dplyr::ungroup() |>
    dplyr::mutate(V9 = paste0('ID=', m_id, ';Parent=', geneid)) |>
    dplyr::select(V1, V2, V3, V4, V5, V6, V7, V8, V9, mRNAid, newid = m_id)
  
  # 提取exon成表，更新parent信息

  gff.exon <- gff |>
    dplyr::filter(V3 == 'exon') |>
    dplyr::left_join(gff.mrna.rename[, 10:11], by = 'mRNAid') |>
    dplyr::mutate(V9 = paste0('Parent=', newid)) |>
    dplyr::select(1:9)
  
  # 提取exon成表，更新parent信息
  
  gff.cds <- gff |>
    dplyr::filter(V3 == 'CDS') |>
    dplyr::left_join(gff.mrna.rename[, 10:11], by = 'mRNAid') |>
    dplyr::mutate(V9 = paste0('Parent=', newid)) |>
    dplyr::select(1:9)
  
  
  # 合并基因，转录本和外显子表，导出
  gff.out <- rbind(gff.gene, gff.mrna.rename[, -10:-11], gff.exon,gff.cds) |>
    dplyr::arrange(V1, V4)
  
  return(gff.out)

}

extract_info_from_gff <- function(gff_file) {# 提取染色体数量，gene数量，转录本数量

 
  chr_num <- length(unique(gff_file$V1)) 
  gene_num <- nrow(gff_file[gff_file[,3]=='gene',]) 
  mrna_num <- nrow(gff_file[gff_file[,3]=='mRNA',]) 
  
  return(c(chr_num,gene_num,mrna_num))
}


add_CDS_to_gff=function(gff,transdecoder){

library(tidyverse)
##读取transdecoder.pep的id信息
cds.info = read.table('transdecoder.pep.info')
cds.info = cds.info |> transmute(
  mRNA = str_extract(V1, 'chr[0-9]{2}g[0-9]{5}.mRNA[0-9]+'),
  start = as.integer(str_extract(V7, '(?<=:)[0-9]+')) ,
  end = as.integer(str_extract(V7, '(?<=-)[0-9]+'))
)

#转录本内外显子排序，不排序会出bug，详见后续
gff.exon=arrange(gff.exon,mRNA,V4)

#########正链外显子与pep的位置关系,识别外显子是否包含CDS
exon_plus = gff.exon |>
filter(V7 == '+') |>
mutate(mRNA = str_extract(V9, 'chr[0-9]{2}g[0-9]{5}.mRNA[0-9]+')) |>
filter(!mRNA %in% rm.mRNA.id) |>
inner_join(cds.info, by = 'mRNA') |>
mutate(lens = V5 - V4 + 1) |>
group_by(mRNA) |>
mutate(sum_s = c(0, cumsum(lens[-length(lens)])) + 1,
         sum_e = cumsum(lens)) |>
mutate(
  trans_start = ifelse(start >= sum_s &
                         start <= sum_e, 'trans_start', ''),
  trans_end = ifelse(end >= sum_s &
                       end <= sum_e, 'trans_end', ''),
  cds = ifelse(start < sum_s & end > sum_e, 'cds', '')
) |>
mutate(type = paste0(trans_start, trans_end, cds))

###过滤掉UTR，计算CDS坐标
exon_plus = exon_plus |> ungroup() |> filter(type != '') |>
mutate(
  cds_start = case_when(
    type == 'cds' ~ V4,
    type == 'trans_end' ~ V4,
    type == 'trans_start' ~ as.integer(V4 + start - sum_s) ,
    type == 'trans_starttrans_end' ~ as.integer(V4 + start - sum_s)
  ),
  cds_end = case_when(
    type == 'cds' ~ V5,
    type == 'trans_end' ~ as.integer(V5 + end - sum_e),
    type == 'trans_start' ~ V5,
    type == 'trans_starttrans_end' ~ as.integer(V5 + end - sum_e)
  )
)

###############反链外显子同理，注意方向
exon_reverse = gff.exon |>
filter(V7 == '-') |>
mutate(mRNA = str_extract(V9, 'chr[0-9]{2}g[0-9]{5}.mRNA[0-9]+')) |>
filter(!mRNA %in% rm.mRNA.id) |>
inner_join(cds.info, by = 'mRNA') |>
mutate(lens = V5 - V4 + 1) |>
group_by(mRNA) |> arrange(desc(V4), .by_group = T) |>
mutate(sum_s = c(0, cumsum(lens[-length(lens)])) + 1,
         sum_e = cumsum(lens)) |>
mutate(
  trans_start = ifelse(start >= sum_s &
                         start <= sum_e, 'trans_start', ''),
  trans_end = ifelse(end >= sum_s &
                       end <= sum_e, 'trans_end', ''),
  cds = ifelse(start < sum_s & end > sum_e, 'cds', '')
) |>
mutate(type = paste0(trans_start, trans_end, cds))


exon_reverse = exon_reverse |> ungroup() |> filter(type != '') |>
mutate(
  cds_start = case_when(
    type == 'cds' ~ V4,
    type == 'trans_end' ~ as.integer(V4 - end + sum_e),
    type == 'trans_start' ~ V4 ,
    type == 'trans_starttrans_end' ~ as.integer(V4 - end + sum_e)
  ),
  cds_end = case_when(
    type == 'cds' ~ V5,
    type == 'trans_end' ~ V5,
    type == 'trans_start' ~ as.integer(V5 - start + sum_s),
    type == 'trans_starttrans_end' ~ as.integer(V5 - start + sum_s)
  )
)

###############产生cds位置信息并合并正反链
exon_plus_cds = exon_plus |> transmute(
  V3 = 'CDS',
  V4 = cds_start,
  V5 = cds_end,
  V6 = '.',
  V7 = V7,
  V8 = '.',
  mRNA = mRNA
)
exon_reverse_cds = exon_reverse |> transmute(
  V3 = 'CDS',
  V4 = cds_start,
  V5 = cds_end,
  V6 = '.',
  V7 = V7,
  V8 = '.',
  mRNA = mRNA
)
exon_cds = rbind(exon_plus_cds, exon_reverse_cds) |>
transmute(
  V1 = str_sub(mRNA, 1, 5),
  V2 = 'YN19',
  V3,
  V4,
  V5,
  V6,
  V7,
  V8,
  V9 = paste0('Parent=', mRNA)
)

##############将cds信息加入之前的fix.gff3,注意过滤掉冗余转录本
gff.mrna.rename = gff.mrna.rename |> filter(!newid %in% rm.mRNA.id)
gff.out = gff.exon |> 
mutate(mRNA = str_extract(V9, 'chr[0-9]{2}g[0-9]{5}.mRNA[0-9]+')) |>
filter(!mRNA %in% rm.mRNA.id)
rbind(gff.gene,
      gff.mrna.rename[, -10:-11],
      gff.exon[, -10],
      exon_cds) 
return(gff.out)
}


identity_micro_exon=function(bam,fa){
  bam=bam
  ref.fa=fa
  ##创建sam表格，包含比对的各类信息，包括reference的位置和序列
  sam = tibble::tibble(
    qname = bam[[1]]$qname,
    flag = bam[[1]]$flag,
    rname = bam[[1]]$rname,
    strand = bam[[1]]$strand,
    pos = bam[[1]]$pos,
    qwidth = bam[[1]]$qwidth,
    cigar = bam[[1]]$cigar,
    seq = as.character(bam[[1]]$seq)
  )
  
  ##过滤read unmapped，not primary alignment，supplementary alignment
  sam=sam|>filter(flag==0 | flag==16)
  
  
  ##利用cigar信息筛选含有微外显子的reads
  ##如果外显子正确匹配，则cigar表现为(n1)M(n2)N(n3)M(n4)N(n5)M
  ##常见的微外显子为单-微外显子，cigar表现为(n1)M(n2)N(n3)I(n4)M，n3<=15
  
  microexon.sam=sam|>filter(str_detect(cigar,'[0-9]+M[0-9]+N[1-9]I[0-9]+M'))
  
  
  ##end比对终止位置
  microexon.sam = microexon.sam |> mutate(end = apply(
    microexon.sam[, c('pos', 'cigar')],
    1,
    get_reference_span,
    start = 'pos',
    cigar = 'cigar'
  ))
  ##ref sequence，后期鉴定内含子
  microexon.sam = microexon.sam |> mutate(refseq = apply(
    microexon.sam[, c('rname','pos', 'end')],
    1,
    function(x) get_reference_seq(x, chr='rname', start='pos', end='end', fa=ref.fa)
  ))
  
  ##insertion pos,找到微外显子的绝对位置
  microexon.sam = microexon.sam |> mutate(insert_pos = apply(
    microexon.sam[, c('pos', 'cigar')],
    1,
    get_insert_pos,
    start = 'pos',
    cigar = 'cigar'
  ))
  
  microexon.sam=microexon.sam|>group_by(insert_pos)|>
    filter(n() > 2)|>
    ungroup()
  
  
  ##排序,正反链分开 
  microexon.sam=microexon.sam|>arrange(strand,rname,pos)
  
  microexon.sam=microexon.sam|>dplyr::select(-flag,-insert_pos)
   return(microexon.sam)
}

###自定义函数，从flag判断正负链###
flag_judge_strand = function(flag) {
  #将10进制flag值转为二进制数
  flag_vector = as.integer(intToBits(flag))
  #flag=16  read reverse strand
  #若flag=16，16=2^4，则二进制数的第4+1位数是1
  if (flag_vector[5] == 1) {
    return('-')
  } else{
    return('+')
  }
}

###自定义函数，基于start和cigar计算reference span end###
get_reference_span = function(x, start, cigar) {
  len = str_remove_all(x[cigar], '\\d+[IS]') |>
    #去除insert和soft clipping的长度
    str_replace_all('[A-Z]', ',') |>
    #逗号分隔，保留数字
    str_split(',') |> unlist() |>
    as.integer() |> sum(na.rm = T)
  #得M,N,D的长度
  
  return(as.integer(x[start]) + len - 1)
  
}

###自定义函数，基于reference span 的 start end获得ref seq###
get_reference_seq = function(x, chr, start, end, fa) {
  start = as.integer(x[start])
  end = as.integer(x[end])
  Biostrings::subseq(fa[x[chr]],
                     start,
                     end) |> toString()
  
}

###自定义函数，基于起始位置和cigar获得I(<=15)的位置###
get_insert_pos = function(x,start,cigar) {
  
  start=as.integer(x[start])
  pos=str_extract_all(as.character(x[cigar]),'\\d+')|>unlist()
  pos=as.integer(pos)
  pos=cumsum(pos)
  pos=pos+start
  cigar=str_extract_all(as.character(x[cigar]),'[A-Z]')|>unlist()
  tbl=data.frame(pos=pos,cigar=cigar)
  tbl=tbl[tbl$cigar=='I',]
  insert_pos=tbl$pos|>paste(sep = ',',collapse = ',')
  return(insert_pos)
}

