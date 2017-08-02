x <- read.table("coverage_B100-6599.bed", header = F, sep = "\t", as.is=T) ##lendo um arquivo para usar na função
x <- dir(".", pattern = ".bed") ##ou lendo todos os arquivos em uma pasta
plot.baseCoverage(x, "SOX3", "mean") ##exemplo



plot.baseCoverage <- function(x, geneID, measure){ ##função e seus argumentos padrões
  if(missing(geneID)) { ##se não for colocado algo no argumento do geneID
    stop("A gene identification is necessary for this function") ##dar um aviso que ele é um argumento necessário para que a função rode
  }
  if(class(x) == "data.frame"){ ##testa se o primeiro objeto inserido com um argumento é um dataframe, que é o objeto criado pelo read.table
    genes.ids <- unique(x$V4) ##cria um objeto com todos os nomes de genes como aparecem no arquivo
    if (any(geneID == genes.ids)){ ##se o argumento dado for encontrado no arquivo, continuar
      gene <- x[x$V4 == geneID,] ##salva apenas a tabela do gene de interesse no objeto, dado como argumento no geneID
      gene.min <- gene[1,2] ##determina a primeira base do gene
      gene.max <- gene[nrow(gene),3] ##determina a ultima base do gene
      gene$exact.base <- gene$V2 + gene$V7 ##cria uma coluna na tabela do gene correspondente a cada uma das bases que estão na tabela original
      dataf1 <- cbind(gene$exact.base, gene$V8) ##cria um dataframe com a cobertura em cada base e sua base correspondente
      gene.extension <- seq(gene.min, gene.max) ##cria um vetor com cada base no gene inteiro (incluindo introns, se não estiverem na tabela original)
      fragments <- c(unique(gene$V2), unique(gene$V3)) ##cria um objeto com o começo e o final de todos os fragmentos, ou seja, inicio e fim de todos os exons
      fragments <- sort(fragments) ##ordenar os fragmentos em ordem numérica
      dataf2 <- NULL ##cria um objeto nulo para:
      dataf2 <- cbind(dataf2, gene.extension) ##cria um dataframe com todas as bases do gene inteiro. a ideia é que se tiverem bases que não tiveram leitura, não se tenha um ponto para elas. por isso:
      dataf3 <- merge(dataf2, dataf1, by.x = "gene.extension", by.y = "V1", all = T) ##junta os dois dataframes anteriores criados em um só pela base, colocando NA nas bases que não possuem leitura
      par(las = 1, bty = "l") ##especifica para o gráfico o tipo de caixa, números na horizontal
      plot(dataf3, type = "l", xaxt = "n", ##plota o gráfico final, com uma linha ligando todos os pontos, sem marcas no eixo X
         ylim = range(0:max(gene$V8)), ##padroniza o eixo Y para que seja possível se ver quais bases chegam até 100,
         xlab = geneID, ylab = "Coverage", main = "Sample") ##nomeia os eixos X e Y e o gráfico
      axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
    } else { ##caso o gene dado no argumento não seja encontrado no arquivo
      stop("The geneID was not found in the file! Please make sure it is the same as appears on your file or it is not misspeled!") ##retornar mensagem de erro
    }
  } else { ##caso o primeiro if retorne negativo, entrar na opção de vários arquivos
    ##criação de objetos que serão utilizados após o for
    coverages <- NULL ##para salvar as coberturas
    gene.extension <- NULL ##para salvar a extensão do gene
    exact.base <- NULL ##para salvar a base exata de onde se tem cobertura
  
    for (i in x){ ##para cada um dos arquivos cujo nomes estão em x
      table <- read.table(i, header = F, sep = "\t", as.is = T) ##leitura dos arquivos como dataframes
      genes.ids <- unique(table$V4) ##cria um objeto com todos os nomes de genes como aparecem no arquivo
     if (any(geneID == genes.ids)){ ##se o argumento dado for encontrado no arquivo, continuar
        gene <- table[table$V4 == geneID,] ##salva apenas a tabela do gene de interesse dado no argumento geneID no objeto
        gene.min <- gene[1,2] ##determina a primeira base do gene
        gene.max <- gene[nrow(gene),3] ##determina a ultima base do gene
        exact.base <- gene$V2 + gene$V7 ##cria um vetor que possua cada uma das bases que possuem cobertura
        fragments <- c(unique(gene$V2), unique(gene$V3)) ##cria um objeto com todos os inícios e fins dos fragmentos, ou seja, o inicio e o fim dos exons, que são a parte que terão cobertura no caso do exoma
        fragments <- sort(fragments) ##ordena esse objeto em ordem numérica
        coverages <- cbind(coverages, gene$V8) ##cria dataframe apenas com as coberturas de todos os arquivos
        gene.extension <- seq(gene.min, gene.max) ##cria um vetor com cada base no gene inteiro (incluindo introns, se eles não estiverem na tabela original)
    } else { ##caso o gene dado no argumento não seja encontrado no arquivo
        stop("The geneID was not found in the file! Please make sure it is the same as appears on your file or it is not misspeled!") ##retornar mensagem de erro
      }
        }
    dataf2 <- NULL ##cria um objeto nulo para:
    dataf2 <- cbind(dataf2, gene.extension) ##transforma o gene.extension criado anteriormente em dataframe
    par(las = 1, bty = "l") ##padroniza para todos os gráficos feitos o tipo de caixa e os números nas horizontais nos eixos
    if(missing(measure)){ ##se não tiver o argumento med, retornar os 3 gráficos
        media <- apply(coverages, 1, mean) ##aplicar a função média escolhida no argumento
        dataf1a <- cbind(exact.base, media) #juntar as médias com as bases aonde se encontram em um dataframe
        dataf3a <- merge(dataf2, dataf1a, by.x = "gene.extension", by.y = "exact.base", all = T) ##juntar os dois dataframes criados em um só pela base, afim de se ter um dataframe no qual não se tem nada nas bases que não tiveram cobertura nenhuma
        plot(dataf3a, type = "l", xaxt = "n", ##plota o gráfico com a média das bases e retira as marcas e números no eixo X
             ylim = range(0:max(media)), ##padroniza o eixo Y
             xlab = "geneID", ylab = "Coverage", main = "Mean Coverage") ##nomeia o gráfico e seus eixos
        axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
        mediana <- apply(coverages, 1, median) ##aplicar a função média escolhida no argumento
        dataf1b <- cbind(exact.base, mediana) #juntar as médias com as bases aonde se encontram em um dataframe
        dataf3b <- merge(dataf2, dataf1b, by.x = "gene.extension", by.y = "exact.base", all = T) ##juntar os dois dataframes criados em um só pela base, afim de se ter um dataframe no qual não se tem nada nas bases que não tiveram cobertura nenhuma
        par(ask=T) ##pergunta antes de plotar os gráficos, pedindo para que o enter seja apertado
        plot(dataf3b, type = "l", xaxt = "n", ##plota o gráfico com a mediana e retira os números do eixo X
             ylim = range(0:max(mediana)), ##padroniza o eixo Y
             xlab = geneID, ylab = "Coverage", main = "Median Coverage") ##nomeia o gráfico e seus eixos
        axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
        dataf1c <- cbind(exact.base, coverages) #juntar as médias com as bases aonde se encontram em um dataframe
        plot(NULL, ##plota o gráfico vazio
             ylim = range(0:max(coverages)), xlim = range(gene.min, gene.max), ##especifica os limites dos eixos X e Y
             xaxt = "n", ##retira as marcas do eixo X
             xlab = geneID, ylab = "Coverage", main = "Variance") ##nomeia o gráfico e os dois eixos
        axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
        for(i in 2:ncol(dataf1c)){ ##entra em um for para plotar todos os pontos de cobertura encontrados no dataframe
          points(dataf1c[,1], dataf1c[,i], pch = 20) #plota cada um dos pontos no gráfico
        }
        lines(dataf3a, col = "red") ##plota a linha de média de referência
      }
      else if(measure == "mean"){ ##ou, caso o argumento dado em med for "mean", ou seja, a média
        media <- apply(coverages, 1, mean) ##aplica a função média a cada linha do dataframe, para determinar a média de cobertura de cada base 
        dataf1 <- cbind(exact.base, media) #junta as médias com suas bases em um dataframe
        dataf3 <- merge(dataf2, dataf1, by.x = "gene.extension", by.y = "exact.base", all = T) ##juntar os dois dataframes criados em um só pela base, afim de se ter um dataframe no qual as bases que não tiveram cobertura são NAs
        plot(dataf3, type = "l", xaxt = "n", ##plota o gráfico com a média das bases e retira os números do eixo X
           ylim = range(0:max(media)), ##padroniza o eixo Y
           xlab = geneID, ylab = "Coverage", main = "Mean Coverage") ##nomeia o gráfico e seus eixos
        axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
      }
      else if(measure == "median"){ ##ou, caso o argumento dado em med for "median", ou seja, mediana
        mediana <- apply(coverages, 1, median) ##aplica a função mediana em cada linha do dataframe, e guarda em um objeto com apenas essas medianas
        dataf1 <- cbind(exact.base, mediana) #junta as medianas e as bases correspondentes em um dataframe
        dataf3 <- merge(dataf2, dataf1, by.x = "gene.extension", by.y = "exact.base", all = T) ##junta os dois dataframes criados em um só pela base, afim de se ter um dataframe no qual as bases que não tiveram coberturas apresentam NAs
        plot(dataf3, type = "l", xaxt = "n", ##plota o gráfico com a mediana  e retira os números do eixo X
           ylim = range(0:max(mediana)), ##padroniza o eixo Y
           xlab = geneID, ylab = "Coverage", main = "Median Coverage") ##nomeia o gráfico e seus eixos
        axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
      }
      else if(measure == "variance" || measure == "var"){ ##ou, caso o argumento dado em med for "var", ou seja, variancia
        dataf1 <- cbind(exact.base, coverages) #juntar as coberturas com as bases
        media <- apply(coverages, 1, mean) ##aplica a função média a cada linha do dataframe, para determinar a média de cobertura de cada base 
        datafm <- cbind(exact.base, media) #junta as médias com suas bases em um dataframe
        datafm <- merge(dataf2, datafm, by.x = "gene.extension", by.y = "exact.base", all = T) ##juntar os dois dataframes criados em um só pela base, afim de se ter um dataframe no qual as bases que não tiveram cobertura são NAs
        plot(NULL, ##plota o gráfico vazio
           ylim = range(0:max(coverages)), xlim = range(gene.min, gene.max), ##especifica os limites dos eixos X e Y
           xaxt = "n", ##retira as marcas do eixo X
           xlab = geneID, ylab = "Coverage", main = "Variance") ##nomeia o gráfico e os dois eixos
        axis(1, at = fragments, labels = F) ##coloca marcas em cada começo e fim dos fragmentos no eixo X do gráfico
        for(i in 2:ncol(dataf1)){ ##entra em um for para plotar todos os pontos de cobertura encontrados no dataframe
          points(dataf1[,1], dataf1[,i], pch = 20) #plota cada um dos pontos no gráfico
        }
        lines(datafm, col = "red") ##plota a linha de média de referência
      }
  }
}
