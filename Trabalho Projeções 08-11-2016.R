#Carrega funções auxiliares
{
  projetarTFT = function (k1,k2,TFT0,TFT,T) {
    t = seq(0,T,5)
    k = k1+k2
    a = log((k-TFT0)/(TFT0-k1))
    b = (1/T)*(log((k-TFT)/(TFT-k1))-a)
    TFTt = k1 + (k2/(1+exp(a+b*t)))
    return(TFTt)
  }
  
  projetarTEF = function (TEFInicial,TEFPadrao,TFTProjetada,AnoInicial,AnoFinal,AnoProjetado) {
    #Calcula o acumulado (Fx) das TEF da população analisada
    TEFInicialAcumulada = cumsum(TEFInicial)[1:6]*5
    
    #Calcula o acumulado (Fx) das TEF da população escolhida como padrão
    TEFPadraoAcumulada = cumsum(TEFPadrao)[1:6]*5
    
    #Calcula a TFT da população observada e da população padrão atrvés dos somatórios das TEF
    TFTInicial = sum(TEFInicial)*5
    TFTPadrao = sum(TEFPadrao)*5
    
    #Calcula o Vx, V1 e V2 da população observada e padrão, usando o método de Gompertz
    VxInicial = log(-log(TEFInicialAcumulada[1:6]/TFTInicial))
    VxPadrao = log(-log(TEFPadraoAcumulada[1:6]/TFTPadrao))
    V1Inicial = sum(VxInicial[1:3])/3
    V2Inicial = sum(VxInicial[4:6])/3
    V1Padrao = sum(VxPadrao[1:3])/3
    V2Padrao = sum(VxPadrao[4:6])/3
    
    #Calcula os parâmetros Alfa e Beta para a população observada
    BETAi = (V2Inicial-V1Inicial)/(V2Padrao-V1Padrao)
    ALFAi = V1Inicial - BETAi*V1Padrao
    
    #Assume-se que para população padrão, Alfa é igual a 0 e Beta é igual a 1
    ALFAs = 0
    BETAs = 1
    
    #Calcula o Alfa e Beta
    ALFAx = ((AnoProjetado-AnoInicial)*ALFAs + (AnoFinal-AnoProjetado)*ALFAi)/(AnoFinal-AnoInicial)
    
    BETAx = ((AnoProjetado-AnoInicial)*BETAs + (AnoFinal-AnoProjetado)*BETAi)/(AnoFinal-AnoInicial)
    
    #Estima os valores de V
    Vx = ALFAx + BETAx*VxPadrao
    
    #Calcula projeção das TEF
    TEFProjetadaAcumulada = exp(-exp(Vx))*TFTProjetada
    TEFProjetadaAcumulada[7] = TFTProjetada
    TEFProjetada=TEFProjetadaAcumulada[1]/5
    TEFProjetada[2:7] = diff(TEFProjetadaAcumulada)/5
    
    return(TEFProjetada)
  }
  
  projetaMortalidade = function(AnoInicial,AnoFinal,lxInicial,lxFinal) {
    AnosProjetados = seq(AnoInicial+5,AnoFinal-5,5)
    NGrupos = nrow(lxInicial)
    NProjecoes = length(AnosProjetados)
    #Calcula logito inicial e logito final
    LogitoInicial = 0.5*log((1-lxInicial)/lxInicial)
    LogitoFinal = 0.5*log((1-lxFinal)/lxFinal)
    
    LogitoProjetado = matrix (0,nrow=NGrupos,ncol=NProjecoes)
    #Estima logito para o tempo projetado
    for (i in 1:NProjecoes) {
      LogitoProjetado[,i] = (((AnosProjetados[i]-AnoInicial)*LogitoFinal+(AnoFinal-AnosProjetados[i])*LogitoInicial)/(AnoFinal-AnoInicial))[,1]
    }
    #Projeta lx para o ano escolhido
    lxProjetado = 1/(1+exp(2*LogitoProjetado))
    
    #Acrescenta lx inicial e lx padrão na matriz
    lxProjetado = cbind(lxInicial,lxProjetado,lxFinal)
    
    #Transforma em um data.frame
    lxProjetado = as.data.frame(lxProjetado)
    names(lxProjetado) = as.character(c(AnoInicial,AnosProjetados,AnoFinal))
    GruposEtarios = seq(0,(NGrupos-1)*5,5)
    row.names(lxProjetado) = as.character(paste(GruposEtarios,'a',GruposEtarios+4))
    
    return(lxProjetado)
  }
  
  #Função de calcular tábua de vida precisa de ajustes
  calculaTabuaDeVida = function(lx) {
    #Transforma lx em um vetor
    lx = lx[,1]
    
    #Quantidade de grupos etários para a tábua de vida
    NGrupos = length(lx)
    
    #Calcula px e qx (probabilidade de sobrevivência e morte)
    px = 0
    for (i in 1:(NGrupos-1)) {
      px[i] = lx[i+1]/lx[i]
    }
    px[NGrupos]=0
    qx = 1-px
    
    #Calcula a taxa de mortalidade (mx). **OBS: FALTA UM TRATAMENTO ESPECIAL PARA O CÁLCULO DO ÚLTIMO GRUPO ETÁRIO!**
    mx = -log(1-qx)
    mx[NGrupos] = mx[NGrupos-1]
    
    #Calcula dx
    dx = rev(diff(rev(lx)))
    dx[NGrupos] = lx[NGrupos]
    
    #Cria vetor com valores de n e estima ax. **OBS: FALTA UM TRATAMENTO ESPECIAL PARA O CÁLCULO DO ax[1]**
    n = rep(5,NGrupos)
    n[NGrupos] = Inf
    ax = rep(2.5,NGrupos)
    ax[1] = 1.651-(2.816*(mx[1]/5))
    
    #Calcula Lx
    Lx = rep(0,NGrupos)
    Lx[1:(NGrupos-1)] = lx[2:NGrupos]*n[1:(NGrupos-1)] + dx[1:(NGrupos-1)]*ax[1:(NGrupos-1)]
    ################################FALTA AJUSTES NO LX PARA O PRIMEIRO E ULTIMO GRUPO ETARIO
    Lx[NGrupos] = 2*lx[NGrupos]/qx[NGrupos-1]
    
    #Calcula Tx
    Tx = rep(0,NGrupos)
    Tx[NGrupos] = Lx[NGrupos]
    for (i in (NGrupos-1):1) {
      Tx[i] = Tx[i+1]+Lx[i]
    }
    
    #Calcula esperança de vida
    ex = Tx/lx
    
    #Cria data.frame com todas as informações
    x = seq(0,(NGrupos*5)-1,5)
    TabuaDeVida = data.frame(x,n,qx,lx,ax,dx,Lx,Tx,ex)
    
    return(TabuaDeVida)
  }
  
  projetaPopulacao5Anos = function (PopulacaoInicial,Lx,TEFProjetada) {
    #Quantidade de Grupos Etarios
    NGrupos = length(Lx)
    
    #Calcula razão de sobrevivência (Sx)
    Sx = rep(0,(NGrupos-1))
    Sx[1:(NGrupos-1)] = Lx[2:NGrupos]/Lx[1:(NGrupos-1)]
    #Cria matriz de sobrevivência
    MatrizSobrevivencia = matrix(0,nrow=NGrupos,ncol=NGrupos)
    diag(MatrizSobrevivencia[2:NGrupos,]) = Sx
    MatrizSobrevivencia[1,4:10] = TEFProjetada*5
    #if(b1==FALSE) {  #Ignorar essa parte. Apenas para testes.
    #  m1 <<- MatrizSobrevivencia
    #  b1 <<- TRUE
    #}
    #Projeta a população 5 anos adiante, multiplicando a matriz sobreviência pelo vetor de número de pessoas
    PopulacaoProjetada = round(MatrizSobrevivencia%*%as.matrix(PopulacaoInicial))
    
    return(PopulacaoProjetada)
  }
  
  #Função Principal
  projetaPopulacao = function (k1,k2,TFTInicial,TFTFinal,TEFInicial,TEFPadrao,AnoInicial,AnoFinal,
                               lxInicialMasc,lxFinalMasc,lxInicialFem,lxFinalFem,PopulacaoInicialMasc,
                               PopulacaoInicialFem,RazaoDeSexo) {
    
    TempoProjetado = AnoFinal - AnoInicial
    
    #Projeta nível da fecundidade
    TFTProjetada = projetarTFT(k1,k2,TFTInicial,TFTFinal,TempoProjetado)
    
    #Projeta estrutura da fecundidade
    TEFProjetada = matrix(NA,nrow=7,ncol=0)
    AnosProjetados = seq(AnoInicial+5,AnoFinal,5)
    j=2
    for (i in AnosProjetados) {
      TEFProjetada = cbind(TEFProjetada, projetarTEF(TEFInicial,TEFPadrao,TFTProjetada[j],AnoInicial,AnoFinal,i))
      j=j+1
    }
    TEFProjetada = as.data.frame(TEFProjetada)
    names(TEFProjetada) = as.character(AnosProjetados)
    GruposEtarios = seq(15,45,5)
    row.names(TEFProjetada) = as.character(paste(GruposEtarios,'a',GruposEtarios+4))
    
    
    #Projeta mortalidade
    lxMasc = projetaMortalidade(AnoInicial,AnoFinal,(lxInicialMasc),(lxFinalMasc))
    lxFem = projetaMortalidade(AnoInicial,AnoFinal,(lxInicialFem),(lxFinalFem))
    
    #Cria tábua de vida para cada lx e obtém separadamente cada Lx
    TabuaDeVidaMasc = TabuaDeVidaFem = list()
    LxMasc = LxFem = matrix(NA,nrow=nrow(lxMasc),ncol=0)
    for (i in 1:length(lxMasc)) {
      TabuaDeVidaMasc[[i]] = calculaTabuaDeVida(lxMasc[i])
      LxMasc = cbind(LxMasc,TabuaDeVidaMasc[[i]]$Lx)
      TabuaDeVidaFem[[i]] = calculaTabuaDeVida(lxFem[i])
      LxFem = cbind(LxFem,TabuaDeVidaFem[[i]]$Lx)
    }
    LxMasc = as.data.frame(LxMasc)
    LxFem = as.data.frame(LxFem)
    names(TabuaDeVidaMasc) = names(TabuaDeVidaFem) = names(LxMasc) = names(LxFem) = names(lxMasc)
    GruposEtarios = seq(0,(nrow(LxMasc)-1)*5,5)
    row.names(LxMasc) =row.names(LxFem) = as.character(paste(GruposEtarios,'a',GruposEtarios+4))
    
    #Projeta população
    PopulacaoProjetadaMasc = PopulacaoProjetadaFem = matrix(NA,nrow=nrow(LxMasc),ncol=0)
    PopulacaoProjetadaMasc = cbind(PopulacaoProjetadaMasc,PopulacaoInicialMasc)
    PopulacaoProjetadaFem = cbind(PopulacaoProjetadaFem,PopulacaoInicialFem)
    NascimentosTotais = NascimentosMasc = NascimentosFem = 0
    for (i in 1:(length(LxMasc)-1)) {
      PopulacaoProjetadaFem = cbind(PopulacaoProjetadaFem,projetaPopulacao5Anos(PopulacaoProjetadaFem[i],LxFem[,i+1],TEFProjetada[,i]))
      PopulacaoProjetadaMasc = cbind(PopulacaoProjetadaMasc,projetaPopulacao5Anos(PopulacaoProjetadaMasc[i],LxMasc[,i+1],rep(0,nrow(TEFProjetada))))
      
      #Calcula quantos dos nascimentos foram masculinos e femininos, e preenche na matriz
      NascimentosTotais = PopulacaoProjetadaFem[1,i+1]
      NascimentosMasc = (RazaoDeSexo/(1+RazaoDeSexo))*NascimentosTotais
      NascimentosFem = NascimentosTotais-NascimentosMasc
      PopulacaoProjetadaMasc[1,i+1] = NascimentosMasc
      PopulacaoProjetadaFem[1,i+1] = NascimentosFem
    }
    
    #Transforma a matriz de populações projetadas em um data.frame
    PopulacaoProjetadaMasc = as.data.frame(PopulacaoProjetadaMasc)
    PopulacaoProjetadaFem = as.data.frame(PopulacaoProjetadaFem)
    names(PopulacaoProjetadaMasc) = names(PopulacaoProjetadaFem) = as.character(seq(AnoInicial,AnoFinal,5))
    row.names(PopulacaoProjetadaMasc) = row.names(PopulacaoProjetadaFem) = as.character(paste(GruposEtarios,'a',GruposEtarios+4))
    row.names(PopulacaoProjetadaMasc)[nrow(PopulacaoProjetadaMasc)] = row.names(PopulacaoProjetadaFem)[nrow(PopulacaoProjetadaFem)] = paste(GruposEtarios[nrow(PopulacaoProjetadaMasc)],'+')
    
    return(list(TFTProjetada=TFTProjetada,TEFProjetada=TEFProjetada,TabuaDeVidaMasc=TabuaDeVidaMasc,
                TabuaDeVidaFem=TabuaDeVidaFem,PopulacaoProjetadaMasc=PopulacaoProjetadaMasc,
                PopulacaoProjetadaFem=PopulacaoProjetadaFem))
  }
}

#Carrega base de dados já inseridos via dput (RECOMENDADO. Neste caso, não utilize o código de leitura de arquivos posterior)
{  
  TEFSP2010 = c(0.0538, 0.085, 0.0836, 0.069, 0.0384, 0.01, 6e-04)
  TEFCanada2000 = c(0.016962, 0.058912, 0.098526, 0.086446, 0.034296, 0.00579, 0.00023)
  TEFTaiwan2001 = c(0.006674, 0.050656, 0.101596, 0.087024, 0.02916, 0.004406, 2e-04)
  lxSP2010Masc = structure(c(1, 0.9965724442, 0.995157311329236, 0.990101912187684, 
                             0.982458325425595, 0.974785325904021, 0.965183690443866, 0.952337095524058, 
                             0.934461728241072, 0.908109907504674, 0.869606047426475, 0.817647086092744, 
                             0.750158495606648, 0.660116971378983, 0.545573474505301, 0.404842796756659, 
                             0.260613501984132, 0.131130289658336), .Dim = c(18L, 1L))
  lxSP2010Fem = structure(c(1, 0.986480184, 0.98547397421232, 0.983611428401059, 
                            0.981595024972836, 0.979141037410404, 0.975655295317223, 0.970533105016808, 
                            0.962293278955215, 0.948965517041685, 0.929891310149147, 0.901669108886121, 
                            0.859858712307071, 0.798206842634654, 0.708711891438457, 0.583808507691343, 
                            0.421404657021765, 0.240799049115377), .Dim = c(18L, 1L))
  lxCanada2000Masc = structure(c(1, 0.9944046716, 0.993718532376596, 0.992724813844219, 
                                 0.989419040214118, 0.985372316339642, 0.981558925475408, 0.976710024383559, 
                                 0.970566518330187, 0.961744068678565, 0.948943255124453, 0.92802854578151, 
                                 0.892949066750969, 0.836541474204311, 0.748763177316052, 0.623233030639016, 
                                 0.459553339802291, 0.274261433194008), .Dim = c(18L, 1L))
  lxCanada2000Fem = structure(c(1, 0.9953033616, 0.994726085650272, 0.994079513694599, 
                                0.992687802375427, 0.991357600720244, 0.989781342135099, 0.987762188197143, 
                                0.984630982060558, 0.979520747263664, 0.971312363401594, 0.958763007666446, 
                                0.937420943115791, 0.903439433927843, 0.850922499633618, 0.770033806818446, 
                                0.647028606517267, 0.470985063256049), .Dim = c(18L, 1L))
  PopSP2010Masc = structure(list(V1 = c(1361616L, 1457203L, 1687826L, 1667482L, 
                                        1835222L, 1881495L, 1741346L, 1549270L, 1444230L, 1308853L, 1149501L, 
                                        930303L, 705940L, 499180L, 371655L, 246532L, 150452L, 89767L)), .Names = "V1", class = "data.frame", row.names = c(NA,-18L))
  PopSP2010Fem = structure(list(V1 = c(1313756L, 1403430L, 1637087L, 1636426L, 
                                       1802466L, 1908294L, 1815101L, 1634851L, 1536444L, 1444270L, 1286603L, 
                                       1057688L, 831069L, 609906L, 484550L, 354796L, 246113L, 181476L)), .Names = "V1", class = "data.frame", row.names = c(NA, -18L))
  lxSP2000Masc = structure(c(1, 0.9763877196, 0.974073680704548, 0.962433500220129, 
                             0.945812273671327, 0.928428244081248, 0.909609003573721, 0.888287768529953, 
                             0.861727964250907, 0.826931389054456, 0.782161323651048, 0.721160562019502, 
                             0.642604541998718, 0.542338955310658, 0.42500934571875, 0.293723958826228, 
                             0.166412246112588, 0.0714457696235174), .Dim = c(18L, 1L))
  lxSP2000Fem = structure(c(1, 0.9800132748, 0.978651056348028, 0.976165282664904, 
                            0.973139170288643, 0.96898386603151, 0.963383139285848, 0.955743510991311, 
                            0.944494409866944, 0.927701299259509, 0.903172876907088, 0.867172406033571, 
                            0.814595743055756, 0.741689424052266, 0.642392043960148, 0.505748832289385, 
                            0.334739979627375, 0.170479724224426), .Dim = c(18L, 1L))
}

#Carregar base de dados pela leitura de arquivos (Só utilize esse código abaixo caso queira carregar dados diferentes do utilizado no trabalho)
{
  #Função auxiliar para rearranjar dados
  convertePxParalx = function (px) {
    Npx = length(px)
    px = c(px[1:18],px[Npx])
    Npx = length(px)
    px = c(px[1]*px[2],px[3:Npx])
    Npx = length(px)
    lx = 1
    for (i in 2:(Npx)){
      lx[i] = lx[i-1]*px[i-1]
    }
    
    return(as.matrix(lx))
  }
  TEFCanada2000 = read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/TEFCanada2000.txt")/5
  TEFTaiwan2001TodasIdades = read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/TEFTaiwan2001TodasIdades.txt")
  TEFTaiwan2001 = c(sum(TEFTaiwan2001TodasIdades[1:4,]),sum(TEFTaiwan2001TodasIdades[5:9,]),sum(TEFTaiwan2001TodasIdades[10:14,]),sum(TEFTaiwan2001TodasIdades[15:19,]),sum(TEFTaiwan2001TodasIdades[20:24,]),sum(TEFTaiwan2001TodasIdades[25:29,]),sum(TEFTaiwan2001TodasIdades[30:34,]))/5
  TEFSP2010 = read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/TEFSP2010.txt")
  TabuaDeVidaCanada2000 = read.csv("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/canada.txt")
  lxCanada2000Masc = convertePxParalx(px = 1-TabuaDeVidaCanada2000$q.x.[TabuaDeVidaCanada2000$Sex==1 & TabuaDeVidaCanada2000$TypeLT==5]) #Faz um pequeno ajuste para ajustar os grupos quinquenais
  lxCanada2000Fem = convertePxParalx(px = 1-TabuaDeVidaCanada2000$q.x.[TabuaDeVidaCanada2000$Sex==2 & TabuaDeVidaCanada2000$TypeLT==5]) #Faz um pequeno ajuste para ajustar os grupos quinquenais
  lxSP2010Masc = convertePxParalx(read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/pxSP2010Masc.txt")[,1])
  lxSP2010Fem = convertePxParalx(read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/pxSP2010Fem.txt")[,1])
  PopSP2010Masc = read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/PopSP2010Masc.txt")
  PopSP2010Fem = read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/PopSP2010Fem.txt")
  
  lxSP2000Masc = convertePxParalx(read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/pxSP2000Masc.txt")[,1])
  lxSP2000Fem = convertePxParalx(read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/pxSP2000Fem.txt")[,1])
  lxBrasil2060Masc = convertePxParalx(read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/pxBrasil2060Masc.txt")[,1])
  lxBrasil2060Fem = convertePxParalx(read.table("C:/Users/PC/OneDrive/UFRN/2016.2/Projecoes Populacionais/Trabalho Projecao Populacao no R/pxBrasil2060Fem.txt")[,1])
}

#Salva resultado da projeção, em uma lista, onde o primeiro elemento é TFT Projetada,
#o segundo a TEF projetada, o terceiro é uma lista com todas as tábuas de vida masculina projetadas,
#o quarto é uma lista com todas as tábuas de vida femininas projetadas, o quinto é
#a população masculina projetada, e o sexto é a população feminina projetada.
#Nas entradas, os dados TEFInicial e TEFPadrao devem estar no formato de um vetor simples,
#k1, k2, TFTInicial, TFTFinal, AnoInicial, AnoFinal e RazaoDeSexo são valores únicos e simples,
#os demais argumentos devem ser passaados em formato de matriz coluna.
ResultadosProjecoes = projetaPopulacao(k1=1.3, k2=3.39, TFTInicial=1.66, TFTFinal=1.476, TEFInicial=TEFSP2010, TEFPadrao=TEFCanada2000,
                 AnoInicial=2010, AnoFinal=2030, lxInicialMasc=lxSP2010Masc, lxFinalMasc=lxCanada2000Masc,
                 lxInicialFem=lxSP2010Fem,  lxFinalFem=lxCanada2000Fem, PopulacaoInicialMasc=PopSP2010Masc,
                 PopulacaoInicialFem=PopSP2010Fem, RazaoDeSexo=0.9478)

#Projeção apenas da mortalidade (lx) de São Paulo em 2000 para 2030, utilizando como padrão
#a mortalidade do Brasil projetada para 2060 segundo dados do IBGE
{
  lxProjetadoSPMasc = projetaMortalidade(AnoInicial = 2000,AnoFinal = 2030,lxInicial = lxSP2000Masc,lxFinal = lxBrasil2060Masc)
  lxProjetadoSPFem = projetaMortalidade(AnoInicial = 2000,AnoFinal = 2030,lxInicial = lxSP2000Fem,lxFinal = lxBrasil2060Fem)
  
  calculaTabuaDeVida(lxProjetadoSPMasc[6])
}

#Plota gráficos das projeções da estrutura de fecundidade
{
  png(filename = "Projecao Fecundidade em São Paulo.png")
  plot(ResultadosProjecoes$TEFProjetada$`2030`,xaxt='n',type='n',main="Projeção da estrutura da fecundidade \n para o estado de São Paulo",xlab="Idade",ylab="Taxa Específica de Fecundidade")
  lines(TEFSP2010,type='l',col='blue',lwd=2)
  lines(ResultadosProjecoes$TEFProjetada$`2015`,type='l',col='red',lwd=2)
  lines(ResultadosProjecoes$TEFProjetada$`2020`,type='l',col='green',lwd=2)
  lines(ResultadosProjecoes$TEFProjetada$`2025`,type='l',col='blueviolet',lwd=2)
  lines(ResultadosProjecoes$TEFProjetada$`2030`,type='l',col='cyan',lwd=2)
  axis(1,at=1:7,labels=seq(15,45,5))
  abline(v=c(1:7),lty=10,col='gray',pch=23)
  abline(h=c(seq(0.02,0.08,0.02)),lty=10,col='gray',pch=23)
  legend(5,0.09,legend=c(as.character(seq(2010,2030,5))),col=c('blue','red','green','blueviolet','cyan'),lty=1,cex=1,lwd=3)
  dev.off()
}

attach(ResultadosProjecoes)

#Gera pirâmides etárias
{
  library(plotrix)
  mcol<-color.gradient(c(0,0,0.5,1),c(0,0,0.5,1),c(1,1,0.5,1),18)
  fcol<-color.gradient(c(1,1,0.5,1),c(0.5,0.5,0.5,1),c(0.5,0.5,0.5,1),18)
  
  for (i in 1:length(PopulacaoProjetadaMasc)) {
    png(filename=paste("Pirâmide Etária de SP em",names(PopulacaoProjetadaMasc)[i],".png"))
    PopulacaoTotal = sum(c(PopulacaoProjetadaMasc[,i],PopulacaoProjetadaFem[,i]))/100
    par(mar=pyramid.plot(PopulacaoProjetadaMasc[,i]/PopulacaoTotal,PopulacaoProjetadaFem[,i]/PopulacaoTotal,
                         main=paste("Pirâmide Etária do Estado de São Paulo em",names(PopulacaoProjetadaFem)[i]),
                         top.labels = c('Homens','Idade','Mulheres'),labels=row.names(PopulacaoProjetadaMasc),
                         lxcol=mcol,rxcol=fcol,gap=1,show.values=TRUE))
    dev.off()
  }
}

