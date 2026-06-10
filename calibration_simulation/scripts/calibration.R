# Carregar o pacote necessário
library(ape)

# 1. CRIAR A ÁRVORE REAL (O Cenário Verdadeiro na Natureza)
# Vamos definir uma árvore onde a separação ocorreu há exatamente 10 milhões de anos (Ma).
# Na escala de tempo real (Cronograma), os ramos têm comprimentos iguais até o presente.
crono = read.tree(text = "((A:2, B:2):8, (C:5, D:5):5) ;")

# 2. SIMULAR O FILOGRAMA (O que o software lê a partir do DNA)
# Imaginemos que a linhagem A evoluiu 3 vezes mais rápido que as outras.
# O filograma resultante (comprimento = tempo * taxa) terá o ramo A muito mais longo.
filo = crono
filo$edge.length = c(2, 2, 4, 1, 1, 2) # Ramo de A virou 30 devido à taxa alta

# 3. DATAR USANDO RELÓGIO ESTRITO (Strict Clock)
# O relógio estrito assume taxa constante. Vamos usar o método K槍sch (chronos) 
# forçando uma rigidez total (model = "strict") e calibrando a raiz em 10 Ma.
relogio_estrito = chronos(phy = filo, 
                           model = "clock", 
                           calibration = makeChronosCalib(
                             phy = filo, 
                             node = "root", 
                             age.min = 10, 
                             age.max = 10)
                           )

# 4. DATAR USANDO RELÓGIO RELAXADO (Relaxed Clock)
# Vamos permitir que as taxas variem entre os ramos usando o modelo "correlated" 
# (similar ao espírito de um relógio relaxado no BEAST, onde as taxas variam)
relogio_relaxado = chronos(phy = filo, 
                            model = "correlated", 
                            calibration = makeChronosCalib(
                              phy = filo, 
                              node = "root", 
                              age.min = 10, 
                              age.max = 10)
                            )

# 5. PLOTAR OS RESULTADOS PARA COMPARAR
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot(crono, main = "1. Tempo Real (Alvo)\nNó (A,B) = 10 Ma", srv = FALSE)
axisPhylo()

plot(relogio_estrito, main = "2. Relógio Estrito\nNó (A,B) Distorcido", srv = FALSE)
axisPhylo()

plot(relogio_relaxado, main = "3. Relógio Relaxado\nNó (A,B) Corrigido", srv = FALSE)
axisPhylo()
