import numpy as np
print(np.random.random_integers(0, 1))

# arr = np.zeros(9)
# arr.shape = (3,3)
# print(arr)
# from keras.models import Sequential
# from keras.layers import Dense

# model = Sequential()
# model.add(Dense(units=64,activation='relu', input_dim=2))
# model.add(Dense(units=10, activation='softmax'))

# model.compile(loss='categorical_crossentropy',
#               optimizer='sgd',
#               metrics=['accuracy'])

# # model.fit(x_train, y_train, epochs=5, batch_size=32)
# ALGORITMO GENÉTICO
GERACAO = 1
INDIVIDUOS = 60
GENE = 516
NINTER = 10000

# REDE NEURAL
TIPO = 4
N = 4
NINT = 4
NINT1 = 7
NINT2 = 7
NENT = 5
NSAI = 1
NAPD = 900
TARP3 = 0.01
TARP4 = 0.01
MAXINTER = 5000000
ACEITAVEL = 0.001
CONTID = 0
CONTREGRASVAL = 0
INDIC = 0
SUP = 0
NEURIT = 0

# TAXAS
TX_TRANS = 0.001
ERR = 0.0
EMQ = 0.0
FITNESS = 0.0
AUX = 0
SOMA_NINT_TOTAL = 0

# CONSTANTES
FIT = 0
HIST_FIT = 0

tab_converte = ['f', 'F', 'n', '.', 'n', '.', 'f', 'F', 'F', 'f', 'B', 'f', '[', 'n', '[', '.',
                'f', ']', 'n', '*', '.', 'F', 'f', 'F', ']', '.', '[', 'f', 'f', '*', 'B', ']',
                '.', ']', 'n', 'F', 'f', 'B', 'f', 'B', 'F', '[', 'B', 'n', '*', 'f', '.', ']',
                ']', '[', 'n', 'F', 'n', 'B', '[', '.', 'f', ']', 'B', 'F', 'B', 'f', '*', '[']


def zerar_fitness(fit):
    fit = np.zeros(INDIVIDUOS * GENE)
    fit.shape = (INDIVIDUOS, GENE)
    return fit


def imprime_hist_fitness(hist_fit, file, contador):
    # imprimir no txt
    print("\n", contador)
    print("\t\t", hist_fit[contador-1][0])
    print("\t\t", hist_fit[contador-1][1])


def genotipo_estatico(individuos, gene):
    gen = np.zeros(individuos * gene)
    gen.shape = (individuos, gene)
    for i in range(0, gene):
        for j in range(0, individuos):
            gen[i][j] = np.random.random_integers(0, 1)
    return gen


def gen_bin_dec_genotipo_dinamico(individuos, gene_dec):
    gen_dec = np.zeros(individuos * gene_dec)
    gen_dec.shape = (individuos, gene_dec)
    return gen_dec


def genotipo_dinamico_string(individuos, gene_dec):
    gen_string = np.zeros(individuos * (gene_dec+1))
    return (gen_string)


def imprime_genbin(gen, individuos, gene, file):
    i, j = 0, 0
    print("\n")
    for i in range(0, individuos):
        for j in range(0, gene):
            # fprintf(pFile,"encerrado")
            print(gen[i][j])
    return gen


def imprime_genbindec(gen_bin_dec, individuos, gene_dec):
    i, j = 0, 0
    print("\n")
    for i in range(0, individuos):
        print("-Individuo["+str(i+1)+"]-\n")

        for j in range(0, gene_dec):
            print("%d\t", gen_bin_dec[i][j])
    print("\n\n\n")
    return gen_bin_dec


def imprime_genstring(gen_string, individuos, gene_dec, file):

    for i in range(0, individuos):
        print("\n-Individuo["+str(i+1)+"]-\n")
        # fprintf(pFile, "\n-Individuo[%d]-\n", i + 1);
        for j in range(0, gene_dec):
            print(gen_string[i][j])
            # fprintf(pFile, "%c", gen_string[i][j])
        CONTID += 1
        if gen_string[i][gene_dec + 1] is 'V':

            print("\t STRING VALIDA")
            # fprintf(pFile, "\t STRING VALIDA");
            mapeamento_genotipo_fenotipo(NENT, NSAI, aleatorio, TIPO, file)
        else:
            FIT[INDIC][0] = (INDIC + 1)
            FIT[INDIC][1] = 0.0
            FIT[INDIC][2] = FIT[INDIC][0]
            aux = indic + 1
            FIT[INDIC][3] = (INDIC + 1)
            if INDIC >= 1:
                FIT[INDIC][3] = INDIC + 1 + FIT[INDIC - 1][3]
            else:
                FIT[INDIC][3] = INDIC + 1
            FIT[INDIC][4] = 0.0
            FIT[INDIC][5] = 0.0
            FIT[INDIC][6] = 0.0
            INDIC += 1
            printf("\t STRING INVALIDA")
            # fprintf(pFile, "\t STRING INVALIDA"); }
    return FIT


def legenes_genbin(gen, gen_bin_dec, individuos, gene):

  # gene1, gene2, gene3, gene4, gene5, gene6, i;
    start = 0
    j = 0
    compactador = 0
    stop = (gene / 6)

    while (j < individuos):
        for i in range(0, stop):
            gene1 = gen[j][start]
            gene2 = gen[j][start + 1]
            gene3 = gen[j][start + 2]
            gene4 = gen[j][start + 3]
            gene5 = gen[j][start + 4]
            gene6 = gen[j][start + 5]
            decimal = converte_genbindec(gene1, gene2, gene3, gene4, gene5, gene6, j, compactador)
            gen_bin_dec[j][compactador] = decimal
            compactador+=1
            start = start + 6
        j+=1
        start = 0
        compactador = 0
  


def mapeamento_genotipo_fenotipo(NENT, NSAI, aleatorio,  tipo, ile):
    pass

def converte_genbindec(gene1, gene2, gene3, gene4, gene5, gene6, j, compactador):
    decimal = gene6 * 32 + gene5 * 16 + gene4 * 8 + gene3 * 4 + gene2 * 2 + gene1
    return decimal

def legenes_genbindec_string(gen_bin_dec, gen_string, individuos, gene_dec):
    for j in range(0, individuos):
        for i in range(0, gene_dec):
            gen_string[j][i] = tab_converte[gen_bin_dec[j][i]]

    return gen_string

