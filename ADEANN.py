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
    i, j = 0, 0
    fit = np.zeros(INDIVIDUOS * GENE)
    fit.shape = (INDIVIDUOS, GENE)
    return fit


def imprime_hist_fitness(hist_fit, file, contador):
    # imprimir no txt
    print("\n", contador)
    print("\t\t", hist_fit[contador-1][0])
    print("\t\t", hist_fit[contador-1][1])


def genotipo_estatico(individuos, gene):
    i, j = 0, 0
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

