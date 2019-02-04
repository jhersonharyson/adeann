# from keras.models import Sequential 
# from keras.layers import Dense

# model = Sequential()
# model.add(Dense(units=64,activation='relu', input_dim=2))
# model.add(Dense(units=10, activation='softmax'))

# model.compile(loss='categorical_crossentropy',
#               optimizer='sgd',
#               metrics=['accuracy'])

# # model.fit(x_train, y_train, epochs=5, batch_size=32)

class AG:
    __init__(self):

        # ALGORITMO GENÉTICO
        self.GERACAO    = 1
        self.INDIVIDUOS = 60
        self.GENE       = 516
        self.NINTER     = 10000

        # REDE NEURAL
        self.TIPO = 4
        self.N = 4
        self.NINT = 4
        self.NINT1 = 7
        self.NINT2 = 7
        self.NENT = 5
        self.NSAI = 1
        self.NAPD = 900
        self.TARP3 = 0.01
        self.TARP4 = 0.01
        self.MAXINTER = 5000000
        self.ACEITAVEL = 0.001
        self.CONTID = 0
        self.CONTREGRASVAL = 0
        self.INDIC = 0
        self.SUP = 0
        self.NEURIT = 0
        
        # TAXAS 
        self.TX_TRANS = 0.001
        self.ERR = 0.0
        self.EMQ = 0.0
        self.FITNESS = 0.0
        self.AUX = 0
        self.SOMA_NINT_TOTAL = 0

        # CONSTANTES
        self.FIT = [[]]
        self.HIST_FIT = [[]]
    
    def zerar_fitness(fit):
        i, j = 0, 0
        for i in enumerate(range(0, self.INDIVIDUOS)):
             for j in enumerate(range(0, self.GENE)):
                fit[i][j] = 0

    def imprime_hist_fitness(hist_fit, file, contador):
        # imprimir no txt
        print("\n", contador) 
        print("\t\t", hist_fit[contador-1][0])
        print("\t\t", hist_fit[contador-1][1])

