from keras.models import Sequential
from keras.layers import Dense
import numpy as np

# ALGORITMO GENÉTICO
GERACAO = 10
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
MAXITER = 5000000
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
HIST_FIT = [(None, [])] * GERACAO

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
    print("\t\t", hist_fit[contador - 1][0])
    print("\t\t", hist_fit[contador - 1][1])


def genotipo_estatico(individuos, gene):
    gen = np.zeros(individuos * gene)
    gen.shape = (individuos, gene)
    for i in range(0, gene):
        for j in range(0, individuos):
            gen[i][j] = np.random.random_integers(0, 1, 1)[0]
    return gen


def genotipo_dinamico(individuos, gene):
    print(gene, " - ", int(gene))
    gene = int(gene)
    gen = np.zeros(individuos * gene)
    gen.shape = (individuos, gene)

    for i in range(0, individuos):
        for j in range(0, gene):
            gen[i][j] = np.random.random_integers(0, 1, 1)[0]

    return gen


def gen_bin_dec_genotipo_dinamico(individuos, gene_dec):
    gen_dec = np.zeros(individuos * gene_dec)
    gen_dec.shape = (individuos, gene_dec)
    return gen_dec


def genotipo_dinamico_string(individuos, gene_dec):
    print(gene_dec)
    gene_dec = int(gene_dec)
    gen_string = np.chararray(individuos * (gene_dec + 1))
    gen_string.shape = (individuos, gene_dec + 1)
    return gen_string


def imprime_genbin(gen, individuos, gene, file):
    print("\n")
    for i in range(0, individuos):
        for j in range(0, gene):
            # fprintf(pFile,"encerrado")
            print(gen[i][j])
    return gen


def imprime_genbindec(gen_bin_dec, individuos, gene_dec):
    print("\n")
    gene_dec = int(gene_dec)
    for i in range(0, individuos):
        print("-Individuo[" + str(i + 1) + "]-\n")

        for j in range(0, gene_dec):
            print("\t", gen_bin_dec[i][j])
    print("\n\n\n")
    return gen_bin_dec


def imprime_genstring(gen_string, individuos, gene_dec, file):
    global CONTID, INDIC, FIT
    for i in range(0, individuos):
        print("\n-Individuo[" + str(i + 1) + "]-\n")
        # fprintf(pFile, "\n-Individuo[%d]-\n", i + 1);
        for j in range(0, gene_dec):
            print(gen_string[i][j])
            # fprintf(pFile, "%c", gen_string[i][j])
        CONTID += 1
        if gen_string[i][gene_dec] is b'V':

            print("\t STRING VALIDA")
            # fprintf(pFile, "\t STRING VALIDA");
            mapeamento_genotipo_fenotipo(NENT, NSAI, 0, TIPO, file)  # 0 -> aleatorio
        else:
            FIT[INDIC][0] = INDIC + 1
            FIT[INDIC][1] = 0.0
            FIT[INDIC][2] = FIT[INDIC][0]
            FIT[INDIC][3] = INDIC + 1
            if INDIC >= 1:
                FIT[INDIC][3] = INDIC + 1 + FIT[INDIC - 1][3]
            else:
                FIT[INDIC][3] = INDIC + 1
            FIT[INDIC][4] = 0.0
            FIT[INDIC][5] = 0.0
            FIT[INDIC][6] = 0.0
            INDIC += 1
            print("\t STRING INVALIDA")
            # fprintf(pFile, "\t STRING INVALIDA"); }
    return FIT


def legenes_genbin(gen, gen_bin_dec, individuos, gene):
    # gene1, gene2, gene3, gene4, gene5, gene6, i;
    start = 0
    j = 0
    compactador = 0
    stop = int(gene / 6)

    while j < individuos:
        for i in range(0, stop):
            gene1 = gen[j][start]
            gene2 = gen[j][start + 1]
            gene3 = gen[j][start + 2]
            gene4 = gen[j][start + 3]
            gene5 = gen[j][start + 4]
            gene6 = gen[j][start + 5]
            decimal = converte_genbindec(gene1, gene2, gene3, gene4, gene5, gene6, j, compactador)
            gen_bin_dec[j][compactador] = decimal
            compactador += 1
            start = start + 6
        j += 1
        start = 0
        compactador = 0


def converte_genbindec(gene1, gene2, gene3, gene4, gene5, gene6, j, compactador):
    decimal = gene6 * 32 + gene5 * 16 + gene4 * 8 + gene3 * 4 + gene2 * 2 + gene1
    return decimal


def legenes_genbindec_string(gen_bin_dec, gen_string, individuos, gene_dec):
    for j in range(0, individuos):
        for i in range(0, gene_dec):
            gen_string[j][i] = tab_converte[int(gen_bin_dec[j][i])]

    return gen_string


def avalia_regras_gen_string(gen_string, individuos, gene_dec):
    string_val = [b'.', b'f', b'[', b'F', b'f', b'n', b'B', b']']
    for j in range(0, individuos-1):
        for i in range(0, gene_dec - 7):
            # Le do inicio para o final com um passo de um bit.
            if (gen_string[j][i] is b'.') and (gen_string[j][i + 1] is b'f') and (gen_string[j][i + 2] is b'[') and (
                    gen_string[j][i + 3] is b'F') and (gen_string[j][i + 4] is b'f') and (
                    gen_string[j][i + 5] is b'n') and (gen_string[j][i + 6] is b'B') and (gen_string[j][i + 7] is b']'):
                gen_string[j][gene_dec] = b'V'

    for j in range(0, individuos-1):
        # Le do final para o inicio com um passo de um bit.
        for i in range(gene_dec-1, 6, -1):
            if (gen_string[j][i] is b'.') and (gen_string[j][i + 1] is b'f') and (gen_string[j][i + 2] is b'[') and (
                    gen_string[j][i + 3] is b'F') and (gen_string[j][i + 4] is b'f') and (
                    gen_string[j][i + 5] is b'n') and (gen_string[j][i + 6] is b'B') and (gen_string[j][i + 7] is b']'):
                gen_string[j][gene_dec] = b'V'

    # mudei aqui baixo todo la?o for

    for j in range(0, individuos-1):
        # Le de forma sequencial bit a bit do inicio para o fim.
        cont_elem = 0
        pos_string = 0

        for i in range(0, gene_dec-1):
            carac = gen_string[j][i]

            if not (cont_elem >= 8):
                if string_val[pos_string] is carac:
                    cont_elem += 1
                    pos_string += 1

            if cont_elem is 8:
                gen_string[j][gene_dec] = b'V'
                # //printf("j=%d",j);
                # //printf("cont_elem=%d",cont_elem);
                # //getch();
    #print("\n2222222222222\n")
    #print(gen_string)
    return gen_string


def mapeamento_genotipo_fenotipo(NENT, NSAI, aleatorio, TIPO, file):
    global NINT, NINT1, NINT2
    # if NR is NR1:
    #     NR1 = np.random.random_integers(2, 20)
    N_TIPO = np.random.random_integers(2, 20, 1)[0]

    NINT = N_TIPO
    NINT1 = N_TIPO
    NINT2 = np.random.random_integers(2, 20, 1)[0]

    NINT_N3 = np.zeros(1)
    NINT_N3[0] = NINT1
    SIZE_N3 = len(NINT_N3)

    NINT_N4 = np.zeros(2)
    NINT_N4[0] = NINT1
    NINT_N4[1] = NINT2
    SIZE_N4 = len(NINT_N4)
    ############################################################################################
    # n == 3 ? Mapeamento(NENT, NSAI, NINT_N3, SIZE_N3, "1.4", pFile) : Mapeamento(NENT, NSAI, NINT_N4, SIZE_N4, "1.4", pFile); // (ENTRADA, SAIDA, NR, TIPO)

    if N is 3:
        treina_rede(CONTID, file, NINT)
    if N is 4:
        treina_rede_(CONTID, file, NINT1, NINT2)


def NR_RAND(int):
    pass


def treina_rede(individuos, file, NINT):
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
    print("Treina rede")


def treina_rede_(contind, pFile, NINT1, NINT2):
    print("\nNENT=(" + str(NENT - 1) + " Entradas + 1 Bias)")
    print("\nNINT1=(" + str(NINT1 - 1) + " Int + 1 Bias)")
    print("\nNINT2=(" + str(NINT2 - 1) + " Int + 1 Bias)")
    print("\nNSAI=", NSAI)
    print("\n\nTreinamento do Individuo=", INDIVIDUOS)

    x_train = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    y_train = np.array([[0], [1], [1], [0]])

    model = Sequential()
    model.add(Dense(units=5, activation='relu', input_dim=2))
    model.add(Dense(units=5, activation='softmax'))
    model.add(Dense(units=1, activation='softmax'))

    model.compile(loss='binary_crossentropy',
                  optimizer='sgd',
                  metrics=['accuracy'])

    model.fit(x_train, y_train, epochs=5, batch_size=32)
    predictions = model.predict(x_train)

    rounded = [x for x in predictions]

    for i, predic in enumerate(rounded):
        print("\nPadrao>>", i)
        print("\ncalculado>>" + str(predic) + "   \tdesejado>>" + str(y_train[i]) + "  \tErro>>", y_train[i] - predic)
    print("\n\nTreina rede 4")


def zera_fitness(fit):
    fit = np.zeros(INDIVIDUOS * GENE)
    fit.shape = (INDIVIDUOS, GENE)
    return fit


def ordena(fit, file):
    min = 0
    for i in range(0, INDIVIDUOS - 1):
        min = i
        for j in range(i + 1, INDIVIDUOS):
            AUX1 = fit[j][1]
            AUX2 = fit[min][1]
            if AUX1 > AUX2:
                min = j

        ch = fit[i][0]
        ch1 = fit[i][1]
        ch2 = fit[i][4]
        ch3 = fit[i][5]
        ch4 = fit[i][6]
        ch5 = fit[min][0]
        fit[i][0] = ch5
        fit[i][1] = fit[min][1]
        fit[i][4] = fit[min][4]
        fit[i][5] = fit[min][5]
        fit[i][6] = fit[min][6]
        fit[min][0] = ch
        fit[min][1] = ch1
        fit[min][4] = ch2
        fit[min][5] = ch3
        fit[min][6] = ch4


def main():
    file = ''
    global FIT, HIST_FIT, MAXITER
    # pFile = fopen("relat_autoMaasd2asdasd.txt", "w");
    for n in range(4, 4):

        if (n is 3):
            pass
            # pFile = fopen("relat_N3asdas.txt", "w");
        else:
            pass
            # fprintf(pFile, "\nTARP1 %.2f\nTARP2 %.2f\n\n\n");

    # gene : represneta o numero de genes por inidividuo, nesse caso 516;
    gene_dec = int(GENE / 6)  # //43+1=44 ultimo elemento armazena status da string valida % = valida;

    # //////    printf("Gerando Genotipo Aleatoriamente!\n");
    gen_bin = genotipo_dinamico(INDIVIDUOS, GENE)
    gen_bin_dec = genotipo_dinamico(INDIVIDUOS, gene_dec)
    gen_string = genotipo_dinamico_string(INDIVIDUOS, gene_dec)

    FIT = zera_fitness(FIT)
    contador = 1  # //contador do numero de geracoes
    global contador1
    contador1 = 0
    # //int pontos_corte=(gene_dec*0.9);

    while contador <= GERACAO:

        imprime_genbin(gen_bin, INDIVIDUOS, GENE, file)
        legenes_genbin(gen_bin, gen_bin_dec, INDIVIDUOS, GENE)
        imprime_genbindec(gen_bin_dec, INDIVIDUOS, gene_dec)
        legenes_genbindec_string(gen_bin_dec, gen_string, INDIVIDUOS, gene_dec)
        contador1 = 1  # //contador do sumero de cruzamentos por gera??o
        FIT = zera_fitness(FIT)
        gen_string = avalia_regras_gen_string(gen_string, INDIVIDUOS, gene_dec)  # //Avalia Regras Validas
        imprime_genstring(gen_string, INDIVIDUOS, gene_dec, file)  # //Imprime Individuo[i] e DNA
        # //Ordenar Individuos
        ordena(FIT, file)
        # //Imprimir Individuos

        imprime_fitness(FIT, file, contador)
        contador += 1
        INDIC = 0
        while contador1 <= (gene_dec * 0.8):
            # realiza (individuos/2) cruzamentos
            gen_bin = selecao(gen_bin, gen_string, gene_dec)
            contador1 += 1

        contador1 = 1
        while contador1 <= (gene_dec * 0.8):
            # //516-86-14
            gen_bin = mutacao(gen_bin)
            contador1 += 1
        # // transloca(gen_string,gene_dec);

        if GERACAO < 20:
            MAXITER = MAXITER + 25000
        elif (GERACAO >= 20) and (GERACAO < 40):
            MAXITER = MAXITER + 80000
        else:
            MAXITER = MAXITER + 150000

        contador = 1

        imprime_cabec(file)
        while GERACAO >= contador:
            imprime_hist_fitness(HIST_FIT, file, contador)
            contador += 1

        # //End Do

        print("\n\n<<Simulacao Concluida - Relatorio Gerado!!>>")
        # fprintf(pFile, "\nsimulacao concluida");


def selecao(gen, gen_string, gene_dec):
    global e5
    e1 = 0
    e = 0
    j = INDIVIDUOS - 1

    while e == e1:
        e = np.random.random_integers(0, 1, 1)[0]
        e2 = int(FIT[e][1])

        e2 -= 1
        if e2 >= (INDIVIDUOS / 2):
            e2 = 1
        e1 = np.random.random_integers(0, 1, 1)[0]
        e3 = int(FIT[e1][1])
        e3 -= 1
        if e3 >= (INDIVIDUOS / 2):
            e3 = 0

    e4 = np.random.random_integers(0, 1, 1)[0]
    e5 = int(FIT[e4][1])
    e5 -= 1
    if e5 <= (INDIVIDUOS / 2) - 1:
        e5 = j - e5

    return cruzamento(e2, e3, e3, gen, gen_string, gene_dec)


def cruzamento(i, j, pos, gen, gen_string, gene_dec):
    pc = np.random.random_integers(0, INDIVIDUOS, 1)[0]
    pc = (pc / 100)
    ale = (GENE - 12)
    gen_ale = np.random.random_integers(0, ale, 1)[0]
    gen_ale1 = np.random.random_integers(0, ale, 1)[0]
    if pc <= 0.6:
        gen[pos][gen_ale] = gen[i][gen_ale]
        gen[pos][gen_ale1] = gen[j][gen_ale1]
        gen[pos][gen_ale + 1] = gen[i][gen_ale + 1]
        gen[pos][gen_ale1 + 1] = gen[j][gen_ale1 + 1]
        # //
        gen[pos][gen_ale + 2] = gen[i][gen_ale + 2]
        gen[pos][gen_ale1 + 2] = gen[j][gen_ale1 + 2]
        # //
        gen[pos][gen_ale + 3] = gen[i][gen_ale + 3]
        gen[pos][gen_ale1 + 3] = gen[j][gen_ale1 + 3]
        # //
        gen[pos][gen_ale + 4] = gen[i][gen_ale + 4]
        gen[pos][gen_ale1 + 4] = gen[j][gen_ale1 + 4]
        # //
        gen[pos][gen_ale + 5] = gen[i][gen_ale + 5]
        gen[pos][gen_ale1 + 5] = gen[j][gen_ale1 + 5]
    return gen


def mutacao(gen):
    pm = np.random.random_integers(0, 10, 1)[0]
    pm = (pm / 100)

    e6 = np.random.random_integers(0, 1, 1)[0]
    e7 = int(FIT[e6][0])
    e7 -= 1
    if e7 <= (INDIVIDUOS / 2) - 1:
        e7 = (INDIVIDUOS - 1) - e6

    ale = (GENE - 1)
    gen_ale = np.random.random_integers(0, ale, 1)[0]

    if (gen[e7][gen_ale] is 0) and (pm <= 0.1):
        gen[e7][gen_ale] = 1
    else:
        gen[e7][gen_ale] = 0
    return gen


def imprime_fitness(fit, file, contador):
    global HIST_FIT, SOMA_NINT_TOTAL
    soma_fit = 0.0
    soma_nint = 0.0
    contador_val = 0
    var_nint = 0.0

    # fprintf(pFile, "\n\nResumo da Geracao:%d", contador);
    # fprintf(pFile, "\n_______________________________________________________________");
    # n == 3 ? fprintf(pFile, "\nIND  |Fitness      |Posto  |Acum    |EMQ      |NINT  |NT") : fprintf(pFile, "\nIND  |Fitness      |Posto  |Acum    |EMQ      |NINT1 |NINT2  |NT");
    # fprintf(pFile, "\n_______________________________________________________________");
    print("\n\nResumo da Geracao:%d", contador)
    print("\n_______________________________________________________________")
    print(
        "\nIND  |Fitness      |Posto  |Acum    |EMQ      |NINT  |NT" if N == 3 else "\nIND  |Fitness      |Posto  |Acum    |EMQ      |NINT1 |NINT2  |NT")
    print("\n_______________________________________________________________")

    for i in range(0, INDIVIDUOS):

        print("\n", fit[i][0])
        print("    ", fit[i][1])
        print("     ", fit[i][2])
        print("     ", fit[i][3])
        print("      ", fit[i][4])
        print("   ", fit[i][5])
        # fprintf(pFile, "\n%3.0f", fit[i][0]);
        # fprintf(pFile, "    %3.4f", fit[i][1]);
        # fprintf(pFile, "     %3.0f", fit[i][2]);
        # fprintf(pFile, "     %3.0f", fit[i][3]);
        # fprintf(pFile, "      %3.6f", fit[i][4]);
        # fprintf(pFile, "   %3.0f", fit[i][5]);
        if (N == 4):
            print("      ", fit[i][7])
            # fprintf(pFile, "      %d", fit[i][7]);
        print("   ", fit[i][6])
        # fprintf(pFile, "   %3.0f", fit[i][6]);
        soma_fit = soma_fit + fit[i][1]
        if 0 < fit[i][5]:
            soma_nint = soma_nint + fit[i][5]
            contador_val += 1
        global SUP
        SUP = fit[i][3]

    print("\n===============================================================\n")
    # fprintf(pFile, "\n===============================================================\n");
    if 0 < contador_val:
        SOMA_NINT_TOTAL = SOMA_NINT_TOTAL + (soma_nint / contador_val)
        print("\nMedia do Fitness=\n", soma_fit / contador_val)
        print("Media de Neuronios na Camada Intermediaria na Geracao=\n", soma_nint / contador_val)
        # fprintf(pFile, "\nMedia do Fitness=%3.4f\n", soma_fit / contador_val);
        # fprintf(pFile, "Media de Neuronios na Camada Intermediaria na Geracao=%3.2f\n", soma_nint / contador_val);
    # [(None, []), (None, []), (None, []), (None, []), (None, []), (None, []), (None, []), (None, []), (None, []), (None, [])]
    HIST_FIT[contador - 1][1].append(fit[0][1])

    # HIST_FIT[contador - 1][0] = fit[0][1];

    if 0 < contador_val:
        # HIST_FIT[contador - 1][1] = soma_fit / contador_val;
        HIST_FIT[contador - 1][1].append(soma_fit / contador_val)
    else:
        # HIST_FIT[contador - 1][1] = 0.0;
        HIST_FIT[contador - 1][1].append(0.0)

    for i in range(0, INDIVIDUOS):

        if 0 < fit[i][5]:
            var_nint = var_nint + pow((fit[i][5] - (soma_nint / contador_val)), 2) / contador_val

    print("Variancia de Neuronios na Camada Intermediaria na Geracao=\n", var_nint)
    print("Desvio Padrao de Neuronios na Camada Intermediaria na Geracao=\n", (var_nint) ** (1 / 2))
    # fprintf(pFile, "Variancia de Neuronios na Camada Intermediaria na Geracao=%3.2f\n", var_nint);
    # fprintf(pFile, "Desvio Padrao de Neuronios na Camada Intermediaria na Geracao=%3.2f\n", sqrt(var_nint));
    if 0 < contador_val:
        print("Erro Padrao de Neuronios na Camada Intermediaria na Geracao=\n",
              var_nint ** (-1 / 2) / contador_val ** (1 / 2))
        print("Coeficiente de Variacao de Neuronios na Camada Intermediaria na Geracao=\n",
              var_nint ** (1 / 2) / (soma_nint / contador_val))

        # fprintf(pFile, "Erro Padrao de Neuronios na Camada Intermediaria na Geracao=%3.2f\n", sqrt(var_nint) / sqrt(contador_val));
        # fprintf(pFile, "Coeficiente de Variacao de Neuronios na Camada Intermediaria na Geracao=%3.2f\n", sqrt(var_nint) / (soma_nint / contador_val));

    if contador == GERACAO:
        # fprintf(pFile, "\n\nMedia de Neuronios na Camada Intermediaria na Simulacao=%3.2f\n", (SOMA_NINT_TOTAL / contador));
        print("\n\nMedia de Neuronios na Camada Intermediaria na Simulacao=\n", (SOMA_NINT_TOTAL / contador))


def imprime_cabec(file):
    print("\nPercentagem de Regras Validas=\n %", (CONTREGRASVAL / (INDIVIDUOS * GERACAO)) * 100.0)
    print("\n\nHistorico do Fitness na Simulacao")
    print("\n________________________________________________________________")
    print("\n Geracao  |Melhor Fitness         |Fitness Medio")
    print("\n________________________________________________________________")


main()
