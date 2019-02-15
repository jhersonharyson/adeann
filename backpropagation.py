import numpy as np
from random import randint
import math


class backpropagation:
    def __init__(self):
        self.NINT = 5

        self.NINTER = 50
        self.NENT = 3
        self.NSAI = 1
        self.NPAD = 4
        self.TAPR = 0.9
        MAXITER = 100000
        ACEITAVEL = 0.0001
        err = 0.0
        emq = 0.0

        self.pFile = open("relat.txt", "w")
        self.x = [] #matriz numero de padroes com entrada
        self.h = [None] * self.NINTER #NINTER
        self.o = [0.0] * self.NSAI #NSAI
        self.y = [] #matriz numero de padroes com saida
        self.delta1 = [0.0] * self.NINTER #NINTER
        self.delta2 = [0.0] * self.NSAI #NSAI
        self.w1 = [] #entrada p intermediaria
        self.w2 = [] #intermediaria p saida

        #system("cls");
        self.pFile.write("Aprendizagem da Rede\n")
        print("Aprendizagem da Rede\n")
        self.conjunto_treinamento1()
        self.inicializa_pesos()
        erromax = ACEITAVEL * 2.0
        m = 0
        l = 0

        while (erromax > ACEITAVEL and l < MAXITER):
            if(m == self.NPAD):
                m = 0
            self.intermediaria(m)
            self.saida()
            err = self.erro_saida(m)

            if (err > erromax):
                erromax = err

            self.erro2(m)
            self.erro1()
            self.ajusta2()
            self.ajusta1(m)

            l = l + 1
            print("\nPadrao>>" + str(m), end="")
            m = m + 1
            print("\t\t Epoca:"+ str(l) +"\t\t Erro:%.9f" % err)

        self.pFile.write("\n\nVERIFICACAO\n")
        print("\n\nVERIFICACAO\n")
        self.verifica1()
        self.pFile.close()

    def conjunto_treinamento1(self):
        self.y = [
            [0],
            [0],
            [0],
            [1]
        ]
        self.x = np.array([
            [1.0, 0, 0],
            [1.0, 0, 1],
            [1.0, 1, 0],
            [1.0, 1, 1]
        ])

    def inicializa_pesos(self):
        self.w1 = np.zeros(self.NENT * self.NINT)
        self.w1.shape = (self.NENT, self.NINT)

        self.w2 = np.zeros(self.NINT * self.NSAI)
        self.w2.shape = (self.NINT, self.NSAI)


        for i in range(0, self.NENT):
           for j in range(1, self.NINT):
                aleatorio = (randint(0, 100))
                self.w1[i][j] = (1.0 - 2.0 * aleatorio) / 100.0

        for j in range(0, self.NINT):
            for k in range(0, self.NSAI):
                aleatorio = (randint(0, 100))
                self.w2[j][k] = (1.0 - 2.0 * aleatorio) / 100.0

    def intermediaria(self, m):
        self.h[0] = 1.0
        for j in range(1, self.NINT):
            somatorio = 0.0

            for i in range(0, self.NENT):
                somatorio = somatorio + self.x[m][i] * self.w1[i][j]

            somatorio = - somatorio
            self.h[j] = 1.0 / (1.0 + math.exp(somatorio))


    def saida(self):
        for k in range(0, self.NSAI):
            somatorio = 0.0

            for j in range(0, self.NINT):
                somatorio = somatorio + self.h[j] * self.w2[j][k]

            somatorio = - somatorio
            self.o[k] = 1.0 / (1.0 + math.exp(somatorio))


    def erro_saida(self, m):
        somatorio = 0.0

        for k in range(0, self.NSAI):
            somatorio = somatorio + (self.o[k] - self.y[m][k]) * (self.o[k] - self.y[m][k])
        erro = 0.5 * somatorio

        return erro



    def erro2(self, m):
        for k in range(0, self.NSAI):
            self.delta2[k] = self.o[k] * (1.0 - self.o[k]) * (self.y[m][k] - self.o[k])


    def erro1(self):
        for j in range(1, self.NINT):
            somatorio = 0.0

            for k in range(0, self.NSAI):
                somatorio = somatorio + self.delta2[k] * self.w2[j][k]

                self.delta1[j] = self.h[j] * (1 - self.h[j]) * somatorio


    def ajusta2(self):
        for j in range(0, self.NINT):
            for k in range(0, self.NSAI):
                self.w2[j][k] = self.w2[j][k] + self.TAPR * self.delta2[k] * self.h[j]


    def ajusta1(self, m):
        for i in range(0, self.NENT):
            for j in range(1, self.NINT):
                self.w1[i][j] = self.w1[i][j] + self.TAPR * self.delta1[j] * self.x[m][i]


    def verifica1(self):
        acerto = 0
        cont = 0
        err = 0.0
        emq = 0.0

        for m in range(0, self.NPAD):
            self.intermediaria(m)
            self.saida()
            err = self.erro_saida(m)
            self.pFile.write("\nPadrao Generalizacao>>"+ str(m))
            self.pFile.write("\ncalculado>>"+ str(self.o[0]) +"  \tdesejado>>"+ str(self.y[m][0]) +"  \tErro>>%.9f" % err)
            print("\nPadrao Generalizacao>>"+ str(m))
            print("\ncalculado>>"+ str(self.o[0]) +"  \tdesejado>>"+ str(self.y[m][0]) +"  \tErro>>%.9f" % err)
            emq = emq + err

            if (err <= 0.0005):
                acerto += 1

        cont = self.NPAD - acerto
        emq = emq / self.NPAD

        print("\nemq>>"+ str(emq))
        self.pFile.write("\n\nErro Medio Quadratico>> %.9f" % emq)
        self.pFile.write("\n\n<<Pesos Camada Entrada Oculta>>")

        for i in range(0, self.NENT):
            for j in range(1, self.NINT):
                self.pFile.write( "\n")
                self.pFile.write( " w1["+ str(i) +"]["+ str(j) +"]="+ str(self.w1[i][j]))

            self.pFile.write("\n")

        self.pFile.write("\n\n<<Pesos Camada Oculta Saida>>")
        for i in range(0, self.NINT):
            for j in range(0, self.NSAI):
                self.pFile.write("\n")
                self.pFile.write(" w2["+ str(i) +"]["+ str(j) +"]="+ str(self.w2[i][j]))

            self.pFile.write("\n")

        for m in range(0, self.NPAD):
            self.intermediaria(m)
            self.saida()
            err = self.erro_saida(m)
            self.pFile.write("\nPadrao Testes>>"+ str(m))
            self.pFile.write("\ncalculado>>"+ str(self.o[0]) +"   \tdesejado>>"+ str(self.y[m][0]) +"  \tErro>>%.9f" % err)
            print("\nPadrao Testes>>"+ str(m))
            print("\ncalculado>>"+ str(self.o[0]) +"   \tdesejado>>"+ str(self.y[m][0]) +"  \tErro>>%.9f" % err)
            emq = emq + err

        self.pFile.write("\n Contagem de erro:>>"+ str(cont))
        self.pFile.write("\n Contagem de acertos:>>"+ str(acerto))
        porcAcerto = acerto / self.NPAD * 100.0
        self.pFile.write("\n Porcentagem de acertos:>> %.2f" % porcAcerto)

NN = backpropagation()