#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int static NENT = 2; // Numero de Entrada
int static NSAI = 1; // Numero de Saída
int NINT = 1;
int NINT1 = 1;
int NINT2 = 1;
FILE * pFile;
int q = 0;
int n = 4;
int tipo = 4;


int NR_RAND()
{ // Neuronios Por Ramifica??o 1-10
  srand((unsigned)time(NULL));
  int x = 2 + (rand() % 8);
  return x;
}
int NR_RAND_N(int n)
{ // Neuronios Por Ramifica??o 1-10
  srand((unsigned)time(NULL));
  int x = n + (rand() % (10 - n));
  return x;
}

int main(int argc, char *argv[ ]){
	
  	if(argc > 1){
	
  	//for(cont=0; cont < argc; cont++){
//    printf("NOME PARAMETRO: %s\n", argv[0]);
//    printf("NENT PARAMETRO: %s\n", argv[1]);
//    printf("NSAI PARAMETRO: %s\n", argv[2]);
//    printf("N    PARAMETRO: %s\n", argv[3]);
//    if(strcmp(argv[3],"3") == 0){
//    	printf("NENT PARAMETRO: %s\n", argv[4]);
//	}else if(strcmp(argv[3],"4") == 0){
//		printf("NENT1 PARAMETRO: %s\n", argv[4]);
//		printf("NENT2 PARAMETRO: %s\n", argv[5]);
//	}

	//}
	//pFile = fopen("resultado.txt", "w");
	NENT = atoi(argv[1]);
	NSAI = atoi(argv[2]);	
	///////////////////////////////////////////////////////////////////////////////////
  	int NR = atoi(argv[4]);
  	int NR1 = atoi(argv[5]);
	n = atoi(argv[3]);
  	NINT = NR;
  	NINT1 = NR;
  	NINT2 = NR1;

  	int NINT_N3[1];
  	NINT_N3[0] = NINT1;
  	int SIZE_N3 = sizeof(NINT_N3) / sizeof(NINT_N3[0]);

  	int NINT_N4[2];
  	NINT_N4[0] = NINT1;
  	NINT_N4[1] = NINT2;
  	int SIZE_N4 = sizeof(NINT_N4) / sizeof(NINT_N4[0]);
	///////////////////////////////////////////////////////////////////////////////////
	//int SIZE = sizeof(NINT)/sizeof(NINT[0]);
    n == 3 ? Mapeamento(NENT, NSAI, NINT_N3, SIZE_N3, "1.4") : Mapeamento(NENT, NSAI, NINT_N4, SIZE_N4, "1.4"); // (ENTRADA, SAIDA, NR, TIPO)
}
return 0;
}// fim main

// MAPEAMENTO
int Mapeamento(int NENT, int NSAI, int NR[], int SIZE, char TIPO[])
{
  int q = 0;
  //FILE *txt = pFile; // apaga possiveis lixos no arquivo
  //fclose(txt);

  ///////////////////////////////////////////////
  ////////////////// STRUCT /////////////////////
  ///////////////////////////////////////////////
  struct ED
  { // estrutura Lista

    char ele;        // Dados
    struct ED *Prox; // Proximo
    struct ED *Ant;  // Anterior
  };

  struct ED *Primeiro[1000];
  struct ED *Ultimo[1000]; // Ultimo elemento
  /////////////////////////////////////////////////
  /////////////////////////////////////////////////

  int size = 0;

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////
  ///////////////// METODOS ED /////////////////
  //////////////////////////////////////////////

  int pList = 0;

  void QueueNull()
  { // inicaliza fila NULL

    struct ED *aux;

    aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED)); // aoca??o de memoria
    Primeiro[q] = aux;
    Ultimo[q] = Primeiro[q];
    Primeiro[q]->Ant = NULL;
  }

  void InsertInFronte(char elemento)
  { // insere dentro da ED no inicio

    struct ED *aux;
    aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED)); // aloca??o de memoria

    aux->ele = elemento;
    aux->Ant = NULL;
    aux->Prox = NULL;

    if (!Primeiro[q]->ele)
    {
      //Primeiro[q] = Primeiro[q]->Prox;
      Primeiro[q] = Ultimo[q] = aux;
    }

    else
    {
      aux->Prox = Primeiro[q];
      Primeiro[q]->Ant = aux;
      Primeiro[q] = aux;
    }
    size++;
  }

  void ReplaceOnPosition(char elemento, int pos)
  {

    int p = 1;

    struct ED *aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;
    while (aux != NULL)
    {
      if (p == pos)
      {
        aux->ele = elemento;
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        p++;
      }
    }
    if (flag == 1)
    {
      //////	     printf("O Elemento %c N?o existe na ED\n", elemento);
    }
  }

  void InsertMiddle(char elemento, int posicao)
  { // insere em qualquer lugar
    char count = posicao;
    count--;
    struct ED *aux;
    struct ED *aux2;
    aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED));  //aloca??o de memoria para aux
    aux2 = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED)); //aloca??o de memoria para aux2

    if (count < 0)
    {
//      printf("\nn?o Existe Esta Posi??o");
    }
    else
    {
      if (count == 0)
      {
        aux->ele = elemento;
        aux->Prox = Primeiro[q];
        aux->Ant = NULL;
        Primeiro[q]->Ant = aux;

        aux = Primeiro[q];
        size++;
      }
      else
      {

        aux->Prox = Primeiro[q]->Prox;
        aux->Ant = Primeiro[q]->Ant;
        aux = Primeiro[q];
      }
      while (count != 0)
      {

        aux = aux->Prox;
        count--;
      }
      aux2->ele = elemento;

      aux->Ant->Prox = aux2;
      aux2->Prox = aux;
      aux2->Ant = aux->Ant;
      aux->Ant = aux2;
      size++;
    }
  }

  void InsertLast(char elemento)
  { // insere pelo final
    struct ED *aux;

    aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED));
    aux->ele = elemento;

    if (!Primeiro[q]->ele)
    {
      Primeiro[q] = Primeiro[q]->Prox;
      Primeiro[q] = Ultimo[q] = aux;
      size = 0;
      size++;
    }
    else
    {

      Ultimo[q]->Prox = aux;
      aux->Ant = Ultimo[q];
      aux->Prox = NULL;
      Ultimo[q] = aux;
      size++;
    }
  }

  void ShowStack()
  { // Imprime (PILHA)
    //FILE *save;
    //save = pFile;
    struct ED *aux;
    aux = Primeiro[q]; //->Prox;
    if (!aux->ele)
    {
      printf("A ED est? Vazia !");
    }
    else
    {
    	printf("  ||  ");
//      fprintf(pFile, "  ||  ");
      while (aux != NULL)
      {

        //printf("\nElemento: %c \n", aux->ele);
        	printf(" %c ",aux->ele);
        //fprintf(pFile, " %c ", aux->ele);
        aux = aux->Prox;
      }
    }
    if (q == NENT - 1)
    {
      	  	printf("  ||");
//      fprintf(pFile, "  ||");
      	   	printf("\n\n");
//      fprintf(pFile, "\n\n");
    }
  }
  
   int SearchLastF(char elemento)
  { // pesquisa por elemntos dentro da ED

    int pos = 1;
    struct ED *aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;

    while ((Primeiro[q]->ele == 'f') && (Primeiro[q]->Prox) && Primeiro[q]->ele == 'f')
    {

      aux = Primeiro[q]->Prox;
      pos++;
    }

    while (aux != NULL)
    {
      if ((aux->ele == 'f' && aux->Ant->ele == ']'))
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        pos++;
      }
    }

    if (flag == 1)
    {
      // printf("O Elemento %c N?o existe na ED\n\n", elemento);
      pos = 0;
    }

    return pos;
  }
  
  int Search(char elemento)
  { // pesquisa por elemntos dentro da ED

    int pos = 1;
    struct ED *aux = (struct ED *)calloc(sizeof(struct ED), sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;

    while ((Primeiro[q]->ele == 'f') && (Primeiro[q]->Prox) && Primeiro[q]->ele == 'f')
    {

      aux = Primeiro[q]->Prox;
      pos++;
    }

    while (aux != NULL)
    {
      if (aux->ele == elemento || (aux->ele == 'f' && aux->Ant->Ant->ele == ']'))
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        pos++;
      }
    }

    if (flag == 1)
    {
      // printf("O Elemento %c N?o existe na ED\n\n", elemento);
      pos = 0;
    }

    return pos;
  }

  int SearchR36(char elemento, int pos_start)
  { // pesquisa por elemntos dentro da ED
    int pos = 1;
    int c = 1;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;

    while (pos_start >= c && aux != NULL)
    {
      aux = aux->Prox;
      pos++;
      c++;
    }

    while (aux != NULL)
    {
      if (aux->ele == elemento && aux->Prox->ele == ']' && aux->Prox != NULL)
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        pos++;
        if (aux->ele == ']' && aux->Prox->ele == 'f')
        {
          flag = 1;
          aux = NULL;
          pos = 0;
        }
      }
    }

    if (flag == 1 || pos >= size)
    {
      //printf("O Elemento %c N?o existe na ED\n\n", elemento);
      pos = 0;
    }
    return pos;
  }


  int SearchR322(char elemento, int pos_start)
  { // pesquisa por elemntos dentro da ED
    int pos = 1;
    int c = 1;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;

    while (pos_start >= c && aux != NULL)
    {
      aux = aux->Prox;
      pos++;
      c++;
    }

    while (aux != NULL)
    {
      if (aux->ele == elemento && aux->Ant->Ant->ele == ']' || (aux->ele == elemento && aux->Ant->Ant->ele == '['))
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        pos++;
      }
    }

    if (flag == 1 || pos >= size)
    {
      ////////   printf("O Elemento %c N?o existe na ED\n\n", elemento);
      pos = 0;
    }
    return pos;
  }

  int SearchR341(char elemento, int pos_start)
  { // pesquisa por elemntos dentro da ED
    int pos = 1;
    int c = 1;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1, flag2 = 0;

    while (pos_start >= c && aux != NULL)
    {
      aux = aux->Prox;
      pos++;
      c++;
    }

    while (aux != NULL)
    {
      if (!aux->Prox)
      {
        flag = 1;
        break;
      }
      if ((aux->ele == elemento && aux->Prox->ele == 'B' && aux->ele != 'B'))
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        // printf("|%c|",aux->ele);
        flag = 0;
        if (!aux->Prox)
        {
          flag2 = 1;
          break;
        }
        aux = NULL;
      }
      else
      {
        flag = 1;
        aux = aux->Prox;
        pos++;
      }
    }

    if (flag == 1 || flag2 == 1)
    {
      //printf("O Elemento %c N?o existe na ED\n\n", elemento);
      if (Ultimo[q]->ele == 'f')
      {
        Ultimo[q]->ele = 'n';
      }
      // printf("pas");
      pos = 0;
    }
    return pos;
  }

  int SearchfF(char elemento, int pos_start)
  { // pesquisa por elemntos dentro da ED
    int pos = 1;
    int c = 1;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1, flag2 = 0;

    while (pos_start >= c && aux != NULL)
    {
      aux = aux->Prox;
      pos++;
      c++;
    }

    while (aux != NULL)
    {
      if (!aux->Prox)
      {
        flag = 1;
        break;
      }
      if (aux->ele == elemento && (aux->Ant->Ant->ele == ']' || aux->Ant->Ant->ele == '[') && aux->Prox->ele != 'B')
      {
        flag = 0;
        if (!aux->Prox)
        {
          break;
        }
        aux = NULL;
      }
      else
      {
        flag = 1;
        aux = aux->Prox;
        pos++;
      }
    }
    if (flag == 1)
    {
      pos = 0;
    }
    return pos;
  }

  
  int SearchR33(char elemento, int pos_start)
  { // pesquisa por elemntos dentro da ED
    int pos = 1;
    int c = 1;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;

    while (pos_start >= c && aux != NULL)
    {
      aux = aux->Prox;
      pos++;
      c++;
    }

    while (aux != NULL)
    {
      if (aux->ele == elemento && aux->Ant->Ant->ele == ']')
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        pos++;
      }
    }

    if (flag == 1)
    {
      //printf("O Elemento %c N?o existe na ED\n\n", elemento);
      pos = 0;
    }
    return pos;
  }

  int SearchAllf(char elemento, int pos_start)
  { // pesquisa por elemntos dentro da ED
    int pos = 1;
    int c = 1;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[q];
    int flag = 1;

    while (pos_start >= c && aux != NULL)
    {
      aux = aux->Prox;
      pos++;
      c++;
    }

    while (aux != NULL)
    {
      if (aux->ele == elemento && !(aux->Ant->ele == ']'))
      {
        //printf("O Elemento %c est? contido na ED \n", elemento);
        flag = 0;
        aux = NULL;
      }
      else
      {
        aux = aux->Prox;
        pos++;
      }
    }

    if (flag == 1 || pos >= size)
    {
      ////////   printf("O Elemento %c N?o existe na ED\n\n", elemento);
      pos = 0;
    }
    return pos;
  }
  
  int nf(int NENT, int NSAI, int NR, int tipo)
  {
    int c = 0;
    struct ED *aux = (struct ED *)calloc(NENT, sizeof(struct ED));
    aux = Primeiro[0];
    if (!aux->ele)
    {
      printf("A ED est? Vazia !");
    }
    else
      while (aux != NULL)
      {

        if (aux->ele == 'f' && tipo == 1)
          c++;
        else if (aux->ele == 'n')
          c++;
        aux = aux->Prox;
      }
    return (c - 1 - NSAI); // "1" por causa que os neuronios s? s?o gerados na 1 string e s?o confirmados nela msm
  }
  //////////////////////////////////////////////
  //////////////////////////////////////////////

  void SetQ(int x)
  {
    q = x;
  }

  void R1(int NENT)
  { // S -> *
    int x = Search('S');
    if (x != 0)
    {
      ReplaceOnPosition('.', 1);
    }
  }

  void R2(int NENT)
  { // * -> (f..f)*NENT
    int x = Search('.');

    if (x != 0)
    {
      ReplaceOnPosition('f', 1);
    }
  }

  void R31()
  { // f -> [f
    int x = Search('f');
    if (x != 0)
    {
      InsertInFronte('[');
    }
  }

  void R32(int NR)
  { // f -> f F f
    int c = NR;
    int x = 1;
    if (NSAI == 1)
    {
      x = Search('f');
      InsertMiddle('F', x + 1);
      InsertMiddle('f', x + 2);
    }
    else if (NSAI > 1)
    {
      while (x != 0 && c > 0)
      {
        x = SearchR322('f', x);
        if (x == 0)
          break;
        InsertMiddle('F', x + 1);
        InsertMiddle('f', x + 2);
        c--;
      }
    }
  } // fim R32

  void R322(int SIZE, int MAX, int NSAI)
  { // f -> f F f
    int Bigger;
    if (NSAI >= MAX)
    {
      Bigger = NSAI;
    }
    else
    {
      Bigger = MAX;
    }
    int c = Bigger;
    int x = 1;

    int i;
    for (i = 0; i < SIZE; i++)
    {
      x = 1;
      c = Bigger;
      while (x != 0 && c > 0)
      {
        x = SearchR322('f', x);
        //	  printf("%d while\n",x);
        //	  system("pause");
        if (x == 0)
          break;
        InsertMiddle('F', x + 1);
        InsertMiddle('f', x + 2);
        c--;
        ShowStack();
      }
      //	printf("i>>%d for\n",i);
      //	system("pause");
    }

  } // fim R32

  void R33(int NSAI)
  { // f -> f F
    int c = 1;
    int x = 1;
    while (c <= NSAI)
    {
      while (x != 0)
      {
        x = SearchR33('f', x);
        if (x == 0)
          break;
        InsertMiddle('F', x + 1);
        //cout << endl << endl << x  <<" < - "<< endl;
        c++;
      }
    }
  }

  //////////////////////////////////////////////////
  void R34(int NENT, int NSAI, int NR)
  {
    int nrr = NR;
    int x = 1;
    while (x != 0)
    {

      x = SearchR341('f', x); // coloca n B
      if (x == 0)
      {
        break;
      }
      else
      {
        ReplaceOnPosition('n', x);
      }
    }
    x = 1;
    while (x != 0 && nrr > 0)
    {

      x = SearchfF('f', x);
      if (x == 0)
      {

        break;
      }
      else
      {
      	
        ReplaceOnPosition('n', x);
        nrr--;
      }
    }
    if (NENT != 1)
    { // verifica se existe mais de 1 string v?lida
      while (q < NENT)
      {
        Ultimo[q]->Ant->ele = 'n';
        q++;
      }
    }
  }

  void R344(int NENT, int NSAI, int N, int SIZE, int MAX, int NR[])
  {

    int Bigger;
    int Vet[SIZE + 1];
    Vet[SIZE] = NSAI;

    int k;

    for (k = 0; k < SIZE; k++)
    {
      Vet[k] = NR[k];
      //		printf("k%d\n",Vet[k]);
    }
    //	printf("k%d\n",Vet[SIZE]);

    if (NSAI >= MAX)
    {
      Bigger = NSAI;
    }
    else
    {
      Bigger = MAX;
    }

    int nrr = N;
    int x = 1;
    int i = 0;
    int j = 0;

    //		printf("%d Size", SIZE);
    //for(i = 0 ; i < SIZE ; i++){
    //	printf("%d\n\n",NR[i]);
    //}

    while (x != 0)
    {
      x = SearchR341('f', x); // coloca n B
      if (x == 0)
      {
        break;
      }
      else
      {
        if (Vet[j] - i > 0)
        {
          ReplaceOnPosition('n', x);
        }
        else
          ReplaceOnPosition('f', x); //
          
      }
//      ReplaceOnPosition('0', x+1); //////////////////////////////////////////////////////////////////
		ReplaceOnPosition('n', SearchLastF('f'));
      ////////				   printf("i: %d  j: %d x: %d\n",i,j,x);
      j++;
      if (j >= SIZE + 1)
      {
        i++;
        j = 0;
      }
      if (i >= Bigger)
      {
        i = 0;
        return;
      }
      //        system("pause");
      //        ShowStack();    XX
    }
  }
  /////////////////////////////////////////////////////////////

  void R36(int NENT, int NSAI)
  {
    int x = 1;
    while (x != 0)
    {
      x = SearchR36('f', x);
      if (x != 0)
      {
        InsertMiddle('B', x + 1);
      }
    }
    x = 1;
    if (NENT > 1)
    {

      InsertLast('B');
    }
  }

  void R366(int NENT, int NSAI)
  {
    int x = 1;

    while (x != 0)
    {
      int x = SearchAllf('f', x);
      if (x != 0)
      {
        InsertMiddle('B', x + 1);
      }
      else
      {
        return;
      }
    }

    x = 1;
    if (NENT > 1)
    {
      InsertLast('B');
    }
  }

  void R4(int NSAI, int MAX)
  { // [ -> [Ff]
    int c = 1;
    int Bigger;
    if (NSAI >= MAX)
    {
      Bigger = NSAI;
    }
    else
    {
      Bigger = MAX;
    }
    while (c <= Bigger)
    {
      int x = Search('[');
      InsertMiddle('F', x + 1);
      InsertMiddle('f', x + 2);
      InsertMiddle(']', x + 3);
      c++;
    }
  } // fim R4

  ////////////////////////////////////////////
  ////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  //////////////////////// FIM REGRAS /////////////////////////
  /////////////////////////////////////////////////////////////

      printf("\n\n\n\nMAPEAMENTO GENOTIPO / FENOTIPO");
//  fprintf(pFile, "\n\n\n\nMAPEAMENTO GENOTIPO / FENOTIPO\n\n\n");

  int i;
  int MAX = NR[0];
  for (i = 0; i < SIZE; i++)
  {
    if (MAX < NR[i])
      MAX = NR[i];
  }

  QueueNull(); // Inicializa a Lista

  InsertLast('S'); // COLOCA 'S' NA PILHA

  // 1to1
  //////////////////////////////////////////////////////////////
  if (NENT == 1 && NSAI == 1 && tipo == 1)
  { // Type 1 to 1
     printf("\n\nTIPO 1x1\n\n");
//    fprintf(pFile, "\n\nTIPO 1x1\n\n");
    q = 0;
    if (NR[0] == 1)
    {
      // CASO ESPECIAL NENT = 1 | NSAI = 1 | NR = 1

      QueueNull();
      InsertLast('S'); // COLOCA 'S' NA PILHA

      R1(NENT);
      ShowStack();

      R2(NENT);
      ShowStack();

      R31();
      ShowStack();

      R4(1, MAX);
      ShowStack();

      R32(MAX);
      ShowStack();
    }
    else
    {
      // CASO NENT = 1 | NSAI = 1 | NR = n

      QueueNull();
      InsertLast('S'); // COLOCA 'S' NA PILHA

      R1(NENT);
      ShowStack();

      R2(NENT);
      ShowStack();

      R31();
      ShowStack();

      R4(1, MAX);
      ShowStack();

      R32(NSAI);
      ShowStack();

      R33(NSAI);
      ShowStack();

    } // fim else

  } // Fim Type 1 to 1
  //1toN
  //////////////////////////////////////////////////////////////////
  if (NENT == 1 && NSAI >= 1 && tipo == 2)
  {
    	printf("\n\nTIPO 1xN\n\n");
//    fprintf(pFile, "\n\nTIPO 1xN\n\n");

    q = 0;
    QueueNull();
    InsertLast('S'); // COLOCA 'S' NA PILHA

    R1(NENT);
    ShowStack();

    R2(NENT);
    ShowStack();

    R31();
    ShowStack();

    R4(NSAI, MAX);
    ShowStack();

    R32(NSAI);
    ShowStack();

    R36(NENT, NSAI);
    ShowStack();

    R34(NENT, NSAI, MAX);
    ShowStack();

  } // fim tipo 1.2

  if (NENT >= 1 && NSAI == 1 && tipo == 3)
  {
    	printf("\n\nTIPO Nx1\n\n");
//    fprintf(pFile, "\n\nTIPO Nx1\n\n");

    // para o Primeiro[0]

    q = 0;
    //  QueueNull();
    //  InsertLast('S'); // COLOCA 'S' NA PILHA
    while (q + 1 <= NENT)
    {
      //q++;
      QueueNull();
      InsertLast('S'); // COLOCA 'S' NA PILHA DAS OUTRAS STRINGS CASO EXISTAM
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R1(NENT); // S -> .  | R1 ? APLICADO NENT VEZES
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R2(NENT); // . -> f  | R2 ? APLICADO NENT VEZES
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R31(); // f -> [f  | R31 ? APLICADO EM CADA STRING
             //	printf("  ||  ");
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R4(NSAI, MAX); // [ -> [ F f  | NR > NSAI ENT?O ? APLICADO NR VEZES, SE N?O, ? APLICADO NSAI VEZES EM CADA STRING
      ShowStack();
      q++;
    }

    q = 0;
    R32(NSAI); // f -> f F f  | APLICADO NSAI VEZES SOMENTE PARA A PRIMEIRA STRING
    ShowStack();
    q++;
    while (q + 1 <= NENT)
    {
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R36(NENT, NSAI); // f -> f B  | APLICADO NR VEZES PARA O ULTIMO F DE CADA BLOCO E APLICADO PARA O PRIMEIRO f DE
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      if (q == 0)
      {
        R34(NENT, NSAI, MAX); // f -> n  | APLICADO PARA OS PRIMEIROS f DE CADA BLOCO E APLICADO NSAI VEZES PARA O PRIMEIRO f
        q = 0;
      }
      ShowStack();
      q++;
    }
  } //  fim tipo 1.3

  if (NENT >= 1 && NSAI >= 1 && tipo == 4)
  {
    	printf("\n\nTIPO MxN\n\n");
//    fprintf(pFile, "\n\nTIPO MxN\n\n");
    // para o Primeiro[0]
    q = 0;
    QueueNull();
    InsertLast('S'); // COLOCA 'S' NA PILHA

    while (q + 1 <= NENT)
    {
      q++;
      QueueNull();
      InsertLast('S'); // COLOCA 'S' NA PILHA DAS OUTRAS STRINGS CASO EXISTAM
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R1(NENT); // S -> .  | R1 ? APLICADO NENT VEZES
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R2(NENT); // . -> f  | R2 ? APLICADO NENT VEZES
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R31(); // f -> [f  | R31 ? APLICADO EM CADA STRING
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R4(NSAI, MAX); // [ -> [ F f  | NR > NSAI ENT?O ? APLICADO NR VEZES, SE N?O, ? APLICADO NSAI VEZES EM CADA STRING
      ShowStack();
      q++;
    }

    q = 0;
    while (q + 1 <= NENT)
    {
      R322(SIZE, MAX, NSAI); // f -> f F f  | APLICADO NSAI VEZES SOMENTE PARA A PRIMEIRA STRING /////// alterado
      //R32(NSAI); <-Antigo
      ShowStack();
      q++;
    }

    //			q = 0;
    //			while(q+1 <= NENT){
    //				R32(NSAI);  // f -> f F f  | APLICADO NSAI VEZES SOMENTE PARA A PRIMEIRA STRING
    //				ShowStack();
    //		    	q++;
    //			}

    q = 0;
    while (q + 1 <= NENT)
    {
      R366(NENT, NSAI); // f -> f B  | APLICADO NR VEZES PARA O ULTIMO F DE CADA BLOCO E APLICADO PARA O PRIMEIRO f DE
                        // R36; <-Antigo
      ShowStack();
      q++;
    }

    //			q = 0;
    //			while(q+1 <= NENT){
    //				R36(NENT, NSAI); // f -> f B  | APLICADO NR VEZES PARA O ULTIMO F DE CADA BLOCO E APLICADO PARA O PRIMEIRO f DE
    //		     	ShowStack();
    //			 	q++;
    //			}

    q = 0;
    while (q + 1 <= NENT)
    {
      if (q == 0)
      {
        R344(NENT, NSAI, 5, SIZE, MAX, NR);
        //R34(NENT,NSAI,NR); // f -> n  | APLICADO PARA OS PRIMEIROS f DE CADA BLOCO E APLICADO NSAI VEZES PARA O PRIMEIRO f
      }
      Ultimo[q]->ele = 'n';
      ShowStack();
      q++;
    }
  
    
    

  } 
  
  
  
  	printf("\n| Numero de Neuronios Camada Intermediaria 1 -> %d \n", NR[0]);
//  fprintf(pFile, "| Numero de Neuronios Camada Intermediaria 1 -> %d \n", NR[0]);
  if (n == 4)
//    fprintf(pFile, "| Numero de Neuronios Camada Intermediaria 2 -> %d \n", NR[1]);
    printf("| Numero de Neuronios Camada Intermediaria 2 -> %d \n", NR[1]);
  //    fprintf(pFile,"| Numero de Neuronios por Ramificacao -> %d \n", NR);
//  fprintf(pFile, "| NENT -> %d \n", NENT);
//  fprintf(pFile, "| NSAI -> %d \n\n", NSAI);
  printf("| NENT -> %d \n", NENT);
  printf("| NSAI -> %d \n\n", NSAI);

  return 0;

} // Fim Mapeamento

