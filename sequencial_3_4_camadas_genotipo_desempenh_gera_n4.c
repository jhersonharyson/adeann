#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//
#include<string.h>
#include <windows.h>
//
//50
//Parametros de Inicializacao
int geracao=1;//numero de rega??es
int individuos=60;//epresenta o numero de individuos do genotipo nesse caso 100
int gene=516;
//gene/pop=4
const int NINTER = 10000;
int tipo=4;
//int subtipo=1;
int NINT=4; //E3uivale a 2 neuronios na intermediaria 1 e bias assumira NENT*aleatorio do outro modulo variaveis.

int n=4;

int NINT1 = 7;
int NINT2 = 7;
//***PC Neuronios de entrada saida e numero de padroes, taxa de aprendizado, numero de epcocas e erro aceitavel
const int NENT = 5;// 3; //Equivale a 2 neuronios na entradas e 1 bias
const int NSAI = 1;//Equilave a 1 saida
int NPAD = 900;//numero de padroes
double TAPR3 = 0.01;
double TAPR4 = 0.01; // 7
unsigned long MAXITER = 500000;//1800000//00;//numero de epocas;
const double ACEITAVEL = 0.001;
int  contind=0;
double contregrasval=0;
int indic=0;
int sup=0;
int neuint=0;
//***PC tax de translocacao
double tx_trans=0.001;
//
double err=0.0;
double emq=0.0;
double fitness=0.0;
int aux=0;
double soma_nint_total=0.0;

//MODELO DE DESENVOLVIMENTO SISTEMA-L
void mapeamento_genotipo_fenotipo(int NENT,int NSAI,int aleatorio,int tipo,FILE *pFile);
void genotipo_estatico(int individuos,int gene);
int** genotipo_dinamico(int individuos,int gene);
char** genotipo_dinamico_string(int individuos,int gene_dec);
void imprime_genbin(int** gen,int individuos,int gene,FILE *pfile);
void legenes_genbin(int** gen,int** gen_bin_dec,int individuos,int gene);
int converte_genbindec(int gene1,int gene2,int gene3,int gene4,int gene5,int gene6,int j,int compactador);
void imprime_genbindec(int** gen_dec,int individuos,int gene_dec);
void imprime_regras_genotipo_string(char** gen_string,int individuos,int gene_dec,int aleatorio);
void libera_gen_bin(int** gen,int individuos);
void libera_gen_bin_dec(int** gen_bin_dec,int individuos);
void legenes_genbindec_string(int** gen_bin_dec,char** gen_string,int individuos,int gene_dec);
void avalia_regras_gen_string(char** gen_string,int individuos,int gene_dec);
void imprime_genstring(char** gen_string,int individuos,int gene_dec,FILE *Pfile);
void libera_gen_string(char** gen_string,int individuos);

//REDE NEURAL
void conjunto_treinamento1(double x[][NENT],double y[][NSAI]);
void inicializa_pesos(double w1[][NINTER],double w2[][NSAI]);
void intermediaria(double x[][NENT],double w1[][NINTER],double h[],int m);
void treina_rede(int individuos,FILE *pFile,int NINT);
void saida(double h[],double w2[][NSAI],double o[]);
double erro_saida(double o[],double y[][NSAI],int m);
void erro2(double o[],double y[][NSAI],int m,double delta2[]);
void erro1(double h[],double delta2[],double w2[][NSAI],double delta1[]);
void ajusta2(double w2[][NSAI],double delta2[],double h[]);
void ajusta1(double w1[][NINTER],double delta1[],double x[][NENT],int m);
void verifica1(double x[][NENT],double w1[][NINTER],double w2[][NSAI],double y[][NSAI],FILE *pFile);


//REDE NEURAL VERIFICA 2
void conjunto_treinamento1(double x[][NENT],double y[][NSAI]);

void treina_rede_(int individuos,FILE *pFile,int NINT1, int NINT2);

void inicializa_pesos_(double w1[][NINTER],double w2[][NSAI], double w3[NINT1][NINT2]);
void intermediaria_(double x[][NENT],double w1[][NINTER],double h1[],int m);
void intermediaria2_(double w3[NINT1][NINT2], double h1[],double h2[] ,int m);
void saida_(double h2[],double w2[][NSAI],double o[]);
double erro_saida_(double o[],double y[][NSAI],int m);

void ajusta_(double w3[][NINTER],double delta1[],double h1[],int m);
void ajusta1_(double w1[][NINTER],double delta[],double x[][NENT],int m);
void ajusta2_(double w2[][NSAI],double delta2[],double h2[]);
void erro2_(double o[],double y[][NSAI],int m,double delta2[]);
void erro1_(double h2[],double delta1[],double w2[NINT1][NINT2],double delta[]);
void erroN_(double h1[],double delta2[],double w3[][NSAI],double delta1[]);
void ajusta3_(double w2[][NSAI],double delta2[],double h1[]);
void verifica1_(double x[][NENT],double w1[][NINTER],double w2[][NSAI], double w3[NINT1][NINT2],double y[][NSAI],FILE *pFile);







//ALGORITMO GEN?TICO
double fit[100][10];
double hist_fit[100][2];
void imprime_fitness(double fit[][10],FILE* pFile,int contador);
void zera_fitness(double fit[][10]);
void imprime_hist_fitness(double hist_fit[][2],FILE* pFile,int contador);
void imprime_cabec(FILE* pFile);
void ordena(double fit[][10],FILE* pFile);
void cruzamento(int i,int j,int pos,int** gen,char** gen_string,int gene_dec);
void mutacao(int** gen);
void selecao(int** gen,char** gen_string,int gene_dec);
void transloca(char** gen_string,int gene_dec);

// METODOS MAPEAMENTO GENOTIPO FENOTIPO
int Mapeamento(int NENT, int NSAI, int NR[],int SIZE, char TIPO[], FILE *pFile);
int NR_RAND();
int NR_RAND_N(int n);
//char *TIPO = "1.3";
//


unsigned char tab_converte[]={'f','F','n','.','n','.','f','F','F','f','B','f','[','n','[','.',
			'f',']','n','*','.','F','f','F',']','.','[','f','f','*','B',']',
			'.',']','n','F','f','B','f','B','F','[','B','n','*','f','.',']',
			']','[','n','F','n','B','[','.','f',']','B','F','B','f','*','['};
//Alterar o tamanho 100 sempre que mudar numero de individuos
//


int main()
{
    FILE *pFile;
//    	pFile = fopen("relat","w");
pFile = fopen("relat_autoMaasd2asdasd.txt","w");
    for(n = 4; n <= 4; n++){
   // 	for( TAPR3 = 0.1f ; TAPR3 <= 0.9f ;TAPR3 += 0.15f){
    //		for(TAPR4 = 0.1f ; TAPR4 <= 0.9f ;TAPR4 += 0.15f){
//			printf("TARP1 %f\nTARP2 %f\n\n", TAPR3,TAPR4);
	if(n==3)
		pFile = fopen("relat_N3asdas.txt","w");
    else{
//		float valor1 = TAPR3;
//		float valor2 = TAPR4;
//    	char T1[16];
//    	char T2[16];
//    	sprintf(T1, "%.2f", valor1);
//    	sprintf(T2, "%.2f", valor2);
//    	char * arq =  strcat(strcat(strcat(strcat("relat_N4_",T1) ,"_"), T2),".txt");
    	
//		pFile = fopen("relat_auto.txt","w");
		fprintf(pFile,"\nTARP1 %.2f\nTARP2 %.2f\n\n\n");
	}
	// gene : represneta o numero de genes por inidividuo, nesse caso 516;
    int gene_dec=(gene/6); //43+1=44 ultimo elemento armazena status da string valida % = valida;
//////    printf("Gerando Genotipo Aleatoriamente!\n");
    int** gen_bin=genotipo_dinamico(individuos,gene);
    int** gen_bin_dec=genotipo_dinamico(individuos,gene_dec);
    char** gen_string=genotipo_dinamico_string(individuos,gene_dec);
    zera_fitness(fit);
    int contador=1;//contador do numero de geracoes
    int contador1;
    //int pontos_corte=(gene_dec*0.9);

    while(contador<=geracao)
    {

       imprime_genbin(gen_bin,individuos,gene,pFile);
       legenes_genbin(gen_bin,gen_bin_dec,individuos,gene);
       imprime_genbindec(gen_bin_dec,individuos,gene_dec);
       legenes_genbindec_string(gen_bin_dec,gen_string,individuos,gene_dec);
       contador1=1;//contador do sumero de cruzamentos por gera??o
       zera_fitness(fit);
       avalia_regras_gen_string(gen_string,individuos,gene_dec);//Avalia Regras Validas
       imprime_genstring(gen_string,individuos,gene_dec,pFile); //Imprime Individuo[i] e DNA
       //Ordenar Individuos
	   ordena(fit,pFile);
       //Imprimir Individuos
       imprime_fitness(fit,pFile,contador);
       contador++;
       indic=0;
       while(contador1<=((gene_dec*0.8)))
       {
       // printf("contador1=%d",contador1);
       // getch();
        //realiza (individuos/2) cruzamentos
        selecao(gen_bin,gen_string,gene_dec);
        contador1++;
        
       	//printf("%d", contador1);

       }
       contador1=1;
       
       while(contador1<=(gene_dec*0.8))
       {
         //516-86-14
         mutacao(gen_bin);
         contador1++;

       }

      // transloca(gen_string,gene_dec);

       if (geracao<20)
         MAXITER=MAXITER+25000;
       else if((geracao>=20)&&(geracao<40))
               MAXITER=MAXITER+80000;
             else
               MAXITER=MAXITER+150000;
    }
    contador=1;


    imprime_cabec(pFile);
    while(contador<=(geracao))
    {
      imprime_hist_fitness(hist_fit,pFile,contador);
      contador++;
    }
    //End Do
    //LIBERA MEM?RIA ALOCADA DINAMICAMENTE
    libera_gen_bin(gen_bin,individuos);
    libera_gen_bin_dec(gen_bin_dec,individuos);
    libera_gen_string(gen_string,individuos);
//////    printf("\n\n<<Simulacao Concluida - Relatorio Gerado!!>>");
    fprintf(pFile,"\nsimulacao concluida");
    
  
	}
//}
//}
  fclose (pFile);
	return 0;
}

void genotipo_estatico(int individuos,int gene)
{
  //n?o esta sendo usado
  int i,j,aleatorio;
  time_t t;

  int gen[individuos][gene];
  srand((unsigned) time(&t));

  for(i=0;i<gene;i++)
    for(j=0;j<individuos;j++)
    {
      aleatorio=(rand()%2);
      gen[i][j]=aleatorio;
    }

  for(i=0;i<gene;i++)
    for(j=0;j<individuos;j++)
    {
//////      printf("%d",gen[i][j]);
    }
}

int** genotipo_dinamico(int individuos,int gene)
{
  int i,j,aleatorio;
  int** gen;

  time_t t;
  srand((unsigned) time(&t));

  gen=(int**)malloc(individuos*sizeof(int*));
  if(gen==NULL)
   {
//////    printf("Memoria Insuficiente.\n");
     return NULL;
   }
   for(i=0;i<individuos;i++)
   {
    gen[i]=(int*)malloc(gene*sizeof(int));
    if(gen[i]==NULL)
    {
//////    printf("Memoria Insuficiente.\n");
     return NULL;
    }
   }

  for(i=0;i<individuos;i++)
    for(j=0;j<gene;j++)
    {
      aleatorio=(rand()%2);
      gen[i][j]=aleatorio;
    }
  return (gen);

}

int** gen_bin_dec_genotipo_dinamico(individuos,gene_dec)
{
  int i;
  int** gen_dec;

  gen_dec=(int**)malloc(individuos*sizeof(int*));
  if(gen_dec==NULL)
   {
//////    printf("Memoria Insuficiente.\n");
     return NULL;
   }
   for(i=0;i<individuos;i++)
   {
    gen_dec[i]=(int*)malloc(gene_dec*sizeof(int));
    if(gen_dec[i]==NULL)
    {
//////    printf("Memoria Insuficiente.\n");
     return NULL;
    }
   }
   return(gen_dec);
}

char** genotipo_dinamico_string(individuos,gene_dec)
{
  int i;
  char** gen_string;

  gen_string=(char**)malloc(individuos*sizeof(int*));
  if(gen_string==NULL)
   {
//////    printf("Memoria Insuficiente.\n");
     return NULL;
   }
   for(i=0;i<individuos;i++)
   {
    gen_string[i]=(char*)malloc((gene_dec+1)*sizeof(int));
    if(gen_string[i]==NULL)
    {
//////     printf("Memoria Insuficiente.\n");
     return NULL;
    }
   }
   return(gen_string);

}
void imprime_genbin(int** gen,int individuos,int gene,FILE *pFile)
{
  int i,j;
 // printf("\n");
  for(i=0;i<individuos;i++)
    for(j=0;j<gene;j++)
    {
      //fprintf(pFile,"encerrado");
//////      printf("%d",gen[i][j]);
    }

}

void imprime_genbindec(int** gen_bin_dec,int individuos,int gene_dec)
{
  int i,j=0;
  //printf("\n");
  for(i=0;i<individuos;i++)
  {
//////    printf("-Individuo[%d]-\n",i+1);
    while(j<gene_dec)
    {
//////     printf("%d\t",gen_bin_dec[i][j]);
     j++;
    }
  //  printf("\n\n\n");
    j=0;
  }

}

void imprime_genstring(char** gen_string,int individuos,int gene_dec,FILE *pFile)
{

  int i,aleatorio,j=0;
  time_t t;
  srand((unsigned) time(&t));
  //zera_fitness(fit);

  for(i=0;i<individuos;i++)
  {
//////    printf("\n-Individuo[%d]-\n",i+1);
    fprintf(pFile,"\n-Individuo[%d]-\n",i+1);

    while(j<gene_dec)
    {
//     printf("%c",gen_string[i][j]);

     fprintf(pFile,"%c",gen_string[i][j]);
     j++;
      
    }

    if(gen_string[i][gene_dec+1]=='V')
    {
       contind++;
//       printf("\t STRING VALIDA");
       fprintf(pFile,"\t STRING VALIDA");
       //******PC - N?mero de neuronios por ramifica??o
       //aleatorio=(rand()%5);
       //while((aleatorio<NSAI))
       //   aleatorio++;
       mapeamento_genotipo_fenotipo(NENT,NSAI,aleatorio,tipo,pFile);
      

    }
    else
    {
    	
          contind++;
          fit[indic][0]=(indic+1);
          fit[indic][1]=0.0;
          fit[indic][2]=fit[indic][0];
          aux=indic+1;
          fit[indic][3]=aux;
          if(aux>1)
            fit[indic][3]=aux+fit[indic-1][3];
          else
            fit[indic][3]=aux;
          fit[indic][4]=0.0;
          fit[indic][5]=0.0;
          fit[indic][6]=0.0;
          indic++;
//////          printf("\t STRING INVALIDA");
          fprintf(pFile,"\t STRING INVALIDA");
          
     }
     j=0;
     
  }
 
}


void legenes_genbin(int** gen,int** gen_bin_dec,int individuos,int gene)
{
  int gene1,gene2,gene3,gene4,gene5,gene6,i;
  int start=0;
  int j=0;
  int compactador =0;
  int stop=(gene/6);

  while(j<individuos)
  {
   for(i=0;i<stop;i++)
   {
    gene1=gen[j][start];
    gene2=gen[j][start+1];
    gene3=gen[j][start+2];
    gene4=gen[j][start+3];
    gene5=gen[j][start+4];
    gene6=gen[j][start+5];
    int decimal=converte_genbindec(gene1,gene2,gene3,gene4,gene5,gene6,j,compactador);
    gen_bin_dec[j][compactador]=decimal;
    compactador++;
    start=start+6;
   }
   j++;
   start=0;
   compactador=0;
  }
}

int converte_genbindec(int gene1,int gene2,int gene3,int gene4,int gene5,int gene6,int j,int compactador)
{
  int decimal=gene6*32+gene5*16+gene4*8+gene3*4+gene2*2+gene1;
  return decimal;
}

void legenes_genbindec_string(int** gen_bin_dec,char** gen_string,int individuos,int gene_dec)
{
  int i,aux,aux1;
  int j=0;
  while(j<individuos)
  {
   for(i=0;i<gene_dec;i++)
   {
     aux=gen_bin_dec[j][i];
     aux1=tab_converte[aux];
     gen_string[j][i]=aux1;
   }
   j++;
  }
}

void avalia_regras_gen_string(char** gen_string,int individuos,int gene_dec)
{
  int i=0;
  int j=0;
  unsigned char string_val[]={'.','f','[','F','f','n','B',']'};
  char carac;
  int pos_string;
  int cont_elem;

  while(j<individuos)
  {
   for(i=0;i<(gene_dec-6);i++)
   {
    //Le do inicio para o final com um passo de um bit.
    if((gen_string[j][i]=='.')&&(gen_string[j][i+1]=='f')&&(gen_string[j][i+2]=='[')&&(gen_string[j][i+3]=='F')&&(gen_string[j][i+4]=='f')&&(gen_string[j][i+5]=='n')&&(gen_string[j][i+6]=='B')&&(gen_string[j][i+7]==']'))
    {

      gen_string[j][gene_dec+1]='V';

    }
   }
   j++;
  }

  j=0;

  while(j<individuos)
  {
   //Le do final para o inicio com um passo de um bit.
   for(i=gene_dec;i>=5;i--)
   {
    if((gen_string[j][i]=='.')&&(gen_string[j][i+1]=='f')&&(gen_string[j][i+2]=='[')&&(gen_string[j][i+3]=='F')&&(gen_string[j][i+4]=='f')&&(gen_string[j][i+5]=='n')&&(gen_string[j][i+6]=='B')&&(gen_string[j][i+7]==']'))
    {
      gen_string[j][gene_dec+1]='V';

    }
   }
   j++;
  }

  j=0;

  //mudei aqui baixo todo la?o for


  for(j=0;j<individuos;j++)
  {
   //Le de forma sequencial bit a bit do inicio para o fim.
    cont_elem=0;
    pos_string=0;

    for(i=0;i<(gene_dec);i++)
    {
     carac=gen_string[j][i];

     if(string_val[pos_string]==carac)
     {
      cont_elem++;
      pos_string++;
     }

     if(cont_elem==8)
     {
       gen_string[j][gene_dec+1]='V';
       //printf("j=%d",j);
       //printf("cont_elem=%d",cont_elem);
       //getch();
     }
    }
  }
}


void mapeamento_genotipo_fenotipo(int NENT,int NSAI,int aleatorio,int tipo,FILE *pFile)
{
 //contregrasval++;
  int NR = 7;
  int NR1 =7;  
  
  if(NR == NR1)
  	NR1 = NR_RAND(3);
  
  NINT = NR;
  NINT1 = NR;
  NINT2 = NR1;
  
  int NINT_N3[1];
  NINT_N3[0] = NINT1;
  int SIZE_N3 = sizeof(NINT_N3)/sizeof(NINT_N3[0]);
  
  int NINT_N4[2];
  NINT_N4[0] = NINT1;
  NINT_N4[1] = NINT2;
  int SIZE_N4 = sizeof(NINT_N4)/sizeof(NINT_N4[0]);
  //printf("\n\n%d %d %d",NENT,NSAI,NINT);
  //printf("%d %d %d %d", NENT, NSAI, NR, tipo);
  //Mapeamento( NENT, NSAI, NR, tipo, pFile);// NENT = 1 NSAI = 1 NINT = 4 NR = NR = NR_RAND();
    n == 3 ?
	Mapeamento(NENT, NSAI,  NINT_N3, SIZE_N3, "1.4",pFile):
	Mapeamento(NENT, NSAI,  NINT_N4, SIZE_N4, "1.4",pFile);  // (ENTRADA, SAIDA, NR, TIPO)
  //system("pause");
  //Sleep(2000);
  //system("pause");
  
  if(n==3)
  treina_rede(contind,pFile,NINT);
  if(n==4)
  treina_rede_(contind,pFile,NINT1,NINT2);
}

void libera_gen_bin(int** gen,int individuos)
{
 int i;
 for(i=0;i<individuos;i++)
 	free(gen[i]);
  free(gen);
}

void libera_gen_bin_dec(int** gen_bin_dec,int individuos)
{
 int i;
 for(i=0;i<individuos;i++)
 	free(gen_bin_dec[i]);
  free(gen_bin_dec);
}

void libera_gen_string(char** gen_string,int individuos)
{
 int i;
 for(i=0;i<individuos;i++)
 	free(gen_string[i]);
  free(gen_string);
}

 void treina_rede(int individuos,FILE *pFile,int NINT)
 {
//////   printf("\nNENT=(%d Entradas + 1 Bias)",NENT);
   fprintf(pFile,"\nNENT=(%d Entradas + 1 Bias)",NENT-1);
//////   printf("\nNINT=(%d Int + 1 Bias)",NINT-1);
   fprintf(pFile,"\nNINT=(%d Int + 1 Bias)",NINT-1);
//////   printf("\nNSAI=%d",NSAI);
   fprintf(pFile,"\nNSAI=%d",NSAI);
//////   printf("\n\nTreinamento do Individuo=%d",individuos);
   fprintf(pFile,"\n\nTreinamento do Individuo=%d",individuos);

   double x[NPAD][NENT],h[NINTER],o[NSAI],y[NPAD][NSAI],delta1[NINTER],delta2[NSAI],w1[NENT][NINTER],w2[NINTER][NSAI],erromax;
   int m;
   int l;
   /////////////////////////////////////system("cls");
   fprintf(pFile,"\n<<Aprendizagem da Rede>>");
//////   printf("\n<<Aprendizagem da Rede\n>>");
   conjunto_treinamento1(x,y);
   inicializa_pesos(w1,w2);
   erromax = ACEITAVEL*2.0;
   m = 0;
   l = 0;
   while(erromax > ACEITAVEL && l < MAXITER)
   {
      if(m == NPAD) m = 0;
      intermediaria(x,w1,h,m);
      saida(h,w2,o);
      err=erro_saida(o,y,m);
      if(err > erromax) erromax = err;
      erro2(o,y,m,delta2);
      erro1(h,delta2,w2,delta1);
      ajusta2(w2,delta2,h);
      ajusta1(w1,delta1,x,m);
      l = l + 1;
//////      printf("\nPadrao>>%d",m);
      m = m + 1;
//////      printf("\t\t Epoca:%d \t\t Erro:%2.9f",l,err);
    }
   // getch();
    fprintf(pFile,"\n<<Verificacao da Aprendizagem>>\n");
//////    printf("\n<<Verificacao da Aprendizagem>>\n");
    verifica1(x,w1,w2,y,pFile);
   // getch();
 }

void conjunto_treinamento1(double x[][NENT],double y[][NSAI])
{ float a = 1, b = 0;
{

//TREINO

 x[0][0] = 0.001;
 x[0][1] = 3.6216;
 x[0][2] = 8.6661;
 x[0][3] = -2.8073;
 x[0][4] = -0.44699;
 y[0][0] = a;

 x[1][0] = 0.001;
 x[1][1] = 4.5459;
 x[1][2] = 8.1674;
 x[1][3] = -2.4586;
 x[1][4] = -1.4621;
 y[1][0] = a;

 x[2][0] = 0.001;
 x[2][1] = 3.866;
 x[2][2] = -2.6383;
 x[2][3] = 1.9242;
 x[2][4] = 0.10645;
 y[2][0] = a;

 x[3][0] = 0.001;
 x[3][1] = 3.4566;
 x[3][2] = 9.5228;
 x[3][3] = -4.0112;
 x[3][4] = -3.5944;
 y[3][0] = a;

 x[4][0] = 0.001;
 x[4][1] = 0.32924;
 x[4][2] = -4.4552;
 x[4][3] = 4.5718;
 x[4][4] = -0.9888;
 y[4][0] = a;

 x[5][0] = 0.001;
 x[5][1] = 4.3684;
 x[5][2] = 9.6718;
 x[5][3] = -3.9606;
 x[5][4] = -3.1625;
 y[5][0] = a;

 x[6][0] = 0.001;
 x[6][1] = 3.5912;
 x[6][2] = 3.0129;
 x[6][3] = 0.72888;
 x[6][4] = 0.56421;
 y[6][0] = a;

 x[7][0] = 0.001;
 x[7][1] = 2.0922;
 x[7][2] = -6.81;
 x[7][3] = 8.4636;
 x[7][4] = -0.60216;
 y[7][0] = a;

 x[8][0] = 0.001;
 x[8][1] = 3.2032;
 x[8][2] = 5.7588;
 x[8][3] = -0.75345;
 x[8][4] = -0.61251;
 y[8][0] = a;

 x[9][0] = 0.001;
 x[9][1] = 1.5356;
 x[9][2] = 9.1772;
 x[9][3] = -2.2718;
 x[9][4] = -0.73535;
 y[9][0] = a;

 x[10][0] = 0.001;
 x[10][1] = 1.2247;
 x[10][2] = 8.7779;
 x[10][3] = -2.2135;
 x[10][4] = -0.80647;
 y[10][0] = a;

 x[11][0] = 0.001;
 x[11][1] = 3.9899;
 x[11][2] = -2.7066;
 x[11][3] = 2.3946;
 x[11][4] = 0.86291;
 y[11][0] = a;

 x[12][0] = 0.001;
 x[12][1] = 1.8993;
 x[12][2] = 7.6625;
 x[12][3] = 0.15394;
 x[12][4] = -3.1108;
 y[12][0] = a;

 x[13][0] = 0.001;
 x[13][1] = -1.5768;
 x[13][2] = 10.843;
 x[13][3] = 2.5462;
 x[13][4] = -2.9362;
 y[13][0] = a;

 x[14][0] = 0.001;
 x[14][1] = 3.404;
 x[14][2] = 8.7261;
 x[14][3] = -2.9915;
 x[14][4] = -0.57242;
 y[14][0] = a;

 x[15][0] = 0.001;
 x[15][1] = 4.6765;
 x[15][2] = -3.3895;
 x[15][3] = 3.4896;
 x[15][4] = 1.4771;
 y[15][0] = a;

 x[16][0] = 0.001;
 x[16][1] = 2.6719;
 x[16][2] = 3.0646;
 x[16][3] = 0.37158;
 x[16][4] = 0.58619;
 y[16][0] = a;

 x[17][0] = 0.001;
 x[17][1] = 0.80355;
 x[17][2] = 2.8473;
 x[17][3] = 4.3439;
 x[17][4] = 0.6017;
 y[17][0] = a;

 x[18][0] = 0.001;
 x[18][1] = 1.4479;
 x[18][2] = -4.8794;
 x[18][3] = 8.3428;
 x[18][4] = -2.1086;
 y[18][0] = a;

 x[19][0] = 0.001;
 x[19][1] = 5.2423;
 x[19][2] = 11.0272;
 x[19][3] = -4.353;
 x[19][4] = -4.1013;
 y[19][0] = a;

 x[20][0] = 0.001;
 x[20][1] = 5.7867;
 x[20][2] = 7.8902;
 x[20][3] = -2.6196;
 x[20][4] = -0.48708;
 y[20][0] = a;

 x[21][0] = 0.001;
 x[21][1] = 0.3292;
 x[21][2] = -4.4552;
 x[21][3] = 4.5718;
 x[21][4] = -0.9888;
 y[21][0] = a;

 x[22][0] = 0.001;
 x[22][1] = 3.9362;
 x[22][2] = 10.1622;
 x[22][3] = -3.8235;
 x[22][4] = -4.0172;
 y[22][0] = a;

 x[23][0] = 0.001;
 x[23][1] = 0.93584;
 x[23][2] = 8.8855;
 x[23][3] = -1.6831;
 x[23][4] = -1.6599;
 y[23][0] = a;

 x[24][0] = 0.001;
 x[24][1] = 4.4338;
 x[24][2] = 9.887;
 x[24][3] = -4.6795;
 x[24][4] = -3.7483;
 y[24][0] = a;

 x[25][0] = 0.001;
 x[25][1] = 0.7057;
 x[25][2] = -5.4981;
 x[25][3] = 8.3368;
 x[25][4] = -2.8715;
 y[25][0] = a;

 x[26][0] = 0.001;
 x[26][1] = 1.1432;
 x[26][2] = -3.7413;
 x[26][3] = 5.5777;
 x[26][4] = -0.63578;
 y[26][0] = a;

 x[27][0] = 0.001;
 x[27][1] = -0.38214;
 x[27][2] = 8.3909;
 x[27][3] = 2.1624;
 x[27][4] = -3.7405;
 y[27][0] = a;

 x[28][0] = 0.001;
 x[28][1] = 6.5633;
 x[28][2] = 9.8187;
 x[28][3] = -4.4113;
 x[28][4] = -3.2258;
 y[28][0] = a;

 x[29][0] = 0.001;
 x[29][1] = 4.8906;
 x[29][2] = -3.3584;
 x[29][3] = 3.4202;
 x[29][4] = 1.0905;
 y[29][0] = a;

 x[30][0] = 0.001;
 x[30][1] = -0.24811;
 x[30][2] = -0.17797;
 x[30][3] = 4.9068;
 x[30][4] = 0.15429;
 y[30][0] = a;

 x[31][0] = 0.001;
 x[31][1] = 1.4884;
 x[31][2] = 3.6274;
 x[31][3] = 3.308;
 x[31][4] = 0.48921;
 y[31][0] = a;

 x[32][0] = 0.001;
 x[32][1] = 4.2969;
 x[32][2] = 7.617;
 x[32][3] = -2.3874;
 x[32][4] = -0.96164;
 y[32][0] = a;

 x[33][0] = 0.001;
 x[33][1] = -0.96511;
 x[33][2] = 9.4111;
 x[33][3] = 1.7305;
 x[33][4] = -4.8629;
 y[33][0] = a;

 x[34][0] = 0.001;
 x[34][1] = -1.6162;
 x[34][2] = 0.80908;
 x[34][3] = 8.1628;
 x[34][4] = 0.60817;
 y[34][0] = a;

 x[35][0] = 0.001;
 x[35][1] = 2.4391;
 x[35][2] = 6.4417;
 x[35][3] = -0.80743;
 x[35][4] = -0.69139;
 y[35][0] = a;

 x[36][0] = 0.001;
 x[36][1] = 2.6881;
 x[36][2] = 6.0195;
 x[36][3] = -0.46641;
 x[36][4] = -0.69268;
 y[36][0] = a;

 x[37][0] = 0.001;
 x[37][1] = 3.6289;
 x[37][2] = 0.81322;
 x[37][3] = 1.6277;
 x[37][4] = 0.77627;
 y[37][0] = a;

 x[38][0] = 0.001;
 x[38][1] = 4.5679;
 x[38][2] = 3.1929;
 x[38][3] = -2.1055;
 x[38][4] = 0.29653;
 y[38][0] = a;

 x[39][0] = 0.001;
 x[39][1] = 3.4805;
 x[39][2] = 9.7008;
 x[39][3] = -3.7541;
 x[39][4] = -3.4379;
 y[39][0] = a;

 x[40][0] = 0.001;
 x[40][1] = 4.1711;
 x[40][2] = 8.722;
 x[40][3] = -3.0224;
 x[40][4] = -0.59699;
 y[40][0] = a;

 x[41][0] = 0.001;
 x[41][1] = -0.2062;
 x[41][2] = 9.2207;
 x[41][3] = -3.7044;
 x[41][4] = -6.8103;
 y[41][0] = a;

 x[42][0] = 0.001;
 x[42][1] = -0.006892;
 x[42][2] = 9.2931;
 x[42][3] = -0.41243;
 x[42][4] = -1.9638;
 y[42][0] = a;

 x[43][0] = 0.001;
 x[43][1] = 0.96441;
 x[43][2] = 5.8395;
 x[43][3] = 2.3235;
 x[43][4] = 0.066365;
 y[43][0] = a;

 x[44][0] = 0.001;
 x[44][1] = 2.8561;
 x[44][2] = 6.9176;
 x[44][3] = -0.79372;
 x[44][4] = 0.48403;
 y[44][0] = a;

 x[45][0] = 0.001;
 x[45][1] = -0.7869;
 x[45][2] = 9.5663;
 x[45][3] = -3.7867;
 x[45][4] = -7.5034;
 y[45][0] = a;

 x[46][0] = 0.001;
 x[46][1] = 2.0843;
 x[46][2] = 6.6258;
 x[46][3] = 0.48382;
 x[46][4] = -2.2134;
 y[46][0] = a;

 x[47][0] = 0.001;
 x[47][1] = -0.7869;
 x[47][2] = 9.5663;
 x[47][3] = -3.7867;
 x[47][4] = -7.5034;
 y[47][0] = a;

 x[48][0] = 0.001;
 x[48][1] = 3.9102;
 x[48][2] = 6.065;
 x[48][3] = -2.4534;
 x[48][4] = -0.68234;
 y[48][0] = a;

 x[49][0] = 0.001;
 x[49][1] = 1.6349;
 x[49][2] = 3.286;
 x[49][3] = 2.8753;
 x[49][4] = 0.087054;
 y[49][0] = a;

 x[50][0] = 0.001;
 x[50][1] = 4.3239;
 x[50][2] = -4.8835;
 x[50][3] = 3.4356;
 x[50][4] = -0.5776;
 y[50][0] = a;

 x[51][0] = 0.001;
 x[51][1] = 5.262;
 x[51][2] = 3.9834;
 x[51][3] = -1.5572;
 x[51][4] = 1.0103;
 y[51][0] = a;

 x[52][0] = 0.001;
 x[52][1] = 3.1452;
 x[52][2] = 5.825;
 x[52][3] = -0.51439;
 x[52][4] = -1.4944;
 y[52][0] = a;

 x[53][0] = 0.001;
 x[53][1] = 2.549;
 x[53][2] = 6.1499;
 x[53][3] = -1.1605;
 x[53][4] = -1.2371;
 y[53][0] = a;

 x[54][0] = 0.001;
 x[54][1] = 4.9264;
 x[54][2] = 5.496;
 x[54][3] = -2.4774;
 x[54][4] = -0.50648;
 y[54][0] = a;

 x[55][0] = 0.001;
 x[55][1] = 4.8265;
 x[55][2] = 0.80287;
 x[55][3] = 1.6371;
 x[55][4] = 1.1875;
 y[55][0] = a;

 x[56][0] = 0.001;
 x[56][1] = 2.5635;
 x[56][2] = 6.7769;
 x[56][3] = -0.61979;
 x[56][4] = 0.38576;
 y[56][0] = a;

 x[57][0] = 0.001;
 x[57][1] = 5.807;
 x[57][2] = 5.0097;
 x[57][3] = -2.2384;
 x[57][4] = 0.43878;
 y[57][0] = a;

 x[58][0] = 0.001;
 x[58][1] = 3.1377;
 x[58][2] = -4.1096;
 x[58][3] = 4.5701;
 x[58][4] = 0.98963;
 y[58][0] = a;

 x[59][0] = 0.001;
 x[59][1] = -0.78289;
 x[59][2] = 11.3603;
 x[59][3] = -0.37644;
 x[59][4] = -7.0495;
 y[59][0] = a;

 x[60][0] = 0.001;
 x[60][1] = 2.888;
 x[60][2] = 0.44696;
 x[60][3] = 4.5907;
 x[60][4] = -0.24398;
 y[60][0] = a;

 x[61][0] = 0.001;
 x[61][1] = 0.49665;
 x[61][2] = 5.527;
 x[61][3] = 1.7785;
 x[61][4] = -0.47156;
 y[61][0] = a;

 x[62][0] = 0.001;
 x[62][1] = 4.2586;
 x[62][2] = 11.2962;
 x[62][3] = -4.0943;
 x[62][4] = -4.3457;
 y[62][0] = a;

 x[63][0] = 0.001;
 x[63][1] = 1.7939;
 x[63][2] = -1.1174;
 x[63][3] = 1.5454;
 x[63][4] = -0.26079;
 y[63][0] = a;

 x[64][0] = 0.001;
 x[64][1] = 5.4021;
 x[64][2] = 3.1039;
 x[64][3] = -1.1536;
 x[64][4] = 1.5651;
 y[64][0] = a;

 x[65][0] = 0.001;
 x[65][1] = 2.5367;
 x[65][2] = 2.599;
 x[65][3] = 2.0938;
 x[65][4] = 0.20085;
 y[65][0] = a;

 x[66][0] = 0.001;
 x[66][1] = 4.6054;
 x[66][2] = -4.0765;
 x[66][3] = 2.7587;
 x[66][4] = 0.31981;
 y[66][0] = a;

 x[67][0] = 0.001;
 x[67][1] = 2.4235;
 x[67][2] = 9.5332;
 x[67][3] = -3.0789;
 x[67][4] = -2.7746;
 y[67][0] = a;

 x[68][0] = 0.001;
 x[68][1] = 1.0009;
 x[68][2] = 7.7846;
 x[68][3] = -0.28219;
 x[68][4] = -2.6608;
 y[68][0] = a;

 x[69][0] = 0.001;
 x[69][1] = 0.12326;
 x[69][2] = 8.9848;
 x[69][3] = -0.9351;
 x[69][4] = -2.4332;
 y[69][0] = a;

 x[70][0] = 0.001;
 x[70][1] = 3.9529;
 x[70][2] = -2.3548;
 x[70][3] = 2.3792;
 x[70][4] = 0.48274;
 y[70][0] = a;

 x[71][0] = 0.001;
 x[71][1] = 4.1373;
 x[71][2] = 0.49248;
 x[71][3] = 1.093;
 x[71][4] = 1.8276;
 y[71][0] = a;

 x[72][0] = 0.001;
 x[72][1] = 4.7181;
 x[72][2] = 10.0153;
 x[72][3] = -3.9486;
 x[72][4] = -3.8582;
 y[72][0] = a;

 x[73][0] = 0.001;
 x[73][1] = 4.1654;
 x[73][2] = -3.4495;
 x[73][3] = 3.643;
 x[73][4] = 1.0879;
 y[73][0] = a;

 x[74][0] = 0.001;
 x[74][1] = 4.4069;
 x[74][2] = 10.9072;
 x[74][3] = -4.5775;
 x[74][4] = -4.4271;
 y[74][0] = a;

 x[75][0] = 0.001;
 x[75][1] = 2.3066;
 x[75][2] = 3.5364;
 x[75][3] = 0.57551;
 x[75][4] = 0.41938;
 y[75][0] = a;

 x[76][0] = 0.001;
 x[76][1] = 3.7935;
 x[76][2] = 7.9853;
 x[76][3] = -2.5477;
 x[76][4] = -1.872;
 y[76][0] = a;

 x[77][0] = 0.001;
 x[77][1] = 0.049175;
 x[77][2] = 6.1437;
 x[77][3] = 1.7828;
 x[77][4] = -0.72113;
 y[77][0] = a;

 x[78][0] = 0.001;
 x[78][1] = 0.24835;
 x[78][2] = 7.6439;
 x[78][3] = 0.9885;
 x[78][4] = -0.87371;
 y[78][0] = a;

 x[79][0] = 0.001;
 x[79][1] = 1.1317;
 x[79][2] = 3.9647;
 x[79][3] = 3.3979;
 x[79][4] = 0.84351;
 y[79][0] = a;

 x[80][0] = 0.001;
 x[80][1] = 2.8033;
 x[80][2] = 9.0862;
 x[80][3] = -3.3668;
 x[80][4] = -1.0224;
 y[80][0] = a;

 x[81][0] = 0.001;
 x[81][1] = 4.4682;
 x[81][2] = 2.2907;
 x[81][3] = 0.95766;
 x[81][4] = 0.83058;
 y[81][0] = a;

 x[82][0] = 0.001;
 x[82][1] = 5.0185;
 x[82][2] = 8.5978;
 x[82][3] = -2.9375;
 x[82][4] = -1.281;
 y[82][0] = a;

 x[83][0] = 0.001;
 x[83][1] = 1.8664;
 x[83][2] = 7.7763;
 x[83][3] = -0.23849;
 x[83][4] = -2.9634;
 y[83][0] = a;

 x[84][0] = 0.001;
 x[84][1] = 3.245;
 x[84][2] = 6.63;
 x[84][3] = -0.63435;
 x[84][4] = 0.86937;
 y[84][0] = a;

 x[85][0] = 0.001;
 x[85][1] = 4.0296;
 x[85][2] = 2.6756;
 x[85][3] = 0.80685;
 x[85][4] = 0.71679;
 y[85][0] = a;

 x[86][0] = 0.001;
 x[86][1] = -1.1313;
 x[86][2] = 1.9037;
 x[86][3] = 7.5339;
 x[86][4] = 1.022;
 y[86][0] = a;

 x[87][0] = 0.001;
 x[87][1] = 0.87603;
 x[87][2] = 6.8141;
 x[87][3] = 0.84198;
 x[87][4] = -0.17156;
 y[87][0] = a;

 x[88][0] = 0.001;
 x[88][1] = 4.1197;
 x[88][2] = -2.7956;
 x[88][3] = 2.0707;
 x[88][4] = 0.67412;
 y[88][0] = a;

 x[89][0] = 0.001;
 x[89][1] = 3.8027;
 x[89][2] = 0.81529;
 x[89][3] = 2.1041;
 x[89][4] = 1.0245;
 y[89][0] = a;

 x[90][0] = 0.001;
 x[90][1] = 1.4806;
 x[90][2] = 7.6377;
 x[90][3] = -2.7876;
 x[90][4] = -1.0341;
 y[90][0] = a;

 x[91][0] = 0.001;
 x[91][1] = 4.0632;
 x[91][2] = 3.584;
 x[91][3] = 0.72545;
 x[91][4] = 0.39481;
 y[91][0] = a;

 x[92][0] = 0.001;
 x[92][1] = 4.3064;
 x[92][2] = 8.2068;
 x[92][3] = -2.7824;
 x[92][4] = -1.4336;
 y[92][0] = a;

 x[93][0] = 0.001;
 x[93][1] = 2.4486;
 x[93][2] = -6.3175;
 x[93][3] = 7.9632;
 x[93][4] = 0.20602;
 y[93][0] = a;

 x[94][0] = 0.001;
 x[94][1] = 3.2718;
 x[94][2] = 1.7837;
 x[94][3] = 2.1161;
 x[94][4] = 0.61334;
 y[94][0] = a;

 x[95][0] = 0.001;
 x[95][1] = -0.64472;
 x[95][2] = -4.6062;
 x[95][3] = 8.347;
 x[95][4] = -2.7099;
 y[95][0] = a;

 x[96][0] = 0.001;
 x[96][1] = 2.9543;
 x[96][2] = 1.076;
 x[96][3] = 0.64577;
 x[96][4] = 0.89394;
 y[96][0] = a;

 x[97][0] = 0.001;
 x[97][1] = 2.1616;
 x[97][2] = -6.8804;
 x[97][3] = 8.1517;
 x[97][4] = -0.081048;
 y[97][0] = a;

 x[98][0] = 0.001;
 x[98][1] = 3.82;
 x[98][2] = 10.9279;
 x[98][3] = -4.0112;
 x[98][4] = -5.0284;
 y[98][0] = a;

 x[99][0] = 0.001;
 x[99][1] = -2.7419;
 x[99][2] = 11.4038;
 x[99][3] = 2.5394;
 x[99][4] = -5.5793;
 y[99][0] = a;

 x[100][0] = 0.001;
 x[100][1] = 3.3669;
 x[100][2] = -5.1856;
 x[100][3] = 3.6935;
 x[100][4] = -1.1427;
 y[100][0] = a;

 x[101][0] = 0.001;
 x[101][1] = 4.5597;
 x[101][2] = -2.4211;
 x[101][3] = 2.6413;
 x[101][4] = 1.6168;
 y[101][0] = a;

 x[102][0] = 0.001;
 x[102][1] = 5.1129;
 x[102][2] = -0.49871;
 x[102][3] = 0.62863;
 x[102][4] = 1.1189;
 y[102][0] = a;

 x[103][0] = 0.001;
 x[103][1] = 3.3397;
 x[103][2] = -4.6145;
 x[103][3] = 3.9823;
 x[103][4] = -0.23751;
 y[103][0] = a;

 x[104][0] = 0.001;
 x[104][1] = 4.2027;
 x[104][2] = 0.22761;
 x[104][3] = 0.96108;
 x[104][4] = 0.97282;
 y[104][0] = a;

 x[105][0] = 0.001;
 x[105][1] = 3.5438;
 x[105][2] = 1.2395;
 x[105][3] = 1.997;
 x[105][4] = 2.1547;
 y[105][0] = a;

 x[106][0] = 0.001;
 x[106][1] = 2.3136;
 x[106][2] = 10.6651;
 x[106][3] = -3.5288;
 x[106][4] = -4.7672;
 y[106][0] = a;

 x[107][0] = 0.001;
 x[107][1] = -1.8584;
 x[107][2] = 7.886;
 x[107][3] = -1.6643;
 x[107][4] = -1.8384;
 y[107][0] = a;

 x[108][0] = 0.001;
 x[108][1] = 3.106;
 x[108][2] = 9.5414;
 x[108][3] = -4.2536;
 x[108][4] = -4.003;
 y[108][0] = a;

 x[109][0] = 0.001;
 x[109][1] = 2.9163;
 x[109][2] = 10.8306;
 x[109][3] = -3.3437;
 x[109][4] = -4.122;
 y[109][0] = a;

 x[110][0] = 0.001;
 x[110][1] = 3.9922;
 x[110][2] = -4.4676;
 x[110][3] = 3.7304;
 x[110][4] = -0.1095;
 y[110][0] = a;

 x[111][0] = 0.001;
 x[111][1] = 1.518;
 x[111][2] = 5.6946;
 x[111][3] = 0.094818;
 x[111][4] = -0.026738;
 y[111][0] = a;

 x[112][0] = 0.001;
 x[112][1] = 3.2351;
 x[112][2] = 9.647;
 x[112][3] = -3.2074;
 x[112][4] = -2.5948;
 y[112][0] = a;

 x[113][0] = 0.001;
 x[113][1] = 4.2188;
 x[113][2] = 6.8162;
 x[113][3] = -1.2804;
 x[113][4] = 0.76076;
 y[113][0] = a;

 x[114][0] = 0.001;
 x[114][1] = 1.7819;
 x[114][2] = 6.9176;
 x[114][3] = -1.2744;
 x[114][4] = -1.5759;
 y[114][0] = a;

 x[115][0] = 0.001;
 x[115][1] = 2.5331;
 x[115][2] = 2.9135;
 x[115][3] = -0.822;
 x[115][4] = -0.12243;
 y[115][0] = a;

 x[116][0] = 0.001;
 x[116][1] = 3.8969;
 x[116][2] = 7.4163;
 x[116][3] = -1.8245;
 x[116][4] = 0.14007;
 y[116][0] = a;

 x[117][0] = 0.001;
 x[117][1] = 2.108;
 x[117][2] = 6.7955;
 x[117][3] = -0.1708;
 x[117][4] = 0.4905;
 y[117][0] = a;

 x[118][0] = 0.001;
 x[118][1] = 2.8969;
 x[118][2] = 0.70768;
 x[118][3] = 2.29;
 x[118][4] = 1.8663;
 y[118][0] = a;

 x[119][0] = 0.001;
 x[119][1] = 0.9297;
 x[119][2] = -3.7971;
 x[119][3] = 4.6429;
 x[119][4] = -0.2957;
 y[119][0] = a;

 x[120][0] = 0.001;
 x[120][1] = 3.4642;
 x[120][2] = 10.6878;
 x[120][3] = -3.4071;
 x[120][4] = -4.109;
 y[120][0] = a;

 x[121][0] = 0.001;
 x[121][1] = 4.0713;
 x[121][2] = 10.4023;
 x[121][3] = -4.1722;
 x[121][4] = -4.7582;
 y[121][0] = a;

 x[122][0] = 0.001;
 x[122][1] = -1.4572;
 x[122][2] = 9.1214;
 x[122][3] = 1.7425;
 x[122][4] = -5.1241;
 y[122][0] = a;

 x[123][0] = 0.001;
 x[123][1] = -1.5075;
 x[123][2] = 1.9224;
 x[123][3] = 7.1466;
 x[123][4] = 0.89136;
 y[123][0] = a;

 x[124][0] = 0.001;
 x[124][1] = -0.91718;
 x[124][2] = 9.9884;
 x[124][3] = 1.1804;
 x[124][4] = -5.2263;
 y[124][0] = a;

 x[125][0] = 0.001;
 x[125][1] = 2.994;
 x[125][2] = 7.2011;
 x[125][3] = -1.2153;
 x[125][4] = 0.3211;
 y[125][0] = a;

 x[126][0] = 0.001;
 x[126][1] = -2.343;
 x[126][2] = 12.9516;
 x[126][3] = 3.3285;
 x[126][4] = -5.9426;
 y[126][0] = a;

 x[127][0] = 0.001;
 x[127][1] = 3.7818;
 x[127][2] = -2.8846;
 x[127][3] = 2.2558;
 x[127][4] = -0.15734;
 y[127][0] = a;

 x[128][0] = 0.001;
 x[128][1] = 4.6689;
 x[128][2] = 1.3098;
 x[128][3] = 0.055404;
 x[128][4] = 1.909;
 y[128][0] = a;

 x[129][0] = 0.001;
 x[129][1] = 3.4663;
 x[129][2] = 1.1112;
 x[129][3] = 1.7425;
 x[129][4] = 1.3388;
 y[129][0] = a;

 x[130][0] = 0.001;
 x[130][1] = 3.2697;
 x[130][2] = -4.3414;
 x[130][3] = 3.6884;
 x[130][4] = -0.29829;
 y[130][0] = a;

 x[131][0] = 0.001;
 x[131][1] = 5.1302;
 x[131][2] = 8.6703;
 x[131][3] = -2.8913;
 x[131][4] = -1.5086;
 y[131][0] = a;

 x[132][0] = 0.001;
 x[132][1] = 2.0139;
 x[132][2] = 6.1416;
 x[132][3] = 0.37929;
 x[132][4] = 0.56938;
 y[132][0] = a;

 x[133][0] = 0.001;
 x[133][1] = 0.4339;
 x[133][2] = 5.5395;
 x[133][3] = 2.033;
 x[133][4] = -0.40432;
 y[133][0] = a;

 x[134][0] = 0.001;
 x[134][1] = -1.0401;
 x[134][2] = 9.3987;
 x[134][3] = 0.85998;
 x[134][4] = -5.3336;
 y[134][0] = a;

 x[135][0] = 0.001;
 x[135][1] = 4.1605;
 x[135][2] = 11.2196;
 x[135][3] = -3.6136;
 x[135][4] = -4.0819;
 y[135][0] = a;

 x[136][0] = 0.001;
 x[136][1] = 5.438;
 x[136][2] = 9.4669;
 x[136][3] = -4.9417;
 x[136][4] = -3.9202;
 y[136][0] = a;

 x[137][0] = 0.001;
 x[137][1] = 5.032;
 x[137][2] = 8.2026;
 x[137][3] = -2.6256;
 x[137][4] = -1.0341;
 y[137][0] = a;

 x[138][0] = 0.001;
 x[138][1] = 5.2418;
 x[138][2] = 10.5388;
 x[138][3] = -4.1174;
 x[138][4] = -4.2797;
 y[138][0] = a;

 x[139][0] = 0.001;
 x[139][1] = -0.2062;
 x[139][2] = 9.2207;
 x[139][3] = -3.7044;
 x[139][4] = -6.8103;
 y[139][0] = a;

 x[140][0] = 0.001;
 x[140][1] = 2.0911;
 x[140][2] = 0.94358;
 x[140][3] = 4.5512;
 x[140][4] = 1.234;
 y[140][0] = a;

 x[141][0] = 0.001;
 x[141][1] = 1.7317;
 x[141][2] = -0.34765;
 x[141][3] = 4.1905;
 x[141][4] = -0.99138;
 y[141][0] = a;

 x[142][0] = 0.001;
 x[142][1] = 4.1736;
 x[142][2] = 3.3336;
 x[142][3] = -1.4244;
 x[142][4] = 0.60429;
 y[142][0] = a;

 x[143][0] = 0.001;
 x[143][1] = 3.9232;
 x[143][2] = -3.2467;
 x[143][3] = 3.4579;
 x[143][4] = 0.83705;
 y[143][0] = a;

 x[144][0] = 0.001;
 x[144][1] = 3.8481;
 x[144][2] = 10.1539;
 x[144][3] = -3.8561;
 x[144][4] = -4.2228;
 y[144][0] = a;

 x[145][0] = 0.001;
 x[145][1] = 0.5195;
 x[145][2] = -3.2633;
 x[145][3] = 3.0895;
 x[145][4] = -0.9849;
 y[145][0] = a;

 x[146][0] = 0.001;
 x[146][1] = 3.8584;
 x[146][2] = 0.78425;
 x[146][3] = 1.1033;
 x[146][4] = 1.7008;
 y[146][0] = a;

 x[147][0] = 0.001;
 x[147][1] = 1.7496;
 x[147][2] = -0.1759;
 x[147][3] = 5.1827;
 x[147][4] = 1.2922;
 y[147][0] = a;

 x[148][0] = 0.001;
 x[148][1] = 3.6277;
 x[148][2] = 0.9829;
 x[148][3] = 0.68861;
 x[148][4] = 0.63403;
 y[148][0] = a;

 x[149][0] = 0.001;
 x[149][1] = 2.7391;
 x[149][2] = 7.4018;
 x[149][3] = 0.071684;
 x[149][4] = -2.5302;
 y[149][0] = a;

 x[150][0] = 0.001;
 x[150][1] = 4.5447;
 x[150][2] = 8.2274;
 x[150][3] = -2.4166;
 x[150][4] = -1.5875;
 y[150][0] = a;

 x[151][0] = 0.001;
 x[151][1] = -1.7599;
 x[151][2] = 11.9211;
 x[151][3] = 2.6756;
 x[151][4] = -3.3241;
 y[151][0] = a;

 x[152][0] = 0.001;
 x[152][1] = 5.0691;
 x[152][2] = 0.21313;
 x[152][3] = 0.20278;
 x[152][4] = 1.2095;
 y[152][0] = a;

 x[153][0] = 0.001;
 x[153][1] = 3.4591;
 x[153][2] = 11.112;
 x[153][3] = -4.2039;
 x[153][4] = -5.0931;
 y[153][0] = a;

 x[154][0] = 0.001;
 x[154][1] = 1.9358;
 x[154][2] = 8.1654;
 x[154][3] = -0.023425;
 x[154][4] = -2.2586;
 y[154][0] = a;

 x[155][0] = 0.001;
 x[155][1] = 2.486;
 x[155][2] = -0.99533;
 x[155][3] = 5.3404;
 x[155][4] = -0.15475;
 y[155][0] = a;

 x[156][0] = 0.001;
 x[156][1] = 2.4226;
 x[156][2] = -4.5752;
 x[156][3] = 5.947;
 x[156][4] = 0.21507;
 y[156][0] = a;

 x[157][0] = 0.001;
 x[157][1] = 3.9479;
 x[157][2] = -3.7723;
 x[157][3] = 2.883;
 x[157][4] = 0.019813;
 y[157][0] = a;

 x[158][0] = 0.001;
 x[158][1] = 2.2634;
 x[158][2] = -4.4862;
 x[158][3] = 3.6558;
 x[158][4] = -0.61251;
 y[158][0] = a;

 x[159][0] = 0.001;
 x[159][1] = 1.3566;
 x[159][2] = 4.2358;
 x[159][3] = 2.1341;
 x[159][4] = 0.3211;
 y[159][0] = a;

 x[160][0] = 0.001;
 x[160][1] = 5.0452;
 x[160][2] = 3.8964;
 x[160][3] = -1.4304;
 x[160][4] = 0.86291;
 y[160][0] = a;

 x[161][0] = 0.001;
 x[161][1] = 3.5499;
 x[161][2] = 8.6165;
 x[161][3] = -3.2794;
 x[161][4] = -1.2009;
 y[161][0] = a;

 x[162][0] = 0.001;
 x[162][1] = 0.17346;
 x[162][2] = 7.8695;
 x[162][3] = 0.26876;
 x[162][4] = -3.7883;
 y[162][0] = a;

 x[163][0] = 0.001;
 x[163][1] = 2.4008;
 x[163][2] = 9.3593;
 x[163][3] = -3.3565;
 x[163][4] = -3.3526;
 y[163][0] = a;

 x[164][0] = 0.001;
 x[164][1] = 4.8851;
 x[164][2] = 1.5995;
 x[164][3] = -0.000291;
 x[164][4] = 1.6401;
 y[164][0] = a;

 x[165][0] = 0.001;
 x[165][1] = 4.1927;
 x[165][2] = -3.2674;
 x[165][3] = 2.5839;
 x[165][4] = 0.21766;
 y[165][0] = a;

 x[166][0] = 0.001;
 x[166][1] = 1.1166;
 x[166][2] = 8.6496;
 x[166][3] = -0.96252;
 x[166][4] = -1.8112;
 y[166][0] = a;

 x[167][0] = 0.001;
 x[167][1] = 1.0235;
 x[167][2] = 6.901;
 x[167][3] = -2.0062;
 x[167][4] = -2.7125;
 y[167][0] = a;

 x[168][0] = 0.001;
 x[168][1] = -1.803;
 x[168][2] = 11.8818;
 x[168][3] = 2.0458;
 x[168][4] = -5.2728;
 y[168][0] = a;

 x[169][0] = 0.001;
 x[169][1] = 0.11739;
 x[169][2] = 6.2761;
 x[169][3] = -1.5495;
 x[169][4] = -2.4746;
 y[169][0] = a;

 x[170][0] = 0.001;
 x[170][1] = 0.5706;
 x[170][2] = -0.0248;
 x[170][3] = 1.2421;
 x[170][4] = -0.5621;
 y[170][0] = a;

 x[171][0] = 0.001;
 x[171][1] = 4.0552;
 x[171][2] = -2.4583;
 x[171][3] = 2.2806;
 x[171][4] = 1.0323;
 y[171][0] = a;

 x[172][0] = 0.001;
 x[172][1] = -1.6952;
 x[172][2] = 1.0657;
 x[172][3] = 8.8294;
 x[172][4] = 0.94955;
 y[172][0] = a;

 x[173][0] = 0.001;
 x[173][1] = -1.1193;
 x[173][2] = 10.7271;
 x[173][3] = 2.0938;
 x[173][4] = -5.6504;
 y[173][0] = a;

 x[174][0] = 0.001;
 x[174][1] = 1.8799;
 x[174][2] = 2.4707;
 x[174][3] = 2.4931;
 x[174][4] = 0.37671;
 y[174][0] = a;

 x[175][0] = 0.001;
 x[175][1] = 3.583;
 x[175][2] = -3.7971;
 x[175][3] = 3.4391;
 x[175][4] = -0.12501;
 y[175][0] = a;

 x[176][0] = 0.001;
 x[176][1] = 0.19081;
 x[176][2] = 9.1297;
 x[176][3] = -3.725;
 x[176][4] = -5.8224;
 y[176][0] = a;

 x[177][0] = 0.001;
 x[177][1] = 3.6582;
 x[177][2] = 5.6864;
 x[177][3] = -1.7157;
 x[177][4] = -0.23751;
 y[177][0] = a;

 x[178][0] = 0.001;
 x[178][1] = -0.13144;
 x[178][2] = -1.7775;
 x[178][3] = 8.3316;
 x[178][4] = 0.35214;
 y[178][0] = a;

 x[179][0] = 0.001;
 x[179][1] = 2.3925;
 x[179][2] = 9.798;
 x[179][3] = -3.0361;
 x[179][4] = -2.8224;
 y[179][0] = a;

 x[180][0] = 0.001;
 x[180][1] = 1.6426;
 x[180][2] = 3.0149;
 x[180][3] = 0.22849;
 x[180][4] = -0.147;
 y[180][0] = a;

 x[181][0] = 0.001;
 x[181][1] = -0.11783;
 x[181][2] = -1.5789;
 x[181][3] = 8.03;
 x[181][4] = -0.028031;
 y[181][0] = a;

 x[182][0] = 0.001;
 x[182][1] = -0.69572;
 x[182][2] = 8.6165;
 x[182][3] = 1.8419;
 x[182][4] = -4.3289;
 y[182][0] = a;

 x[183][0] = 0.001;
 x[183][1] = 2.9421;
 x[183][2] = 7.4101;
 x[183][3] = -0.97709;
 x[183][4] = -0.88406;
 y[183][0] = a;

 x[184][0] = 0.001;
 x[184][1] = -1.7559;
 x[184][2] = 11.9459;
 x[184][3] = 3.0946;
 x[184][4] = -4.8978;
 y[184][0] = a;

 x[185][0] = 0.001;
 x[185][1] = -1.2537;
 x[185][2] = 10.8803;
 x[185][3] = 1.931;
 x[185][4] = -4.3237;
 y[185][0] = a;

 x[186][0] = 0.001;
 x[186][1] = 3.2585;
 x[186][2] = -4.4614;
 x[186][3] = 3.8024;
 x[186][4] = -0.15087;
 y[186][0] = a;

 x[187][0] = 0.001;
 x[187][1] = 1.8314;
 x[187][2] = 6.3672;
 x[187][3] = -0.036278;
 x[187][4] = 0.049554;
 y[187][0] = a;

 x[188][0] = 0.001;
 x[188][1] = 4.5645;
 x[188][2] = -3.6275;
 x[188][3] = 2.8684;
 x[188][4] = 0.27714;
 y[188][0] = a;

 x[189][0] = 0.001;
 x[189][1] = 2.7365;
 x[189][2] = -5.0325;
 x[189][3] = 6.6608;
 x[189][4] = -0.57889;
 y[189][0] = a;

 x[190][0] = 0.001;
 x[190][1] = 0.9297;
 x[190][2] = -3.7971;
 x[190][3] = 4.6429;
 x[190][4] = -0.2957;
 y[190][0] = a;

 x[191][0] = 0.001;
 x[191][1] = 3.9663;
 x[191][2] = 10.1684;
 x[191][3] = -4.1131;
 x[191][4] = -4.6056;
 y[191][0] = a;

 x[192][0] = 0.001;
 x[192][1] = 1.4578;
 x[192][2] = -0.08485;
 x[192][3] = 4.1785;
 x[192][4] = 0.59136;
 y[192][0] = a;

 x[193][0] = 0.001;
 x[193][1] = 4.8272;
 x[193][2] = 3.0687;
 x[193][3] = 0.68604;
 x[193][4] = 0.80731;
 y[193][0] = a;

 x[194][0] = 0.001;
 x[194][1] = -2.341;
 x[194][2] = 12.3784;
 x[194][3] = 0.70403;
 x[194][4] = -7.5836;
 y[194][0] = a;

 x[195][0] = 0.001;
 x[195][1] = -1.8584;
 x[195][2] = 7.886;
 x[195][3] = -1.6643;
 x[195][4] = -1.8384;
 y[195][0] = a;

 x[196][0] = 0.001;
 x[196][1] = 4.1454;
 x[196][2] = 7.257;
 x[196][3] = -1.9153;
 x[196][4] = -0.86078;
 y[196][0] = a;

 x[197][0] = 0.001;
 x[197][1] = 1.9157;
 x[197][2] = 6.0816;
 x[197][3] = 0.23705;
 x[197][4] = -2.0116;
 y[197][0] = a;

 x[198][0] = 0.001;
 x[198][1] = 4.0215;
 x[198][2] = -2.1914;
 x[198][3] = 2.4648;
 x[198][4] = 1.1409;
 y[198][0] = a;

 x[199][0] = 0.001;
 x[199][1] = 5.8862;
 x[199][2] = 5.8747;
 x[199][3] = -2.8167;
 x[199][4] = -0.30087;
 y[199][0] = a;

 x[200][0] = 0.001;
 x[200][1] = -2.0897;
 x[200][2] = 10.8265;
 x[200][3] = 2.3603;
 x[200][4] = -3.4198;
 y[200][0] = a;

 x[201][0] = 0.001;
 x[201][1] = 4.0026;
 x[201][2] = -3.5943;
 x[201][3] = 3.5573;
 x[201][4] = 0.26809;
 y[201][0] = a;

 x[202][0] = 0.001;
 x[202][1] = -0.78689;
 x[202][2] = 9.5663;
 x[202][3] = -3.7867;
 x[202][4] = -7.5034;
 y[202][0] = a;

 x[203][0] = 0.001;
 x[203][1] = 4.1757;
 x[203][2] = 10.2615;
 x[203][3] = -3.8552;
 x[203][4] = -4.3056;
 y[203][0] = a;

 x[204][0] = 0.001;
 x[204][1] = 0.83292;
 x[204][2] = 7.5404;
 x[204][3] = 0.65005;
 x[204][4] = -0.92544;
 y[204][0] = a;

 x[205][0] = 0.001;
 x[205][1] = 4.8077;
 x[205][2] = 2.2327;
 x[205][3] = -0.26334;
 x[205][4] = 1.5534;
 y[205][0] = a;

 x[206][0] = 0.001;
 x[206][1] = 5.3063;
 x[206][2] = 5.2684;
 x[206][3] = -2.8904;
 x[206][4] = -0.52716;
 y[206][0] = a;

 x[207][0] = 0.001;
 x[207][1] = 2.5605;
 x[207][2] = 9.2683;
 x[207][3] = -3.5913;
 x[207][4] = -1.356;
 y[207][0] = a;

 x[208][0] = 0.001;
 x[208][1] = 2.1059;
 x[208][2] = 7.6046;
 x[208][3] = -0.47755;
 x[208][4] = -1.8461;
 y[208][0] = a;

 x[209][0] = 0.001;
 x[209][1] = 2.1721;
 x[209][2] = -0.73874;
 x[209][3] = 5.4672;
 x[209][4] = -0.72371;
 y[209][0] = a;

 x[210][0] = 0.001;
 x[210][1] = 4.2899;
 x[210][2] = 9.1814;
 x[210][3] = -4.6067;
 x[210][4] = -4.3263;
 y[210][0] = a;

 x[211][0] = 0.001;
 x[211][1] = 3.5156;
 x[211][2] = 10.1891;
 x[211][3] = -4.2759;
 x[211][4] = -4.978;
 y[211][0] = a;

 x[212][0] = 0.001;
 x[212][1] = 2.614;
 x[212][2] = 8.0081;
 x[212][3] = -3.7258;
 x[212][4] = -1.3069;
 y[212][0] = a;

 x[213][0] = 0.001;
 x[213][1] = 0.68087;
 x[213][2] = 2.3259;
 x[213][3] = 4.9085;
 x[213][4] = 0.54998;
 y[213][0] = a;

 x[214][0] = 0.001;
 x[214][1] = 4.1962;
 x[214][2] = 0.74493;
 x[214][3] = 0.83256;
 x[214][4] = 0.753;
 y[214][0] = a;

 x[215][0] = 0.001;
 x[215][1] = 6.0919;
 x[215][2] = 2.9673;
 x[215][3] = -1.3267;
 x[215][4] = 1.4551;
 y[215][0] = a;

 x[216][0] = 0.001;
 x[216][1] = 1.3234;
 x[216][2] = 3.2964;
 x[216][3] = 0.2362;
 x[216][4] = -0.11984;
 y[216][0] = a;

 x[217][0] = 0.001;
 x[217][1] = 1.3264;
 x[217][2] = 1.0326;
 x[217][3] = 5.6566;
 x[217][4] = -0.41337;
 y[217][0] = a;

 x[218][0] = 0.001;
 x[218][1] = -0.16735;
 x[218][2] = 7.6274;
 x[218][3] = 1.2061;
 x[218][4] = -3.6241;
 y[218][0] = a;

 x[219][0] = 0.001;
 x[219][1] = -1.3;
 x[219][2] = 10.2678;
 x[219][3] = -2.953;
 x[219][4] = -5.8638;
 y[219][0] = a;

 x[220][0] = 0.001;
 x[220][1] = -2.2261;
 x[220][2] = 12.5398;
 x[220][3] = 2.9438;
 x[220][4] = -3.5258;
 y[220][0] = a;

 x[221][0] = 0.001;
 x[221][1] = 2.4196;
 x[221][2] = 6.4665;
 x[221][3] = -0.75688;
 x[221][4] = 0.228;
 y[221][0] = a;

 x[222][0] = 0.001;
 x[222][1] = 1.0987;
 x[222][2] = 0.6394;
 x[222][3] = 5.989;
 x[222][4] = -0.58277;
 y[222][0] = a;

 x[223][0] = 0.001;
 x[223][1] = 4.6464;
 x[223][2] = 10.5326;
 x[223][3] = -4.5852;
 x[223][4] = -4.206;
 y[223][0] = a;

 x[224][0] = 0.001;
 x[224][1] = -0.36038;
 x[224][2] = 4.1158;
 x[224][3] = 3.1143;
 x[224][4] = -0.37199;
 y[224][0] = a;

 x[225][0] = 0.001;
 x[225][1] = 1.3562;
 x[225][2] = 3.2136;
 x[225][3] = 4.3465;
 x[225][4] = 0.78662;
 y[225][0] = a;

 x[226][0] = 0.001;
 x[226][1] = 0.5706;
 x[226][2] = -0.0248;
 x[226][3] = 1.2421;
 x[226][4] = -0.5621;
 y[226][0] = a;

 x[227][0] = 0.001;
 x[227][1] = -2.6479;
 x[227][2] = 10.1374;
 x[227][3] = -1.331;
 x[227][4] = -5.4707;
 y[227][0] = a;

 x[228][0] = 0.001;
 x[228][1] = 3.1219;
 x[228][2] = -3.137;
 x[228][3] = 1.9259;
 x[228][4] = -0.37458;
 y[228][0] = a;

 x[229][0] = 0.001;
 x[229][1] = 5.4944;
 x[229][2] = 1.5478;
 x[229][3] = 0.041694;
 x[229][4] = 1.9284;
 y[229][0] = a;

 x[230][0] = 0.001;
 x[230][1] = -1.3389;
 x[230][2] = 1.552;
 x[230][3] = 7.0806;
 x[230][4] = 1.031;
 y[230][0] = a;

 x[231][0] = 0.001;
 x[231][1] = -2.3361;
 x[231][2] = 11.9604;
 x[231][3] = 3.0835;
 x[231][4] = -5.4435;
 y[231][0] = a;

 x[232][0] = 0.001;
 x[232][1] = 2.2596;
 x[232][2] = -0.033118;
 x[232][3] = 4.7355;
 x[232][4] = -0.2776;
 y[232][0] = a;

 x[233][0] = 0.001;
 x[233][1] = 0.46901;
 x[233][2] = -0.63321;
 x[233][3] = 7.3848;
 x[233][4] = 0.36507;
 y[233][0] = a;

 x[234][0] = 0.001;
 x[234][1] = 2.7296;
 x[234][2] = 2.8701;
 x[234][3] = 0.51124;
 x[234][4] = 0.5099;
 y[234][0] = a;

 x[235][0] = 0.001;
 x[235][1] = 2.0466;
 x[235][2] = 2.03;
 x[235][3] = 2.1761;
 x[235][4] = -0.083634;
 y[235][0] = a;

 x[236][0] = 0.001;
 x[236][1] = -1.3274;
 x[236][2] = 9.498;
 x[236][3] = 2.4408;
 x[236][4] = -5.2689;
 y[236][0] = a;

 x[237][0] = 0.001;
 x[237][1] = 3.8905;
 x[237][2] = -2.1521;
 x[237][3] = 2.6302;
 x[237][4] = 1.1047;
 y[237][0] = a;

 x[238][0] = 0.001;
 x[238][1] = 3.9994;
 x[238][2] = 0.90427;
 x[238][3] = 1.1693;
 x[238][4] = 1.6892;
 y[238][0] = a;

 x[239][0] = 0.001;
 x[239][1] = 2.3952;
 x[239][2] = 9.5083;
 x[239][3] = -3.1783;
 x[239][4] = -3.0086;
 y[239][0] = a;

 x[240][0] = 0.001;
 x[240][1] = 3.2704;
 x[240][2] = 6.9321;
 x[240][3] = -1.0456;
 x[240][4] = 0.23447;
 y[240][0] = a;

 x[241][0] = 0.001;
 x[241][1] = -1.3931;
 x[241][2] = 1.5664;
 x[241][3] = 7.5382;
 x[241][4] = 0.78403;
 y[241][0] = a;

 x[242][0] = 0.001;
 x[242][1] = 1.6406;
 x[242][2] = 3.5488;
 x[242][3] = 1.3964;
 x[242][4] = -0.36424;
 y[242][0] = a;

 x[243][0] = 0.001;
 x[243][1] = 2.7744;
 x[243][2] = 6.8576;
 x[243][3] = -1.0671;
 x[243][4] = 0.075416;
 y[243][0] = a;

 x[244][0] = 0.001;
 x[244][1] = 2.4287;
 x[244][2] = 9.3821;
 x[244][3] = -3.2477;
 x[244][4] = -1.4543;
 y[244][0] = a;

 x[245][0] = 0.001;
 x[245][1] = 4.2134;
 x[245][2] = -2.806;
 x[245][3] = 2.0116;
 x[245][4] = 0.67412;
 y[245][0] = a;

 x[246][0] = 0.001;
 x[246][1] = 1.6472;
 x[246][2] = 0.48213;
 x[246][3] = 4.7449;
 x[246][4] = 1.225;
 y[246][0] = a;

 x[247][0] = 0.001;
 x[247][1] = 2.0597;
 x[247][2] = -0.99326;
 x[247][3] = 5.2119;
 x[247][4] = -0.29312;
 y[247][0] = a;

 x[248][0] = 0.001;
 x[248][1] = 0.3798;
 x[248][2] = 0.7098;
 x[248][3] = 0.7572;
 x[248][4] = -0.4444;
 y[248][0] = a;

 x[249][0] = 0.001;
 x[249][1] = 1.0135;
 x[249][2] = 8.4551;
 x[249][3] = -1.672;
 x[249][4] = -2.0815;
 y[249][0] = a;

 x[250][0] = 0.001;
 x[250][1] = 4.5691;
 x[250][2] = -4.4552;
 x[250][3] = 3.1769;
 x[250][4] = 0.004296;
 y[250][0] = a;

 x[251][0] = 0.001;
 x[251][1] = 0.57461;
 x[251][2] = 10.1105;
 x[251][3] = -1.6917;
 x[251][4] = -4.3922;
 y[251][0] = a;

 x[252][0] = 0.001;
 x[252][1] = 0.5734;
 x[252][2] = 9.1938;
 x[252][3] = -0.9094;
 x[252][4] = -1.872;
 y[252][0] = a;

 x[253][0] = 0.001;
 x[253][1] = 5.2868;
 x[253][2] = 3.257;
 x[253][3] = -1.3721;
 x[253][4] = 1.1668;
 y[253][0] = a;

 x[254][0] = 0.001;
 x[254][1] = 4.0102;
 x[254][2] = 10.6568;
 x[254][3] = -4.1388;
 x[254][4] = -5.0646;
 y[254][0] = a;

 x[255][0] = 0.001;
 x[255][1] = 4.1425;
 x[255][2] = -3.6792;
 x[255][3] = 3.8281;
 x[255][4] = 1.6297;
 y[255][0] = a;

 x[256][0] = 0.001;
 x[256][1] = 3.0934;
 x[256][2] = -2.9177;
 x[256][3] = 2.2232;
 x[256][4] = 0.22283;
 y[256][0] = a;

 x[257][0] = 0.001;
 x[257][1] = 2.2034;
 x[257][2] = 5.9947;
 x[257][3] = 0.53009;
 x[257][4] = 0.84998;
 y[257][0] = a;

 x[258][0] = 0.001;
 x[258][1] = 3.744;
 x[258][2] = 0.79459;
 x[258][3] = 0.95851;
 x[258][4] = 1.0077;
 y[258][0] = a;

 x[259][0] = 0.001;
 x[259][1] = 3.0329;
 x[259][2] = 2.2948;
 x[259][3] = 2.1135;
 x[259][4] = 0.35084;
 y[259][0] = a;

 x[260][0] = 0.001;
 x[260][1] = 3.7731;
 x[260][2] = 7.2073;
 x[260][3] = -1.6814;
 x[260][4] = -0.94742;
 y[260][0] = a;

 x[261][0] = 0.001;
 x[261][1] = 3.1557;
 x[261][2] = 2.8908;
 x[261][3] = 0.59693;
 x[261][4] = 0.79825;
 y[261][0] = a;

 x[262][0] = 0.001;
 x[262][1] = 1.8114;
 x[262][2] = 7.6067;
 x[262][3] = -0.9788;
 x[262][4] = -2.4668;
 y[262][0] = a;

 x[263][0] = 0.001;
 x[263][1] = 4.988;
 x[263][2] = 7.2052;
 x[263][3] = -3.2846;
 x[263][4] = -1.1608;
 y[263][0] = a;

 x[264][0] = 0.001;
 x[264][1] = 2.483;
 x[264][2] = 6.6155;
 x[264][3] = -0.79287;
 x[264][4] = -0.90863;
 y[264][0] = a;

 x[265][0] = 0.001;
 x[265][1] = 1.594;
 x[265][2] = 4.7055;
 x[265][3] = 1.3758;
 x[265][4] = 0.081882;
 y[265][0] = a;

 x[266][0] = 0.001;
 x[266][1] = -0.016103;
 x[266][2] = 9.7484;
 x[266][3] = 0.15394;
 x[266][4] = -1.6134;
 y[266][0] = a;

 x[267][0] = 0.001;
 x[267][1] = 3.8496;
 x[267][2] = 9.7939;
 x[267][3] = -4.1508;
 x[267][4] = -4.4582;
 y[267][0] = a;

 x[268][0] = 0.001;
 x[268][1] = 0.9297;
 x[268][2] = -3.7971;
 x[268][3] = 4.6429;
 x[268][4] = -0.2957;
 y[268][0] = a;

 x[269][0] = 0.001;
 x[269][1] = 4.9342;
 x[269][2] = 2.4107;
 x[269][3] = -0.17594;
 x[269][4] = 1.6245;
 y[269][0] = a;

 x[270][0] = 0.001;
 x[270][1] = 3.8417;
 x[270][2] = 10.0215;
 x[270][3] = -4.2699;
 x[270][4] = -4.9159;
 y[270][0] = a;

 x[271][0] = 0.001;
 x[271][1] = 5.3915;
 x[271][2] = 9.9946;
 x[271][3] = -3.8081;
 x[271][4] = -3.3642;
 y[271][0] = a;

 x[272][0] = 0.001;
 x[272][1] = 4.4072;
 x[272][2] = -0.070365;
 x[272][3] = 2.0416;
 x[272][4] = 1.1319;
 y[272][0] = a;

 x[273][0] = 0.001;
 x[273][1] = 2.6946;
 x[273][2] = 6.7976;
 x[273][3] = -0.40301;
 x[273][4] = 0.44912;
 y[273][0] = a;

 x[274][0] = 0.001;
 x[274][1] = 5.2756;
 x[274][2] = 0.13863;
 x[274][3] = 0.12138;
 x[274][4] = 1.1435;
 y[274][0] = a;

 x[275][0] = 0.001;
 x[275][1] = 3.4312;
 x[275][2] = 6.2637;
 x[275][3] = -1.9513;
 x[275][4] = -0.36165;
 y[275][0] = a;

 x[276][0] = 0.001;
 x[276][1] = 4.052;
 x[276][2] = -0.16555;
 x[276][3] = 0.45383;
 x[276][4] = 0.51248;
 y[276][0] = a;

 x[277][0] = 0.001;
 x[277][1] = 1.3638;
 x[277][2] = -4.7759;
 x[277][3] = 8.4182;
 x[277][4] = -1.8836;
 y[277][0] = a;

 x[278][0] = 0.001;
 x[278][1] = 0.89566;
 x[278][2] = 7.7763;
 x[278][3] = -2.7473;
 x[278][4] = -1.9353;
 y[278][0] = a;

 x[279][0] = 0.001;
 x[279][1] = 1.9265;
 x[279][2] = 7.7557;
 x[279][3] = -0.16823;
 x[279][4] = -3.0771;
 y[279][0] = a;

 x[280][0] = 0.001;
 x[280][1] = 0.20977;
 x[280][2] = -0.46146;
 x[280][3] = 7.7267;
 x[280][4] = 0.90946;
 y[280][0] = a;

 x[281][0] = 0.001;
 x[281][1] = 4.068;
 x[281][2] = -2.9363;
 x[281][3] = 2.1992;
 x[281][4] = 0.50084;
 y[281][0] = a;

 x[282][0] = 0.001;
 x[282][1] = 2.877;
 x[282][2] = -4.0599;
 x[282][3] = 3.6259;
 x[282][4] = -0.32544;
 y[282][0] = a;

 x[283][0] = 0.001;
 x[283][1] = 0.3223;
 x[283][2] = -0.89808;
 x[283][3] = 8.0883;
 x[283][4] = 0.69222;
 y[283][0] = a;

 x[284][0] = 0.001;
 x[284][1] = -1.3;
 x[284][2] = 10.2678;
 x[284][3] = -2.953;
 x[284][4] = -5.8638;
 y[284][0] = a;

 x[285][0] = 0.001;
 x[285][1] = 1.7747;
 x[285][2] = -6.4334;
 x[285][3] = 8.15;
 x[285][4] = -0.89828;
 y[285][0] = a;

 x[286][0] = 0.001;
 x[286][1] = 1.3419;
 x[286][2] = -4.4221;
 x[286][3] = 8.09;
 x[286][4] = -1.7349;
 y[286][0] = a;

 x[287][0] = 0.001;
 x[287][1] = 0.89606;
 x[287][2] = 10.5471;
 x[287][3] = -1.4175;
 x[287][4] = -4.0327;
 y[287][0] = a;

 x[288][0] = 0.001;
 x[288][1] = 0.44125;
 x[288][2] = 2.9487;
 x[288][3] = 4.3225;
 x[288][4] = 0.7155;
 y[288][0] = a;

 x[289][0] = 0.001;
 x[289][1] = 3.2422;
 x[289][2] = 6.2265;
 x[289][3] = 0.12224;
 x[289][4] = -1.4466;
 y[289][0] = a;

 x[290][0] = 0.001;
 x[290][1] = 2.5678;
 x[290][2] = 3.5136;
 x[290][3] = 0.61406;
 x[290][4] = -0.40691;
 y[290][0] = a;

 x[291][0] = 0.001;
 x[291][1] = -2.2153;
 x[291][2] = 11.9625;
 x[291][3] = 0.078538;
 x[291][4] = -7.7853;
 y[291][0] = a;

 x[292][0] = 0.001;
 x[292][1] = 4.1349;
 x[292][2] = 6.1189;
 x[292][3] = -2.4294;
 x[292][4] = -0.19613;
 y[292][0] = a;

 x[293][0] = 0.001;
 x[293][1] = 1.934;
 x[293][2] = -0.000009;
 x[293][3] = 4.816;
 x[293][4] = -0.33967;
 y[293][0] = a;

 x[294][0] = 0.001;
 x[294][1] = 2.5068;
 x[294][2] = 1.1588;
 x[294][3] = 3.9249;
 x[294][4] = 0.12585;
 y[294][0] = a;

 x[295][0] = 0.001;
 x[295][1] = 2.1464;
 x[295][2] = 6.0795;
 x[295][3] = -0.5778;
 x[295][4] = -2.2302;
 y[295][0] = a;

 x[296][0] = 0.001;
 x[296][1] = 0.051979;
 x[296][2] = 7.0521;
 x[296][3] = -2.0541;
 x[296][4] = -3.1508;
 y[296][0] = a;

 x[297][0] = 0.001;
 x[297][1] = 1.2706;
 x[297][2] = 8.035;
 x[297][3] = -0.19651;
 x[297][4] = -2.1888;
 y[297][0] = a;

 x[298][0] = 0.001;
 x[298][1] = 1.143;
 x[298][2] = 0.83391;
 x[298][3] = 5.4552;
 x[298][4] = -0.56984;
 y[298][0] = a;

 x[299][0] = 0.001;
 x[299][1] = 2.2928;
 x[299][2] = 9.0386;
 x[299][3] = -3.2417;
 x[299][4] = -1.2991;
 y[299][0] = a;

 x[300][0] = 0.001;
 x[300][1] = 0.3292;
 x[300][2] = -4.4552;
 x[300][3] = 4.5718;
 x[300][4] = -0.9888;
 y[300][0] = a;

 x[301][0] = 0.001;
 x[301][1] = -0.10648;
 x[301][2] = -0.76771;
 x[301][3] = 7.7575;
 x[301][4] = 0.64179;
 y[301][0] = a;

 x[302][0] = 0.001;
 x[302][1] = 0.72252;
 x[302][2] = -0.053811;
 x[302][3] = 5.6703;
 x[302][4] = -1.3509;
 y[302][0] = a;

 x[303][0] = 0.001;
 x[303][1] = 4.2475;
 x[303][2] = 1.4816;
 x[303][3] = -0.48355;
 x[303][4] = 0.95343;
 y[303][0] = a;

 x[304][0] = 0.001;
 x[304][1] = 3.9772;
 x[304][2] = 0.33521;
 x[304][3] = 2.2566;
 x[304][4] = 2.1625;
 y[304][0] = a;

 x[305][0] = 0.001;
 x[305][1] = 3.6667;
 x[305][2] = 4.302;
 x[305][3] = 0.55923;
 x[305][4] = 0.33791;
 y[305][0] = a;

 x[306][0] = 0.001;
 x[306][1] = 2.8232;
 x[306][2] = 10.8513;
 x[306][3] = -3.1466;
 x[306][4] = -3.9784;
 y[306][0] = a;

 x[307][0] = 0.001;
 x[307][1] = -1.4217;
 x[307][2] = 11.6542;
 x[307][3] = -0.057699;
 x[307][4] = -7.1025;
 y[307][0] = a;

 x[308][0] = 0.001;
 x[308][1] = 4.2458;
 x[308][2] = 1.1981;
 x[308][3] = 0.66633;
 x[308][4] = 0.94696;
 y[308][0] = a;

 x[309][0] = 0.001;
 x[309][1] = 4.1038;
 x[309][2] = -4.8069;
 x[309][3] = 3.3491;
 x[309][4] = -0.49225;
 y[309][0] = a;

 x[310][0] = 0.001;
 x[310][1] = 1.4507;
 x[310][2] = 8.7903;
 x[310][3] = -2.2324;
 x[310][4] = -0.65259;
 y[310][0] = a;

 x[311][0] = 0.001;
 x[311][1] = 3.4647;
 x[311][2] = -3.9172;
 x[311][3] = 3.9746;
 x[311][4] = 0.36119;
 y[311][0] = a;

 x[312][0] = 0.001;
 x[312][1] = 1.8533;
 x[312][2] = 6.1458;
 x[312][3] = 1.0176;
 x[312][4] = -2.0401;
 y[312][0] = a;

 x[313][0] = 0.001;
 x[313][1] = 3.5288;
 x[313][2] = 0.71596;
 x[313][3] = 1.9507;
 x[313][4] = 1.9375;
 y[313][0] = a;

 x[314][0] = 0.001;
 x[314][1] = 3.9719;
 x[314][2] = 1.0367;
 x[314][3] = 0.75973;
 x[314][4] = 1.0013;
 y[314][0] = a;

 x[315][0] = 0.001;
 x[315][1] = 3.534;
 x[315][2] = 9.3614;
 x[315][3] = -3.6316;
 x[315][4] = -1.2461;
 y[315][0] = a;

 x[316][0] = 0.001;
 x[316][1] = 3.6894;
 x[316][2] = 9.887;
 x[316][3] = -4.0788;
 x[316][4] = -4.3664;
 y[316][0] = a;

 x[317][0] = 0.001;
 x[317][1] = 3.0672;
 x[317][2] = -4.4117;
 x[317][3] = 3.8238;
 x[317][4] = -0.81682;
 y[317][0] = a;

 x[318][0] = 0.001;
 x[318][1] = 2.6463;
 x[318][2] = -4.8152;
 x[318][3] = 6.3549;
 x[318][4] = 0.003003;
 y[318][0] = a;

 x[319][0] = 0.001;
 x[319][1] = 2.2893;
 x[319][2] = 3.733;
 x[319][3] = 0.6312;
 x[319][4] = -0.39786;
 y[319][0] = a;

 x[320][0] = 0.001;
 x[320][1] = 1.5673;
 x[320][2] = 7.9274;
 x[320][3] = -0.056842;
 x[320][4] = -2.1694;
 y[320][0] = a;

 x[321][0] = 0.001;
 x[321][1] = 4.0405;
 x[321][2] = 0.51524;
 x[321][3] = 1.0279;
 x[321][4] = 1.106;
 y[321][0] = a;

 x[322][0] = 0.001;
 x[322][1] = 4.3846;
 x[322][2] = -4.8794;
 x[322][3] = 3.3662;
 x[322][4] = -0.029324;
 y[322][0] = a;

 x[323][0] = 0.001;
 x[323][1] = 2.0165;
 x[323][2] = -0.25246;
 x[323][3] = 5.1707;
 x[323][4] = 1.0763;
 y[323][0] = a;

 x[324][0] = 0.001;
 x[324][1] = 4.0446;
 x[324][2] = 11.1741;
 x[324][3] = -4.3582;
 x[324][4] = -4.7401;
 y[324][0] = a;

 x[325][0] = 0.001;
 x[325][1] = -0.33729;
 x[325][2] = -0.64976;
 x[325][3] = 7.6659;
 x[325][4] = 0.72326;
 y[325][0] = a;

 x[326][0] = 0.001;
 x[326][1] = -2.4604;
 x[326][2] = 12.7302;
 x[326][3] = 0.91738;
 x[326][4] = -7.6418;
 y[326][0] = a;

 x[327][0] = 0.001;
 x[327][1] = 4.1195;
 x[327][2] = 10.9258;
 x[327][3] = -3.8929;
 x[327][4] = -4.1802;
 y[327][0] = a;

 x[328][0] = 0.001;
 x[328][1] = 2.0193;
 x[328][2] = 0.82356;
 x[328][3] = 4.6369;
 x[328][4] = 1.4202;
 y[328][0] = a;

 x[329][0] = 0.001;
 x[329][1] = 1.5701;
 x[329][2] = 7.9129;
 x[329][3] = 0.29018;
 x[329][4] = -2.1953;
 y[329][0] = a;

 x[330][0] = 0.001;
 x[330][1] = 2.6415;
 x[330][2] = 7.586;
 x[330][3] = -0.28562;
 x[330][4] = -1.6677;
 y[330][0] = a;

 x[331][0] = 0.001;
 x[331][1] = 5.0214;
 x[331][2] = 8.0764;
 x[331][3] = -3.0515;
 x[331][4] = -1.7155;
 y[331][0] = a;

 x[332][0] = 0.001;
 x[332][1] = 4.3435;
 x[332][2] = 3.3295;
 x[332][3] = 0.83598;
 x[332][4] = 0.64955;
 y[332][0] = a;

 x[333][0] = 0.001;
 x[333][1] = 1.8238;
 x[333][2] = -6.7748;
 x[333][3] = 8.3873;
 x[333][4] = -0.54139;
 y[333][0] = a;

 x[334][0] = 0.001;
 x[334][1] = 3.9382;
 x[334][2] = 0.9291;
 x[334][3] = 0.78543;
 x[334][4] = 0.6767;
 y[334][0] = a;

 x[335][0] = 0.001;
 x[335][1] = 2.2517;
 x[335][2] = -5.1422;
 x[335][3] = 4.2916;
 x[335][4] = -1.2487;
 y[335][0] = a;

 x[336][0] = 0.001;
 x[336][1] = 5.504;
 x[336][2] = 10.3671;
 x[336][3] = -4.413;
 x[336][4] = -4.0211;
 y[336][0] = a;

 x[337][0] = 0.001;
 x[337][1] = 2.8521;
 x[337][2] = 9.171;
 x[337][3] = -3.6461;
 x[337][4] = -1.2047;
 y[337][0] = a;

 x[338][0] = 0.001;
 x[338][1] = 1.1676;
 x[338][2] = 9.1566;
 x[338][3] = -2.0867;
 x[338][4] = -0.80647;
 y[338][0] = a;

 x[339][0] = 0.001;
 x[339][1] = 2.6104;
 x[339][2] = 8.0081;
 x[339][3] = -0.23592;
 x[339][4] = -1.7608;
 y[339][0] = a;

 x[340][0] = 0.001;
 x[340][1] = 0.32444;
 x[340][2] = 10.067;
 x[340][3] = -1.1982;
 x[340][4] = -4.1284;
 y[340][0] = a;

 x[341][0] = 0.001;
 x[341][1] = 3.8962;
 x[341][2] = -4.7904;
 x[341][3] = 3.3954;
 x[341][4] = -0.53751;
 y[341][0] = a;

 x[342][0] = 0.001;
 x[342][1] = 2.1752;
 x[342][2] = -0.8091;
 x[342][3] = 5.1022;
 x[342][4] = -0.67975;
 y[342][0] = a;

 x[343][0] = 0.001;
 x[343][1] = 1.1588;
 x[343][2] = 8.9331;
 x[343][3] = -2.0807;
 x[343][4] = -1.1272;
 y[343][0] = a;

 x[344][0] = 0.001;
 x[344][1] = 4.7072;
 x[344][2] = 8.2957;
 x[344][3] = -2.5605;
 x[344][4] = -1.4905;
 y[344][0] = a;

 x[345][0] = 0.001;
 x[345][1] = -1.9667;
 x[345][2] = 11.8052;
 x[345][3] = -0.40472;
 x[345][4] = -7.8719;
 y[345][0] = a;

 x[346][0] = 0.001;
 x[346][1] = 4.0552;
 x[346][2] = 0.40143;
 x[346][3] = 1.4563;
 x[346][4] = 0.65343;
 y[346][0] = a;

 x[347][0] = 0.001;
 x[347][1] = 2.3678;
 x[347][2] = -6.839;
 x[347][3] = 8.4207;
 x[347][4] = -0.44829;
 y[347][0] = a;

 x[348][0] = 0.001;
 x[348][1] = 0.33565;
 x[348][2] = 6.8369;
 x[348][3] = 0.69718;
 x[348][4] = -0.55691;
 y[348][0] = a;

 x[349][0] = 0.001;
 x[349][1] = 4.3398;
 x[349][2] = -5.3036;
 x[349][3] = 3.8803;
 x[349][4] = -0.70432;
 y[349][0] = a;

 x[350][0] = 0.001;
 x[350][1] = 1.5456;
 x[350][2] = 8.5482;
 x[350][3] = 0.4187;
 x[350][4] = -2.1784;
 y[350][0] = a;

 x[351][0] = 0.001;
 x[351][1] = 1.4276;
 x[351][2] = 8.3847;
 x[351][3] = -2.0995;
 x[351][4] = -1.9677;
 y[351][0] = a;

 x[352][0] = 0.001;
 x[352][1] = -0.27802;
 x[352][2] = 8.1881;
 x[352][3] = -3.1338;
 x[352][4] = -2.5276;
 y[352][0] = a;

 x[353][0] = 0.001;
 x[353][1] = 0.93611;
 x[353][2] = 8.6413;
 x[353][3] = -1.6351;
 x[353][4] = -1.3043;
 y[353][0] = a;

 x[354][0] = 0.001;
 x[354][1] = 4.6352;
 x[354][2] = -3.0087;
 x[354][3] = 2.6773;
 x[354][4] = 1.212;
 y[354][0] = a;

 x[355][0] = 0.001;
 x[355][1] = 1.5268;
 x[355][2] = -5.5871;
 x[355][3] = 8.6564;
 x[355][4] = -1.722;
 y[355][0] = a;

 x[356][0] = 0.001;
 x[356][1] = 0.95626;
 x[356][2] = 2.4728;
 x[356][3] = 4.4578;
 x[356][4] = 0.21636;
 y[356][0] = a;

 x[357][0] = 0.001;
 x[357][1] = -2.7914;
 x[357][2] = 1.7734;
 x[357][3] = 6.7756;
 x[357][4] = -0.39915;
 y[357][0] = a;

 x[358][0] = 0.001;
 x[358][1] = 5.2032;
 x[358][2] = 3.5116;
 x[358][3] = -1.2538;
 x[358][4] = 1.0129;
 y[358][0] = a;

 x[359][0] = 0.001;
 x[359][1] = 3.1836;
 x[359][2] = 7.2321;
 x[359][3] = -1.0713;
 x[359][4] = -2.5909;
 y[359][0] = a;

 x[360][0] = 0.001;
 x[360][1] = 0.65497;
 x[360][2] = 5.1815;
 x[360][3] = 1.0673;
 x[360][4] = -0.42113;
 y[360][0] = a;

 x[361][0] = 0.001;
 x[361][1] = 5.6084;
 x[361][2] = 10.3009;
 x[361][3] = -4.8003;
 x[361][4] = -4.3534;
 y[361][0] = a;

 x[362][0] = 0.001;
 x[362][1] = 1.105;
 x[362][2] = 7.4432;
 x[362][3] = 0.41099;
 x[362][4] = -3.0332;
 y[362][0] = a;

 x[363][0] = 0.001;
 x[363][1] = 3.9292;
 x[363][2] = -2.9156;
 x[363][3] = 2.2129;
 x[363][4] = 0.30817;
 y[363][0] = a;

 x[364][0] = 0.001;
 x[364][1] = 1.1558;
 x[364][2] = 6.4003;
 x[364][3] = 1.5506;
 x[364][4] = 0.6961;
 y[364][0] = a;

 x[365][0] = 0.001;
 x[365][1] = 2.5581;
 x[365][2] = 2.6218;
 x[365][3] = 1.8513;
 x[365][4] = 0.40257;
 y[365][0] = a;

 x[366][0] = 0.001;
 x[366][1] = 2.7831;
 x[366][2] = 10.9796;
 x[366][3] = -3.557;
 x[366][4] = -4.4039;
 y[366][0] = a;

 x[367][0] = 0.001;
 x[367][1] = 3.7635;
 x[367][2] = 2.7811;
 x[367][3] = 0.66119;
 x[367][4] = 0.34179;
 y[367][0] = a;

 x[368][0] = 0.001;
 x[368][1] = -2.6479;
 x[368][2] = 10.1374;
 x[368][3] = -1.331;
 x[368][4] = -5.4707;
 y[368][0] = a;

 x[369][0] = 0.001;
 x[369][1] = 1.0652;
 x[369][2] = 8.3682;
 x[369][3] = -1.4004;
 x[369][4] = -1.6509;
 y[369][0] = a;

 x[370][0] = 0.001;
 x[370][1] = -1.4275;
 x[370][2] = 11.8797;
 x[370][3] = 0.41613;
 x[370][4] = -6.9978;
 y[370][0] = a;

 x[371][0] = 0.001;
 x[371][1] = 5.7456;
 x[371][2] = 10.1808;
 x[371][3] = -4.7857;
 x[371][4] = -4.3366;
 y[371][0] = a;

 x[372][0] = 0.001;
 x[372][1] = 5.086;
 x[372][2] = 3.2798;
 x[372][3] = -1.2701;
 x[372][4] = 1.1189;
 y[372][0] = a;

 x[373][0] = 0.001;
 x[373][1] = 3.4092;
 x[373][2] = 5.4049;
 x[373][3] = -2.5228;
 x[373][4] = -0.89958;
 y[373][0] = a;

 x[374][0] = 0.001;
 x[374][1] = -0.2361;
 x[374][2] = 9.3221;
 x[374][3] = 2.1307;
 x[374][4] = -4.3793;
 y[374][0] = a;

 x[375][0] = 0.001;
 x[375][1] = 3.8197;
 x[375][2] = 8.9951;
 x[375][3] = -4.383;
 x[375][4] = -4.0327;
 y[375][0] = a;

 x[376][0] = 0.001;
 x[376][1] = -1.1391;
 x[376][2] = 1.8127;
 x[376][3] = 6.9144;
 x[376][4] = 0.70127;
 y[376][0] = a;

 x[377][0] = 0.001;
 x[377][1] = 4.9249;
 x[377][2] = 0.68906;
 x[377][3] = 0.77344;
 x[377][4] = 1.2095;
 y[377][0] = a;

 x[378][0] = 0.001;
 x[378][1] = 2.5089;
 x[378][2] = 6.841;
 x[378][3] = -0.029423;
 x[378][4] = 0.44912;
 y[378][0] = a;

 x[379][0] = 0.001;
 x[379][1] = -0.2062;
 x[379][2] = 9.2207;
 x[379][3] = -3.7044;
 x[379][4] = -6.8103;
 y[379][0] = a;

 x[380][0] = 0.001;
 x[380][1] = 3.946;
 x[380][2] = 6.8514;
 x[380][3] = -1.5443;
 x[380][4] = -0.5582;
 y[380][0] = a;

 x[381][0] = 0.001;
 x[381][1] = -0.278;
 x[381][2] = 8.1881;
 x[381][3] = -3.1338;
 x[381][4] = -2.5276;
 y[381][0] = a;

 x[382][0] = 0.001;
 x[382][1] = 1.8592;
 x[382][2] = 3.2074;
 x[382][3] = -0.15966;
 x[382][4] = -0.26208;
 y[382][0] = a;

 x[383][0] = 0.001;
 x[383][1] = 0.56953;
 x[383][2] = 7.6294;
 x[383][3] = 1.5754;
 x[383][4] = -3.2233;
 y[383][0] = a;

 x[384][0] = 0.001;
 x[384][1] = 3.4626;
 x[384][2] = -4.449;
 x[384][3] = 3.5427;
 x[384][4] = 0.15429;
 y[384][0] = a;

 x[385][0] = 0.001;
 x[385][1] = 3.3951;
 x[385][2] = 1.1484;
 x[385][3] = 2.1401;
 x[385][4] = 2.0862;
 y[385][0] = a;

 x[386][0] = 0.001;
 x[386][1] = 5.0429;
 x[386][2] = -0.52974;
 x[386][3] = 0.50439;
 x[386][4] = 1.106;
 y[386][0] = a;

 x[387][0] = 0.001;
 x[387][1] = 3.7758;
 x[387][2] = 7.1783;
 x[387][3] = -1.5195;
 x[387][4] = 0.40128;
 y[387][0] = a;

 x[388][0] = 0.001;
 x[388][1] = 4.6562;
 x[388][2] = 7.6398;
 x[388][3] = -2.4243;
 x[388][4] = -1.2384;
 y[388][0] = a;

 x[389][0] = 0.001;
 x[389][1] = 4.0948;
 x[389][2] = -2.9674;
 x[389][3] = 2.3689;
 x[389][4] = 0.75429;
 y[389][0] = a;

 x[390][0] = 0.001;
 x[390][1] = 1.8384;
 x[390][2] = 6.063;
 x[390][3] = 0.54723;
 x[390][4] = 0.51248;
 y[390][0] = a;

 x[391][0] = 0.001;
 x[391][1] = 2.0153;
 x[391][2] = 0.43661;
 x[391][3] = 4.5864;
 x[391][4] = -0.3151;
 y[391][0] = a;

 x[392][0] = 0.001;
 x[392][1] = 3.5251;
 x[392][2] = 0.7201;
 x[392][3] = 1.6928;
 x[392][4] = 0.64438;
 y[392][0] = a;

 x[393][0] = 0.001;
 x[393][1] = 3.757;
 x[393][2] = -5.4236;
 x[393][3] = 3.8255;
 x[393][4] = -1.2526;
 y[393][0] = a;

 x[394][0] = 0.001;
 x[394][1] = 2.5989;
 x[394][2] = 3.5178;
 x[394][3] = 0.7623;
 x[394][4] = 0.81119;
 y[394][0] = a;

 x[395][0] = 0.001;
 x[395][1] = 1.8994;
 x[395][2] = 0.97462;
 x[395][3] = 4.2265;
 x[395][4] = 0.81377;
 y[395][0] = a;

 x[396][0] = 0.001;
 x[396][1] = 3.6941;
 x[396][2] = -3.9482;
 x[396][3] = 4.2625;
 x[396][4] = 1.1577;
 y[396][0] = a;

 x[397][0] = 0.001;
 x[397][1] = 4.4295;
 x[397][2] = -2.3507;
 x[397][3] = 1.7048;
 x[397][4] = 0.90946;
 y[397][0] = a;

 x[398][0] = 0.001;
 x[398][1] = 6.8248;
 x[398][2] = 5.2187;
 x[398][3] = -2.5425;
 x[398][4] = 0.5461;
 y[398][0] = a;

 x[399][0] = 0.001;
 x[399][1] = 1.8967;
 x[399][2] = -2.5163;
 x[399][3] = 2.8093;
 x[399][4] = -0.79742;
 y[399][0] = a;

 x[400][0] = 0.001;
 x[400][1] = 2.1526;
 x[400][2] = -6.1665;
 x[400][3] = 8.0831;
 x[400][4] = -0.34355;
 y[400][0] = a;

 x[401][0] = 0.001;
 x[401][1] = 3.3004;
 x[401][2] = 7.0811;
 x[401][3] = -1.3258;
 x[401][4] = 0.22283;
 y[401][0] = a;

 x[402][0] = 0.001;
 x[402][1] = 2.7213;
 x[402][2] = 7.05;
 x[402][3] = -0.58808;
 x[402][4] = 0.41809;
 y[402][0] = a;

 x[403][0] = 0.001;
 x[403][1] = 3.8846;
 x[403][2] = -3.0336;
 x[403][3] = 2.5334;
 x[403][4] = 0.20214;
 y[403][0] = a;

 x[404][0] = 0.001;
 x[404][1] = 4.1665;
 x[404][2] = -0.4449;
 x[404][3] = 0.23448;
 x[404][4] = 0.27843;
 y[404][0] = a;

 x[405][0] = 0.001;
 x[405][1] = 0.94225;
 x[405][2] = 5.8561;
 x[405][3] = 1.8762;
 x[405][4] = -0.32544;
 y[405][0] = a;

 x[406][0] = 0.001;
 x[406][1] = 5.1321;
 x[406][2] = -0.031048;
 x[406][3] = 0.32616;
 x[406][4] = 1.1151;
 y[406][0] = a;

 x[407][0] = 0.001;
 x[407][1] = 0.38251;
 x[407][2] = 6.8121;
 x[407][3] = 1.8128;
 x[407][4] = -0.61251;
 y[407][0] = a;

 x[408][0] = 0.001;
 x[408][1] = 3.0333;
 x[408][2] = -2.5928;
 x[408][3] = 2.3183;
 x[408][4] = 0.303;
 y[408][0] = a;

 x[409][0] = 0.001;
 x[409][1] = 2.9233;
 x[409][2] = 6.0464;
 x[409][3] = -0.11168;
 x[409][4] = -0.58665;
 y[409][0] = a;

 x[410][0] = 0.001;
 x[410][1] = 1.162;
 x[410][2] = 10.2926;
 x[410][3] = -1.2821;
 x[410][4] = -4.0392;
 y[410][0] = a;

 x[411][0] = 0.001;
 x[411][1] = 3.7791;
 x[411][2] = 2.5762;
 x[411][3] = 1.3098;
 x[411][4] = 0.5655;
 y[411][0] = a;

 x[412][0] = 0.001;
 x[412][1] = 0.77765;
 x[412][2] = 5.9781;
 x[412][3] = 1.1941;
 x[412][4] = -0.3526;
 y[412][0] = a;

 x[413][0] = 0.001;
 x[413][1] = -0.38388;
 x[413][2] = -1.0471;
 x[413][3] = 8.0514;
 x[413][4] = 0.49567;
 y[413][0] = a;

 x[414][0] = 0.001;
 x[414][1] = 0.21084;
 x[414][2] = 9.4359;
 x[414][3] = -0.094543;
 x[414][4] = -1.859;
 y[414][0] = a;

 x[415][0] = 0.001;
 x[415][1] = 2.9571;
 x[415][2] = -4.5938;
 x[415][3] = 5.9068;
 x[415][4] = 0.57196;
 y[415][0] = a;

 x[416][0] = 0.001;
 x[416][1] = 4.6439;
 x[416][2] = -3.3729;
 x[416][3] = 2.5976;
 x[416][4] = 0.55257;
 y[416][0] = a;

 x[417][0] = 0.001;
 x[417][1] = 3.3577;
 x[417][2] = -4.3062;
 x[417][3] = 6.0241;
 x[417][4] = 0.18274;
 y[417][0] = a;

 x[418][0] = 0.001;
 x[418][1] = 3.5127;
 x[418][2] = 2.9073;
 x[418][3] = 1.0579;
 x[418][4] = 0.40774;
 y[418][0] = a;

 x[419][0] = 0.001;
 x[419][1] = 2.6562;
 x[419][2] = 10.7044;
 x[419][3] = -3.3085;
 x[419][4] = -4.0767;
 y[419][0] = a;

 x[420][0] = 0.001;
 x[420][1] = -1.3612;
 x[420][2] = 10.694;
 x[420][3] = 1.7022;
 x[420][4] = -2.9026;
 y[420][0] = a;

 x[421][0] = 0.001;
 x[421][1] = -0.278;
 x[421][2] = 8.1881;
 x[421][3] = -3.1338;
 x[421][4] = -2.5276;
 y[421][0] = a;

 x[422][0] = 0.001;
 x[422][1] = 1.04;
 x[422][2] = -6.9321;
 x[422][3] = 8.2888;
 x[422][4] = -1.2991;
 y[422][0] = a;

 x[423][0] = 0.001;
 x[423][1] = 2.1881;
 x[423][2] = 2.7356;
 x[423][3] = 1.3278;
 x[423][4] = -0.1832;
 y[423][0] = a;

 x[424][0] = 0.001;
 x[424][1] = 4.2756;
 x[424][2] = -2.6528;
 x[424][3] = 2.1375;
 x[424][4] = 0.94437;
 y[424][0] = a;

 x[425][0] = 0.001;
 x[425][1] = -0.11996;
 x[425][2] = 6.8741;
 x[425][3] = 0.91995;
 x[425][4] = -0.6694;
 y[425][0] = a;

 x[426][0] = 0.001;
 x[426][1] = 2.9736;
 x[426][2] = 8.7944;
 x[426][3] = -3.6359;
 x[426][4] = -1.3754;
 y[426][0] = a;

 x[427][0] = 0.001;
 x[427][1] = 3.7798;
 x[427][2] = -3.3109;
 x[427][3] = 2.6491;
 x[427][4] = 0.066365;
 y[427][0] = a;

 x[428][0] = 0.001;
 x[428][1] = 5.3586;
 x[428][2] = 3.7557;
 x[428][3] = -1.7345;
 x[428][4] = 1.0789;
 y[428][0] = a;

 x[429][0] = 0.001;
 x[429][1] = 1.8373;
 x[429][2] = 6.1292;
 x[429][3] = 0.84027;
 x[429][4] = 0.55257;
 y[429][0] = a;

 x[430][0] = 0.001;
 x[430][1] = 1.2262;
 x[430][2] = 0.89599;
 x[430][3] = 5.7568;
 x[430][4] = -0.11596;
 y[430][0] = a;

 x[431][0] = 0.001;
 x[431][1] = -0.048008;
 x[431][2] = -0.56078;
 x[431][3] = 7.7215;
 x[431][4] = 0.453;
 y[431][0] = a;

 x[432][0] = 0.001;
 x[432][1] = 0.5706;
 x[432][2] = -0.024841;
 x[432][3] = 1.2421;
 x[432][4] = -0.56208;
 y[432][0] = a;

 x[433][0] = 0.001;
 x[433][1] = 4.3634;
 x[433][2] = 0.46351;
 x[433][3] = 1.4281;
 x[433][4] = 2.0202;
 y[433][0] = a;

 x[434][0] = 0.001;
 x[434][1] = 3.482;
 x[434][2] = -4.1634;
 x[434][3] = 3.5008;
 x[434][4] = -0.078462;
 y[434][0] = a;

 x[435][0] = 0.001;
 x[435][1] = 0.51947;
 x[435][2] = -3.2633;
 x[435][3] = 3.0895;
 x[435][4] = -0.98492;
 y[435][0] = a;

 x[436][0] = 0.001;
 x[436][1] = 2.3164;
 x[436][2] = -2.628;
 x[436][3] = 3.1529;
 x[436][4] = -0.08622;
 y[436][0] = a;

 x[437][0] = 0.001;
 x[437][1] = -1.8348;
 x[437][2] = 11.0334;
 x[437][3] = 3.1863;
 x[437][4] = -4.8888;
 y[437][0] = a;

 x[438][0] = 0.001;
 x[438][1] = 1.3754;
 x[438][2] = 8.8793;
 x[438][3] = -1.9136;
 x[438][4] = -0.53751;
 y[438][0] = a;

 x[439][0] = 0.001;
 x[439][1] = -0.16682;
 x[439][2] = 5.8974;
 x[439][3] = 0.49839;
 x[439][4] = -0.70044;
 y[439][0] = a;

 x[440][0] = 0.001;
 x[440][1] = 0.29961;
 x[440][2] = 7.1328;
 x[440][3] = -0.31475;
 x[440][4] = -1.1828;
 y[440][0] = a;

 x[441][0] = 0.001;
 x[441][1] = 0.25035;
 x[441][2] = 9.3262;
 x[441][3] = -3.6873;
 x[441][4] = -6.2543;
 y[441][0] = a;

 x[442][0] = 0.001;
 x[442][1] = 2.4673;
 x[442][2] = 1.3926;
 x[442][3] = 1.7125;
 x[442][4] = 0.41421;
 y[442][0] = a;

 x[443][0] = 0.001;
 x[443][1] = 0.77805;
 x[443][2] = 6.6424;
 x[443][3] = -1.1425;
 x[443][4] = -1.0573;
 y[443][0] = a;

 x[444][0] = 0.001;
 x[444][1] = 3.4465;
 x[444][2] = 2.9508;
 x[444][3] = 1.0271;
 x[444][4] = 0.5461;
 y[444][0] = a;

 x[445][0] = 0.001;
 x[445][1] = 2.2429;
 x[445][2] = -4.1427;
 x[445][3] = 5.2333;
 x[445][4] = -0.40173;
 y[445][0] = a;

 x[446][0] = 0.001;
 x[446][1] = 3.7321;
 x[446][2] = -3.884;
 x[446][3] = 3.3577;
 x[446][4] = -0.006049;
 y[446][0] = a;

 x[447][0] = 0.001;
 x[447][1] = 4.3365;
 x[447][2] = -3.584;
 x[447][3] = 3.6884;
 x[447][4] = 0.74912;
 y[447][0] = a;

 x[448][0] = 0.001;
 x[448][1] = -2.0759;
 x[448][2] = 10.8223;
 x[448][3] = 2.6439;
 x[448][4] = -4.837;
 y[448][0] = a;

 x[449][0] = 0.001;
 x[449][1] = 4.0715;
 x[449][2] = 7.6398;
 x[449][3] = -2.0824;
 x[449][4] = -1.1698;
 y[449][0] = a;

 x[450][0] = 0.001;
 x[450][1] = 0.76163;
 x[450][2] = 5.8209;
 x[450][3] = 1.1959;
 x[450][4] = -0.64613;
 y[450][0] = a;

 x[451][0] = 0.001;
 x[451][1] = -0.53966;
 x[451][2] = 7.3273;
 x[451][3] = 0.46583;
 x[451][4] = -1.4543;
 y[451][0] = a;

 x[452][0] = 0.001;
 x[452][1] = 2.6213;
 x[452][2] = 5.7919;
 x[452][3] = 0.065686;
 x[452][4] = -1.5759;
 y[452][0] = a;

 x[453][0] = 0.001;
 x[453][1] = 3.0242;
 x[453][2] = -3.3378;
 x[453][3] = 2.5865;
 x[453][4] = -0.54785;
 y[453][0] = a;

 x[454][0] = 0.001;
 x[454][1] = 5.8519;
 x[454][2] = 5.3905;
 x[454][3] = -2.4037;
 x[454][4] = -0.061652;
 y[454][0] = a;

 x[455][0] = 0.001;
 x[455][1] = 0.5706;
 x[455][2] = -0.0248;
 x[455][3] = 1.2421;
 x[455][4] = -0.5621;
 y[455][0] = a;

 x[456][0] = 0.001;
 x[456][1] = 3.9771;
 x[456][2] = 11.1513;
 x[456][3] = -3.9272;
 x[456][4] = -4.3444;
 y[456][0] = a;

 x[457][0] = 0.001;
 x[457][1] = 1.5478;
 x[457][2] = 9.1814;
 x[457][3] = -1.6326;
 x[457][4] = -1.7375;
 y[457][0] = a;

 x[458][0] = 0.001;
 x[458][1] = 0.74054;
 x[458][2] = 0.36625;
 x[458][3] = 2.1992;
 x[458][4] = 0.48403;
 y[458][0] = a;

 x[459][0] = 0.001;
 x[459][1] = 0.49571;
 x[459][2] = 10.2243;
 x[459][3] = -1.097;
 x[459][4] = -4.0159;
 y[459][0] = a;

 x[460][0] = 0.001;
 x[460][1] = 1.645;
 x[460][2] = 7.8612;
 x[460][3] = -0.87598;
 x[460][4] = -3.5569;
 y[460][0] = a;

 x[461][0] = 0.001;
 x[461][1] = 3.6077;
 x[461][2] = 6.8576;
 x[461][3] = -1.1622;
 x[461][4] = 0.28231;
 y[461][0] = a;

 x[462][0] = 0.001;
 x[462][1] = 3.2403;
 x[462][2] = -3.7082;
 x[462][3] = 5.2804;
 x[462][4] = 0.41291;
 y[462][0] = a;

 x[463][0] = 0.001;
 x[463][1] = 3.9166;
 x[463][2] = 10.2491;
 x[463][3] = -4.0926;
 x[463][4] = -4.4659;
 y[463][0] = a;

 x[464][0] = 0.001;
 x[464][1] = 3.9262;
 x[464][2] = 6.0299;
 x[464][3] = -2.0156;
 x[464][4] = -0.065531;
 y[464][0] = a;

 x[465][0] = 0.001;
 x[465][1] = 5.591;
 x[465][2] = 10.4643;
 x[465][3] = -4.3839;
 x[465][4] = -4.3379;
 y[465][0] = a;

 x[466][0] = 0.001;
 x[466][1] = 3.7522;
 x[466][2] = -3.6978;
 x[466][3] = 3.9943;
 x[466][4] = 1.3051;
 y[466][0] = a;

 x[467][0] = 0.001;
 x[467][1] = 1.3114;
 x[467][2] = 4.5462;
 x[467][3] = 2.2935;
 x[467][4] = 0.22541;
 y[467][0] = a;

 x[468][0] = 0.001;
 x[468][1] = 3.7022;
 x[468][2] = 6.9942;
 x[468][3] = -1.8511;
 x[468][4] = -0.12889;
 y[468][0] = a;

 x[469][0] = 0.001;
 x[469][1] = 4.364;
 x[469][2] = -3.1039;
 x[469][3] = 2.3757;
 x[469][4] = 0.78532;
 y[469][0] = a;

 x[470][0] = 0.001;
 x[470][1] = 3.5829;
 x[470][2] = 1.4423;
 x[470][3] = 1.0219;
 x[470][4] = 1.4008;
 y[470][0] = a;

 x[471][0] = 0.001;
 x[471][1] = 4.65;
 x[471][2] = -4.8297;
 x[471][3] = 3.4553;
 x[471][4] = -0.25174;
 y[471][0] = a;

 x[472][0] = 0.001;
 x[472][1] = 5.1731;
 x[472][2] = 3.9606;
 x[472][3] = -1.983;
 x[472][4] = 0.40774;
 y[472][0] = a;

 x[473][0] = 0.001;
 x[473][1] = 3.2692;
 x[473][2] = 3.4184;
 x[473][3] = 0.20706;
 x[473][4] = -0.066824;
 y[473][0] = a;

 x[474][0] = 0.001;
 x[474][1] = 2.4012;
 x[474][2] = 1.6223;
 x[474][3] = 3.0312;
 x[474][4] = 0.71679;
 y[474][0] = a;

 x[475][0] = 0.001;
 x[475][1] = 1.7257;
 x[475][2] = -4.4697;
 x[475][3] = 8.2219;
 x[475][4] = -1.8073;
 y[475][0] = a;

 x[476][0] = 0.001;
 x[476][1] = 4.7965;
 x[476][2] = 6.9859;
 x[476][3] = -1.9967;
 x[476][4] = -0.35001;
 y[476][0] = a;

 x[477][0] = 0.001;
 x[477][1] = 4.0962;
 x[477][2] = 10.1891;
 x[477][3] = -3.9323;
 x[477][4] = -4.1827;
 y[477][0] = a;

 x[478][0] = 0.001;
 x[478][1] = 2.5559;
 x[478][2] = 3.3605;
 x[478][3] = 2.0321;
 x[478][4] = 0.26809;
 y[478][0] = a;

 x[479][0] = 0.001;
 x[479][1] = 3.4916;
 x[479][2] = 8.5709;
 x[479][3] = -3.0326;
 x[479][4] = -0.59182;
 y[479][0] = a;

 x[480][0] = 0.001;
 x[480][1] = 0.5195;
 x[480][2] = -3.2633;
 x[480][3] = 3.0895;
 x[480][4] = -0.9849;
 y[480][0] = a;

 x[481][0] = 0.001;
 x[481][1] = 2.9856;
 x[481][2] = 7.2673;
 x[481][3] = -0.409;
 x[481][4] = -2.2431;
 y[481][0] = a;

 x[482][0] = 0.001;
 x[482][1] = 4.0932;
 x[482][2] = 5.4132;
 x[482][3] = -1.8219;
 x[482][4] = 0.23576;
 y[482][0] = a;

 x[483][0] = 0.001;
 x[483][1] = 1.7748;
 x[483][2] = -0.76978;
 x[483][3] = 5.5854;
 x[483][4] = 1.3039;
 y[483][0] = a;

 x[484][0] = 0.001;
 x[484][1] = 5.2012;
 x[484][2] = 0.32694;
 x[484][3] = 0.17965;
 x[484][4] = 1.1797;
 y[484][0] = a;

 x[485][0] = 0.001;
 x[485][1] = -0.45062;
 x[485][2] = -1.3678;
 x[485][3] = 7.0858;
 x[485][4] = -0.40303;
 y[485][0] = a;

 x[486][0] = 0.001;
 x[486][1] = 4.8451;
 x[486][2] = 8.1116;
 x[486][3] = -2.9512;
 x[486][4] = -1.4724;
 y[486][0] = a;

 x[487][0] = 0.001;
 x[487][1] = 0.74841;
 x[487][2] = 7.2756;
 x[487][3] = 1.1504;
 x[487][4] = -0.5388;
 y[487][0] = a;

 x[488][0] = 0.001;
 x[488][1] = 5.1213;
 x[488][2] = 8.5565;
 x[488][3] = -3.3917;
 x[488][4] = -1.5474;
 y[488][0] = a;

 x[489][0] = 0.001;
 x[489][1] = 3.6181;
 x[489][2] = -3.7454;
 x[489][3] = 2.8273;
 x[489][4] = -0.71208;
 y[489][0] = a;

 x[490][0] = 0.001;
 x[490][1] = 0.040498;
 x[490][2] = 8.5234;
 x[490][3] = 1.4461;
 x[490][4] = -3.9306;
 y[490][0] = a;

 x[491][0] = 0.001;
 x[491][1] = -2.6479;
 x[491][2] = 10.1374;
 x[491][3] = -1.331;
 x[491][4] = -5.4707;
 y[491][0] = a;

 x[492][0] = 0.001;
 x[492][1] = 0.37984;
 x[492][2] = 0.70975;
 x[492][3] = 0.75716;
 x[492][4] = -0.44441;
 y[492][0] = a;

 x[493][0] = 0.001;
 x[493][1] = -0.95923;
 x[493][2] = 0.091039;
 x[493][3] = 6.2204;
 x[493][4] = -1.4828;
 y[493][0] = a;

 x[494][0] = 0.001;
 x[494][1] = 2.8672;
 x[494][2] = 10.0008;
 x[494][3] = -3.2049;
 x[494][4] = -3.1095;
 y[494][0] = a;

 x[495][0] = 0.001;
 x[495][1] = 1.0182;
 x[495][2] = 9.109;
 x[495][3] = -0.62064;
 x[495][4] = -1.7129;
 y[495][0] = a;

 x[496][0] = 0.001;
 x[496][1] = -2.7143;
 x[496][2] = 11.4535;
 x[496][3] = 2.1092;
 x[496][4] = -3.9629;
 y[496][0] = a;

 x[497][0] = 0.001;
 x[497][1] = 3.8244;
 x[497][2] = -3.1081;
 x[497][3] = 2.4537;
 x[497][4] = 0.52024;
 y[497][0] = a;

 x[498][0] = 0.001;
 x[498][1] = 2.7961;
 x[498][2] = 2.121;
 x[498][3] = 1.8385;
 x[498][4] = 0.38317;
 y[498][0] = a;

 x[499][0] = 0.001;
 x[499][1] = 3.5358;
 x[499][2] = 6.7086;
 x[499][3] = -0.81857;
 x[499][4] = 0.47886;
 y[499][0] = a;

 x[500][0] = 0.001;
 x[500][1] = -0.7056;
 x[500][2] = 8.7241;
 x[500][3] = 2.2215;
 x[500][4] = -4.5965;
 y[500][0] = a;

 x[501][0] = 0.001;
 x[501][1] = 4.1542;
 x[501][2] = 7.2756;
 x[501][3] = -2.4766;
 x[501][4] = -1.2099;
 y[501][0] = a;

 x[502][0] = 0.001;
 x[502][1] = 0.92703;
 x[502][2] = 9.4318;
 x[502][3] = -0.66263;
 x[502][4] = -1.6728;
 y[502][0] = a;

 x[503][0] = 0.001;
 x[503][1] = 1.8216;
 x[503][2] = -6.4748;
 x[503][3] = 8.0514;
 x[503][4] = -0.41855;
 y[503][0] = a;

 x[504][0] = 0.001;
 x[504][1] = -2.4473;
 x[504][2] = 12.6247;
 x[504][3] = 0.73573;
 x[504][4] = -7.6612;
 y[504][0] = a;

 x[505][0] = 0.001;
 x[505][1] = 3.5862;
 x[505][2] = -3.0957;
 x[505][3] = 2.8093;
 x[505][4] = 0.24481;
 y[505][0] = a;

 x[506][0] = 0.001;
 x[506][1] = 0.66191;
 x[506][2] = 9.6594;
 x[506][3] = -0.28819;
 x[506][4] = -1.6638;
 y[506][0] = a;

 x[507][0] = 0.001;
 x[507][1] = 4.7926;
 x[507][2] = 1.7071;
 x[507][3] = -0.051701;
 x[507][4] = 1.4926;
 y[507][0] = a;

 x[508][0] = 0.001;
 x[508][1] = 4.9852;
 x[508][2] = 8.3516;
 x[508][3] = -2.5425;
 x[508][4] = -1.2823;
 y[508][0] = a;

 x[509][0] = 0.001;
 x[509][1] = 0.75736;
 x[509][2] = 3.0294;
 x[509][3] = 2.9164;
 x[509][4] = -0.068117;
 y[509][0] = a;

 x[510][0] = 0.001;
 x[510][1] = 4.6499;
 x[510][2] = 7.6336;
 x[510][3] = -1.9427;
 x[510][4] = -0.37458;
 y[510][0] = a;

 x[511][0] = 0.001;
 x[511][1] = -0.023579;
 x[511][2] = 7.1742;
 x[511][3] = 0.78457;
 x[511][4] = -0.75734;
 y[511][0] = a;

 x[512][0] = 0.001;
 x[512][1] = 0.85574;
 x[512][2] = 0.008268;
 x[512][3] = 6.6042;
 x[512][4] = -0.53104;
 y[512][0] = a;

 x[513][0] = 0.001;
 x[513][1] = 0.88298;
 x[513][2] = 0.66009;
 x[513][3] = 6.0096;
 x[513][4] = -0.43277;
 y[513][0] = a;

 x[514][0] = 0.001;
 x[514][1] = 4.0422;
 x[514][2] = -4.391;
 x[514][3] = 4.7466;
 x[514][4] = 1.137;
 y[514][0] = a;

 x[515][0] = 0.001;
 x[515][1] = 2.2546;
 x[515][2] = 8.0992;
 x[515][3] = -0.24877;
 x[515][4] = -3.2698;
 y[515][0] = a;

 x[516][0] = 0.001;
 x[516][1] = 0.38478;
 x[516][2] = 6.5989;
 x[516][3] = -0.3336;
 x[516][4] = -0.56466;
 y[516][0] = a;

 x[517][0] = 0.001;
 x[517][1] = 3.1541;
 x[517][2] = -5.1711;
 x[517][3] = 6.5991;
 x[517][4] = 0.57455;
 y[517][0] = a;

 x[518][0] = 0.001;
 x[518][1] = 2.3969;
 x[518][2] = 0.23589;
 x[518][3] = 4.8477;
 x[518][4] = 1.437;
 y[518][0] = a;

 x[519][0] = 0.001;
 x[519][1] = 4.7114;
 x[519][2] = 2.0755;
 x[519][3] = -0.2702;
 x[519][4] = 1.2379;
 y[519][0] = a;

 x[520][0] = 0.001;
 x[520][1] = 4.0127;
 x[520][2] = 10.1477;
 x[520][3] = -3.9366;
 x[520][4] = -4.0728;
 y[520][0] = a;

 x[521][0] = 0.001;
 x[521][1] = 2.6606;
 x[521][2] = 3.1681;
 x[521][3] = 1.9619;
 x[521][4] = 0.18662;
 y[521][0] = a;

 x[522][0] = 0.001;
 x[522][1] = 3.931;
 x[522][2] = 1.8541;
 x[522][3] = -0.023425;
 x[522][4] = 1.2314;
 y[522][0] = a;

 x[523][0] = 0.001;
 x[523][1] = 0.01727;
 x[523][2] = 8.693;
 x[523][3] = 1.3989;
 x[523][4] = -3.9668;
 y[523][0] = a;

 x[524][0] = 0.001;
 x[524][1] = 3.2414;
 x[524][2] = 0.40971;
 x[524][3] = 1.4015;
 x[524][4] = 1.1952;
 y[524][0] = a;

 x[525][0] = 0.001;
 x[525][1] = 2.2504;
 x[525][2] = 3.5757;
 x[525][3] = 0.35273;
 x[525][4] = 0.2836;
 y[525][0] = a;

 x[526][0] = 0.001;
 x[526][1] = -1.3971;
 x[526][2] = 3.3191;
 x[526][3] = -1.3927;
 x[526][4] = -1.9948;
 y[526][0] = b;

 x[527][0] = 0.001;
 x[527][1] = 0.39012;
 x[527][2] = -0.14279;
 x[527][3] = -0.031994;
 x[527][4] = 0.35084;
 y[527][0] = b;

 x[528][0] = 0.001;
 x[528][1] = -1.6677;
 x[528][2] = -7.1535;
 x[528][3] = 7.8929;
 x[528][4] = 0.96765;
 y[528][0] = b;

 x[529][0] = 0.001;
 x[529][1] = -3.8483;
 x[529][2] = -12.8047;
 x[529][3] = 15.6824;
 x[529][4] = -1.281;
 y[529][0] = b;

 x[530][0] = 0.001;
 x[530][1] = -3.5681;
 x[530][2] = -8.213;
 x[530][3] = 10.083;
 x[530][4] = 0.96765;
 y[530][0] = b;

 x[531][0] = 0.001;
 x[531][1] = -2.2804;
 x[531][2] = -0.30626;
 x[531][3] = 1.3347;
 x[531][4] = 1.3763;
 y[531][0] = b;

 x[532][0] = 0.001;
 x[532][1] = -1.7582;
 x[532][2] = 2.7397;
 x[532][3] = -2.5323;
 x[532][4] = -2.234;
 y[532][0] = b;

 x[533][0] = 0.001;
 x[533][1] = -0.89409;
 x[533][2] = 3.1991;
 x[533][3] = -1.8219;
 x[533][4] = -2.9452;
 y[533][0] = b;

 x[534][0] = 0.001;
 x[534][1] = 0.3434;
 x[534][2] = 0.12415;
 x[534][3] = -0.28733;
 x[534][4] = 0.14654;
 y[534][0] = b;

 x[535][0] = 0.001;
 x[535][1] = -0.9854;
 x[535][2] = -6.661;
 x[535][3] = 5.8245;
 x[535][4] = 0.5461;
 y[535][0] = b;

 x[536][0] = 0.001;
 x[536][1] = -2.4115;
 x[536][2] = -9.1359;
 x[536][3] = 9.3444;
 x[536][4] = -0.65259;
 y[536][0] = b;

 x[537][0] = 0.001;
 x[537][1] = -1.5252;
 x[537][2] = -6.2534;
 x[537][3] = 5.3524;
 x[537][4] = 0.59912;
 y[537][0] = b;

 x[538][0] = 0.001;
 x[538][1] = -0.61442;
 x[538][2] = -0.091058;
 x[538][3] = -0.31818;
 x[538][4] = 0.50214;
 y[538][0] = b;

 x[539][0] = 0.001;
 x[539][1] = -0.36506;
 x[539][2] = 2.8928;
 x[539][3] = -3.6461;
 x[539][4] = -3.0603;
 y[539][0] = b;

 x[540][0] = 0.001;
 x[540][1] = -5.9034;
 x[540][2] = 6.5679;
 x[540][3] = 0.67661;
 x[540][4] = -6.6797;
 y[540][0] = b;

 x[541][0] = 0.001;
 x[541][1] = -1.8215;
 x[541][2] = 2.7521;
 x[541][3] = -0.72261;
 x[541][4] = -2.353;
 y[541][0] = b;

 x[542][0] = 0.001;
 x[542][1] = -0.77461;
 x[542][2] = -1.8768;
 x[542][3] = 2.4023;
 x[542][4] = 1.1319;
 y[542][0] = b;

 x[543][0] = 0.001;
 x[543][1] = -1.8187;
 x[543][2] = -9.0366;
 x[543][3] = 9.0162;
 x[543][4] = -0.12243;
 y[543][0] = b;

 x[544][0] = 0.001;
 x[544][1] = -3.5801;
 x[544][2] = -12.9309;
 x[544][3] = 13.1779;
 x[544][4] = -2.5677;
 y[544][0] = b;

 x[545][0] = 0.001;
 x[545][1] = -1.8219;
 x[545][2] = -6.8824;
 x[545][3] = 5.4681;
 x[545][4] = 0.057313;
 y[545][0] = b;

 x[546][0] = 0.001;
 x[546][1] = -0.3481;
 x[546][2] = -0.38696;
 x[546][3] = -0.47841;
 x[546][4] = 0.62627;
 y[546][0] = b;

 x[547][0] = 0.001;
 x[547][1] = 0.47368;
 x[547][2] = 3.3605;
 x[547][3] = -4.5064;
 x[547][4] = -4.0431;
 y[547][0] = b;

 x[548][0] = 0.001;
 x[548][1] = -3.4083;
 x[548][2] = 4.8587;
 x[548][3] = -0.76888;
 x[548][4] = -4.8668;
 y[548][0] = b;

 x[549][0] = 0.001;
 x[549][1] = -1.6662;
 x[549][2] = -0.30005;
 x[549][3] = 1.4238;
 x[549][4] = 0.024986;
 y[549][0] = b;

 x[550][0] = 0.001;
 x[550][1] = -2.0962;
 x[550][2] = -7.1059;
 x[550][3] = 6.6188;
 x[550][4] = -0.33708;
 y[550][0] = b;

 x[551][0] = 0.001;
 x[551][1] = -2.6685;
 x[551][2] = -10.4519;
 x[551][3] = 9.1139;
 x[551][4] = -1.7323;
 y[551][0] = b;

 x[552][0] = 0.001;
 x[552][1] = -0.47465;
 x[552][2] = -4.3496;
 x[552][3] = 1.9901;
 x[552][4] = 0.7517;
 y[552][0] = b;

 x[553][0] = 0.001;
 x[553][1] = 1.0552;
 x[553][2] = 1.1857;
 x[553][3] = -2.6411;
 x[553][4] = 0.11033;
 y[553][0] = b;

 x[554][0] = 0.001;
 x[554][1] = 1.1644;
 x[554][2] = 3.8095;
 x[554][3] = -4.9408;
 x[554][4] = -4.0909;
 y[554][0] = b;

 x[555][0] = 0.001;
 x[555][1] = -4.4779;
 x[555][2] = 7.3708;
 x[555][3] = -0.31218;
 x[555][4] = -6.7754;
 y[555][0] = b;

 x[556][0] = 0.001;
 x[556][1] = -2.7338;
 x[556][2] = 0.45523;
 x[556][3] = 2.4391;
 x[556][4] = 0.21766;
 y[556][0] = b;

 x[557][0] = 0.001;
 x[557][1] = -2.286;
 x[557][2] = -5.4484;
 x[557][3] = 5.8039;
 x[557][4] = 0.88231;
 y[557][0] = b;

 x[558][0] = 0.001;
 x[558][1] = -1.6244;
 x[558][2] = -6.3444;
 x[558][3] = 4.6575;
 x[558][4] = 0.16981;
 y[558][0] = b;

 x[559][0] = 0.001;
 x[559][1] = 0.50813;
 x[559][2] = 0.47799;
 x[559][3] = -1.9804;
 x[559][4] = 0.57714;
 y[559][0] = b;

 x[560][0] = 0.001;
 x[560][1] = 1.6408;
 x[560][2] = 4.2503;
 x[560][3] = -4.9023;
 x[560][4] = -2.6621;
 y[560][0] = b;

 x[561][0] = 0.001;
 x[561][1] = 0.81583;
 x[561][2] = 4.84;
 x[561][3] = -5.2613;
 x[561][4] = -6.0823;
 y[561][0] = b;

 x[562][0] = 0.001;
 x[562][1] = -5.4901;
 x[562][2] = 9.1048;
 x[562][3] = -0.38758;
 x[562][4] = -5.9763;
 y[562][0] = b;

 x[563][0] = 0.001;
 x[563][1] = -3.2238;
 x[563][2] = 2.7935;
 x[563][3] = 0.32274;
 x[563][4] = -0.86078;
 y[563][0] = b;

 x[564][0] = 0.001;
 x[564][1] = -2.0631;
 x[564][2] = -1.5147;
 x[564][3] = 1.219;
 x[564][4] = 0.44524;
 y[564][0] = b;

 x[565][0] = 0.001;
 x[565][1] = -0.91318;
 x[565][2] = -2.0113;
 x[565][3] = -0.19565;
 x[565][4] = 0.066365;
 y[565][0] = b;

 x[566][0] = 0.001;
 x[566][1] = 0.6005;
 x[566][2] = 1.9327;
 x[566][3] = -3.2888;
 x[566][4] = -0.32415;
 y[566][0] = b;

 x[567][0] = 0.001;
 x[567][1] = 0.91315;
 x[567][2] = 3.3377;
 x[567][3] = -4.0557;
 x[567][4] = -1.6741;
 y[567][0] = b;

 x[568][0] = 0.001;
 x[568][1] = -0.28015;
 x[568][2] = 3.0729;
 x[568][3] = -3.3857;
 x[568][4] = -2.9155;
 y[568][0] = b;

 x[569][0] = 0.001;
 x[569][1] = -3.6085;
 x[569][2] = 3.3253;
 x[569][3] = -0.51954;
 x[569][4] = -3.5737;
 y[569][0] = b;

 x[570][0] = 0.001;
 x[570][1] = -6.2003;
 x[570][2] = 8.6806;
 x[570][3] = 0.009134;
 x[570][4] = -3.703;
 y[570][0] = b;

 x[571][0] = 0.001;
 x[571][1] = -4.2932;
 x[571][2] = 3.3419;
 x[571][3] = 0.77258;
 x[571][4] = -0.99785;
 y[571][0] = b;

 x[572][0] = 0.001;
 x[572][1] = -3.0265;
 x[572][2] = -0.062088;
 x[572][3] = 0.68604;
 x[572][4] = -0.055186;
 y[572][0] = b;

 x[573][0] = 0.001;
 x[573][1] = -1.7015;
 x[573][2] = -0.010356;
 x[573][3] = -0.99337;
 x[573][4] = -0.53104;
 y[573][0] = b;

 x[574][0] = 0.001;
 x[574][1] = -0.64326;
 x[574][2] = 2.4748;
 x[574][3] = -2.9452;
 x[574][4] = -1.0276;
 y[574][0] = b;

 x[575][0] = 0.001;
 x[575][1] = -0.86339;
 x[575][2] = 1.9348;
 x[575][3] = -2.3729;
 x[575][4] = -1.0897;
 y[575][0] = b;

 x[576][0] = 0.001;
 x[576][1] = -2.0659;
 x[576][2] = 1.0512;
 x[576][3] = -0.46298;
 x[576][4] = -1.0974;
 y[576][0] = b;

 x[577][0] = 0.001;
 x[577][1] = -2.1333;
 x[577][2] = 1.5685;
 x[577][3] = -0.084261;
 x[577][4] = -1.7453;
 y[577][0] = b;

 x[578][0] = 0.001;
 x[578][1] = -1.2568;
 x[578][2] = -1.4733;
 x[578][3] = 2.8718;
 x[578][4] = 0.44653;
 y[578][0] = b;

 x[579][0] = 0.001;
 x[579][1] = -3.1128;
 x[579][2] = -6.841;
 x[579][3] = 10.7402;
 x[579][4] = -1.0172;
 y[579][0] = b;

 x[580][0] = 0.001;
 x[580][1] = -4.8554;
 x[580][2] = -5.9037;
 x[580][3] = 10.9818;
 x[580][4] = -0.82199;
 y[580][0] = b;

 x[581][0] = 0.001;
 x[581][1] = -2.588;
 x[581][2] = 3.8654;
 x[581][3] = -0.3336;
 x[581][4] = -1.2797;
 y[581][0] = b;

 x[582][0] = 0.001;
 x[582][1] = 0.24394;
 x[582][2] = 1.4733;
 x[582][3] = -1.4192;
 x[582][4] = -0.58535;
 y[582][0] = b;

 x[583][0] = 0.001;
 x[583][1] = -1.5322;
 x[583][2] = -5.0966;
 x[583][3] = 6.6779;
 x[583][4] = 0.17498;
 y[583][0] = b;

 x[584][0] = 0.001;
 x[584][1] = -4.0025;
 x[584][2] = -13.4979;
 x[584][3] = 17.6772;
 x[584][4] = -3.3202;
 y[584][0] = b;

 x[585][0] = 0.001;
 x[585][1] = -4.0173;
 x[585][2] = -8.3123;
 x[585][3] = 12.4547;
 x[585][4] = -1.4375;
 y[585][0] = b;

 x[586][0] = 0.001;
 x[586][1] = -3.0731;
 x[586][2] = -0.53181;
 x[586][3] = 2.3877;
 x[586][4] = 0.77627;
 y[586][0] = b;

 x[587][0] = 0.001;
 x[587][1] = -1.979;
 x[587][2] = 3.2301;
 x[587][3] = -1.3575;
 x[587][4] = -2.5819;
 y[587][0] = b;

 x[588][0] = 0.001;
 x[588][1] = -0.4294;
 x[588][2] = -0.14693;
 x[588][3] = 0.044265;
 x[588][4] = -0.15605;
 y[588][0] = b;

 x[589][0] = 0.001;
 x[589][1] = -2.234;
 x[589][2] = -7.0314;
 x[589][3] = 7.4936;
 x[589][4] = 0.61334;
 y[589][0] = b;

 x[590][0] = 0.001;
 x[590][1] = -4.211;
 x[590][2] = -12.4736;
 x[590][3] = 14.9704;
 x[590][4] = -1.3884;
 y[590][0] = b;

 x[591][0] = 0.001;
 x[591][1] = -3.8073;
 x[591][2] = -8.0971;
 x[591][3] = 10.1772;
 x[591][4] = 0.65084;
 y[591][0] = b;

 x[592][0] = 0.001;
 x[592][1] = -2.5912;
 x[592][2] = -0.10554;
 x[592][3] = 1.2798;
 x[592][4] = 1.0414;
 y[592][0] = b;

 x[593][0] = 0.001;
 x[593][1] = -2.2482;
 x[593][2] = 3.0915;
 x[593][3] = -2.3969;
 x[593][4] = -2.6711;
 y[593][0] = b;

 x[594][0] = 0.001;
 x[594][1] = -1.4427;
 x[594][2] = 3.2922;
 x[594][3] = -1.9702;
 x[594][4] = -3.4392;
 y[594][0] = b;

 x[595][0] = 0.001;
 x[595][1] = -0.39416;
 x[595][2] = -0.020702;
 x[595][3] = -0.066267;
 x[595][4] = -0.44699;
 y[595][0] = b;

 x[596][0] = 0.001;
 x[596][1] = -1.522;
 x[596][2] = -6.6383;
 x[596][3] = 5.7491;
 x[596][4] = -0.10691;
 y[596][0] = b;

 x[597][0] = 0.001;
 x[597][1] = -2.8267;
 x[597][2] = -9.0407;
 x[597][3] = 9.0694;
 x[597][4] = -0.98233;
 y[597][0] = b;

 x[598][0] = 0.001;
 x[598][1] = -1.7263;
 x[598][2] = -6.0237;
 x[598][3] = 5.2419;
 x[598][4] = 0.29524;
 y[598][0] = b;

 x[599][0] = 0.001;
 x[599][1] = -0.94255;
 x[599][2] = 0.039307;
 x[599][3] = -0.24192;
 x[599][4] = 0.31593;
 y[599][0] = b;

 x[600][0] = 0.001;
 x[600][1] = -0.89569;
 x[600][2] = 3.0025;
 x[600][3] = -3.6067;
 x[600][4] = -3.4457;
 y[600][0] = b;

 x[601][0] = 0.001;
 x[601][1] = -6.2815;
 x[601][2] = 6.6651;
 x[601][3] = 0.52581;
 x[601][4] = -7.0107;
 y[601][0] = b;

 x[602][0] = 0.001;
 x[602][1] = -2.3211;
 x[602][2] = 3.166;
 x[602][3] = -1.0002;
 x[602][4] = -2.7151;
 y[602][0] = b;

 x[603][0] = 0.001;
 x[603][1] = -1.3414;
 x[603][2] = -2.0776;
 x[603][3] = 2.8093;
 x[603][4] = 0.60688;
 y[603][0] = b;

 x[604][0] = 0.001;
 x[604][1] = -2.258;
 x[604][2] = -9.3263;
 x[604][3] = 9.3727;
 x[604][4] = -0.85949;
 y[604][0] = b;

 x[605][0] = 0.001;
 x[605][1] = -3.8858;
 x[605][2] = -12.8461;
 x[605][3] = 12.7957;
 x[605][4] = -3.1353;
 y[605][0] = b;

 x[606][0] = 0.001;
 x[606][1] = -1.8969;
 x[606][2] = -6.7893;
 x[606][3] = 5.2761;
 x[606][4] = -0.32544;
 y[606][0] = b;

 x[607][0] = 0.001;
 x[607][1] = -0.52645;
 x[607][2] = -0.24832;
 x[607][3] = -0.45613;
 x[607][4] = 0.41938;
 y[607][0] = b;

 x[608][0] = 0.001;
 x[608][1] = 0.009661;
 x[608][2] = 3.5612;
 x[608][3] = -4.407;
 x[608][4] = -4.4103;
 y[608][0] = b;

 x[609][0] = 0.001;
 x[609][1] = -3.8826;
 x[609][2] = 4.898;
 x[609][3] = -0.92311;
 x[609][4] = -5.0801;
 y[609][0] = b;

 x[610][0] = 0.001;
 x[610][1] = -2.1405;
 x[610][2] = -0.16762;
 x[610][3] = 1.321;
 x[610][4] = -0.20906;
 y[610][0] = b;

 x[611][0] = 0.001;
 x[611][1] = -2.4824;
 x[611][2] = -7.3046;
 x[611][3] = 6.839;
 x[611][4] = -0.59053;
 y[611][0] = b;

 x[612][0] = 0.001;
 x[612][1] = -2.9098;
 x[612][2] = -10.0712;
 x[612][3] = 8.4156;
 x[612][4] = -1.9948;
 y[612][0] = b;

 x[613][0] = 0.001;
 x[613][1] = -0.60975;
 x[613][2] = -4.002;
 x[613][3] = 1.8471;
 x[613][4] = 0.6017;
 y[613][0] = b;

 x[614][0] = 0.001;
 x[614][1] = 0.83625;
 x[614][2] = 1.1071;
 x[614][3] = -2.4706;
 x[614][4] = -0.062945;
 y[614][0] = b;

 x[615][0] = 0.001;
 x[615][1] = 0.60731;
 x[615][2] = 3.9544;
 x[615][3] = -4.772;
 x[615][4] = -4.4853;
 y[615][0] = b;

 x[616][0] = 0.001;
 x[616][1] = -4.8861;
 x[616][2] = 7.0542;
 x[616][3] = -0.17252;
 x[616][4] = -6.959;
 y[616][0] = b;

 x[617][0] = 0.001;
 x[617][1] = -3.1366;
 x[617][2] = 0.42212;
 x[617][3] = 2.6225;
 x[617][4] = -0.064238;
 y[617][0] = b;

 x[618][0] = 0.001;
 x[618][1] = -2.5754;
 x[618][2] = -5.6574;
 x[618][3] = 6.103;
 x[618][4] = 0.65214;
 y[618][0] = b;

 x[619][0] = 0.001;
 x[619][1] = -1.8782;
 x[619][2] = -6.5865;
 x[619][3] = 4.8486;
 x[619][4] = -0.021566;
 y[619][0] = b;

 x[620][0] = 0.001;
 x[620][1] = 0.24261;
 x[620][2] = 0.57318;
 x[620][3] = -1.9402;
 x[620][4] = 0.44007;
 y[620][0] = b;

 x[621][0] = 0.001;
 x[621][1] = 1.296;
 x[621][2] = 4.2855;
 x[621][3] = -4.8457;
 x[621][4] = -2.9013;
 y[621][0] = b;

 x[622][0] = 0.001;
 x[622][1] = 0.25943;
 x[622][2] = 5.0097;
 x[622][3] = -5.0394;
 x[622][4] = -6.3862;
 y[622][0] = b;

 x[623][0] = 0.001;
 x[623][1] = -5.873;
 x[623][2] = 9.1752;
 x[623][3] = -0.27448;
 x[623][4] = -6.0422;
 y[623][0] = b;

 x[624][0] = 0.001;
 x[624][1] = -3.4605;
 x[624][2] = 2.6901;
 x[624][3] = 0.16165;
 x[624][4] = -1.0224;
 y[624][0] = b;

 x[625][0] = 0.001;
 x[625][1] = -2.3797;
 x[625][2] = -1.4402;
 x[625][3] = 1.1273;
 x[625][4] = 0.16076;
 y[625][0] = b;

 x[626][0] = 0.001;
 x[626][1] = -1.2424;
 x[626][2] = -1.7175;
 x[626][3] = -0.52553;
 x[626][4] = -0.21036;
 y[626][0] = b;

 x[627][0] = 0.001;
 x[627][1] = 0.20216;
 x[627][2] = 1.9182;
 x[627][3] = -3.2828;
 x[627][4] = -0.61768;
 y[627][0] = b;

 x[628][0] = 0.001;
 x[628][1] = 0.59823;
 x[628][2] = 3.5012;
 x[628][3] = -3.9795;
 x[628][4] = -1.7841;
 y[628][0] = b;

 x[629][0] = 0.001;
 x[629][1] = -0.77995;
 x[629][2] = 3.2322;
 x[629][3] = -3.282;
 x[629][4] = -3.1004;
 y[629][0] = b;

 x[630][0] = 0.001;
 x[630][1] = -4.1409;
 x[630][2] = 3.4619;
 x[630][3] = -0.47841;
 x[630][4] = -3.8879;
 y[630][0] = b;

 x[631][0] = 0.001;
 x[631][1] = -6.5084;
 x[631][2] = 8.7696;
 x[631][3] = 0.23191;
 x[631][4] = -3.937;
 y[631][0] = b;

 x[632][0] = 0.001;
 x[632][1] = -4.4996;
 x[632][2] = 3.4288;
 x[632][3] = 0.56265;
 x[632][4] = -1.1672;
 y[632][0] = b;

 x[633][0] = 0.001;
 x[633][1] = -3.3125;
 x[633][2] = 0.10139;
 x[633][3] = 0.55323;
 x[633][4] = -0.2957;
 y[633][0] = b;

 x[634][0] = 0.001;
 x[634][1] = -1.9423;
 x[634][2] = 0.3766;
 x[634][3] = -1.2898;
 x[634][4] = -0.82458;
 y[634][0] = b;

 x[635][0] = 0.001;
 x[635][1] = -0.75793;
 x[635][2] = 2.5349;
 x[635][3] = -3.0464;
 x[635][4] = -1.2629;
 y[635][0] = b;

 x[636][0] = 0.001;
 x[636][1] = -0.95403;
 x[636][2] = 1.9824;
 x[636][3] = -2.3163;
 x[636][4] = -1.1957;
 y[636][0] = b;

 x[637][0] = 0.001;
 x[637][1] = -2.2173;
 x[637][2] = 1.4671;
 x[637][3] = -0.72689;
 x[637][4] = -1.1724;
 y[637][0] = b;

 x[638][0] = 0.001;
 x[638][1] = -2.799;
 x[638][2] = 1.9679;
 x[638][3] = -0.42357;
 x[638][4] = -2.1125;
 y[638][0] = b;

 x[639][0] = 0.001;
 x[639][1] = -1.8629;
 x[639][2] = -0.84841;
 x[639][3] = 2.5377;
 x[639][4] = 0.097399;
 y[639][0] = b;

 x[640][0] = 0.001;
 x[640][1] = -3.5916;
 x[640][2] = -6.2285;
 x[640][3] = 10.2389;
 x[640][4] = -1.1543;
 y[640][0] = b;

 x[641][0] = 0.001;
 x[641][1] = -5.1216;
 x[641][2] = -5.3118;
 x[641][3] = 10.3846;
 x[641][4] = -1.0612;
 y[641][0] = b;

 x[642][0] = 0.001;
 x[642][1] = -3.2854;
 x[642][2] = 4.0372;
 x[642][3] = -0.45356;
 x[642][4] = -1.8228;
 y[642][0] = b;

 x[643][0] = 0.001;
 x[643][1] = -0.56877;
 x[643][2] = 1.4174;
 x[643][3] = -1.4252;
 x[643][4] = -1.1246;
 y[643][0] = b;

 x[644][0] = 0.001;
 x[644][1] = -2.3518;
 x[644][2] = -4.8359;
 x[644][3] = 6.6479;
 x[644][4] = -0.060358;
 y[644][0] = b;

 x[645][0] = 0.001;
 x[645][1] = -4.4861;
 x[645][2] = -13.2889;
 x[645][3] = 17.3087;
 x[645][4] = -3.2194;
 y[645][0] = b;

 x[646][0] = 0.001;
 x[646][1] = -4.3876;
 x[646][2] = -7.7267;
 x[646][3] = 11.9655;
 x[646][4] = -1.4543;
 y[646][0] = b;

 x[647][0] = 0.001;
 x[647][1] = -3.3604;
 x[647][2] = -0.32696;
 x[647][3] = 2.1324;
 x[647][4] = 0.6017;
 y[647][0] = b;

 x[648][0] = 0.001;
 x[648][1] = -1.0112;
 x[648][2] = 2.9984;
 x[648][3] = -1.1664;
 x[648][4] = -1.6185;
 y[648][0] = b;

 x[649][0] = 0.001;
 x[649][1] = 0.030219;
 x[649][2] = -1.0512;
 x[649][3] = 1.4024;
 x[649][4] = 0.77369;
 y[649][0] = b;

 x[650][0] = 0.001;
 x[650][1] = -1.6514;
 x[650][2] = -8.4985;
 x[650][3] = 9.1122;
 x[650][4] = 1.2379;
 y[650][0] = b;

 x[651][0] = 0.001;
 x[651][1] = -3.2692;
 x[651][2] = -12.7406;
 x[651][3] = 15.5573;
 x[651][4] = -0.14182;
 y[651][0] = b;

 x[652][0] = 0.001;
 x[652][1] = -2.5701;
 x[652][2] = -6.8452;
 x[652][3] = 8.9999;
 x[652][4] = 2.1353;
 y[652][0] = b;

 x[653][0] = 0.001;
 x[653][1] = -1.3066;
 x[653][2] = 0.25244;
 x[653][3] = 0.7623;
 x[653][4] = 1.7758;
 y[653][0] = b;

 x[654][0] = 0.001;
 x[654][1] = -1.6637;
 x[654][2] = 3.2881;
 x[654][3] = -2.2701;
 x[654][4] = -2.2224;
 y[654][0] = b;

 x[655][0] = 0.001;
 x[655][1] = -0.55008;
 x[655][2] = 2.8659;
 x[655][3] = -1.6488;
 x[655][4] = -2.4319;
 y[655][0] = b;

 x[656][0] = 0.001;
 x[656][1] = 0.21431;
 x[656][2] = -0.69529;
 x[656][3] = 0.87711;
 x[656][4] = 0.29653;
 y[656][0] = b;

 x[657][0] = 0.001;
 x[657][1] = -0.77288;
 x[657][2] = -7.4473;
 x[657][3] = 6.492;
 x[657][4] = 0.36119;
 y[657][0] = b;

 x[658][0] = 0.001;
 x[658][1] = -1.8391;
 x[658][2] = -9.0883;
 x[658][3] = 9.2416;
 x[658][4] = -0.10432;
 y[658][0] = b;

 x[659][0] = 0.001;
 x[659][1] = -0.63298;
 x[659][2] = -5.1277;
 x[659][3] = 4.5624;
 x[659][4] = 1.4797;
 y[659][0] = b;

 x[660][0] = 0.001;
 x[660][1] = 0.004055;
 x[660][2] = 0.62905;
 x[660][3] = -0.64121;
 x[660][4] = 0.75817;
 y[660][0] = b;

 x[661][0] = 0.001;
 x[661][1] = -0.28696;
 x[661][2] = 3.1784;
 x[661][3] = -3.5767;
 x[661][4] = -3.1896;
 y[661][0] = b;

 x[662][0] = 0.001;
 x[662][1] = -5.2406;
 x[662][2] = 6.6258;
 x[662][3] = -0.19908;
 x[662][4] = -6.8607;
 y[662][0] = b;

 x[663][0] = 0.001;
 x[663][1] = -1.4446;
 x[663][2] = 2.1438;
 x[663][3] = -0.47241;
 x[663][4] = -1.6677;
 y[663][0] = b;

 x[664][0] = 0.001;
 x[664][1] = -0.65767;
 x[664][2] = -2.8018;
 x[664][3] = 3.7115;
 x[664][4] = 0.99739;
 y[664][0] = b;

 x[665][0] = 0.001;
 x[665][1] = -1.5449;
 x[665][2] = -10.1498;
 x[665][3] = 9.6152;
 x[665][4] = -1.2332;
 y[665][0] = b;

 x[666][0] = 0.001;
 x[666][1] = -2.8957;
 x[666][2] = -12.0205;
 x[666][3] = 11.9149;
 x[666][4] = -2.7552;
 y[666][0] = b;

 x[667][0] = 0.001;
 x[667][1] = -0.81479;
 x[667][2] = -5.7381;
 x[667][3] = 4.3919;
 x[667][4] = 0.3211;
 y[667][0] = b;

 x[668][0] = 0.001;
 x[668][1] = 0.50225;
 x[668][2] = 0.65388;
 x[668][3] = -1.1793;
 x[668][4] = 0.39998;
 y[668][0] = b;

 x[669][0] = 0.001;
 x[669][1] = 0.74521;
 x[669][2] = 3.6357;
 x[669][3] = -4.4044;
 x[669][4] = -4.1414;
 y[669][0] = b;

 x[670][0] = 0.001;
 x[670][1] = -2.9146;
 x[670][2] = 4.0537;
 x[670][3] = -0.45699;
 x[670][4] = -4.0327;
 y[670][0] = b;

 x[671][0] = 0.001;
 x[671][1] = -1.3907;
 x[671][2] = -1.3781;
 x[671][3] = 2.3055;
 x[671][4] = -0.021566;
 y[671][0] = b;

 x[672][0] = 0.001;
 x[672][1] = -1.786;
 x[672][2] = -8.1157;
 x[672][3] = 7.0858;
 x[672][4] = -1.2112;
 y[672][0] = b;

 x[673][0] = 0.001;
 x[673][1] = -1.7322;
 x[673][2] = -9.2828;
 x[673][3] = 7.719;
 x[673][4] = -1.7168;
 y[673][0] = b;

 x[674][0] = 0.001;
 x[674][1] = 0.55298;
 x[674][2] = -3.4619;
 x[674][3] = 1.7048;
 x[674][4] = 1.1008;
 y[674][0] = b;

 x[675][0] = 0.001;
 x[675][1] = 2.031;
 x[675][2] = 1.852;
 x[675][3] = -3.0121;
 x[675][4] = 0.003003;
 y[675][0] = b;

 x[676][0] = 0.001;
 x[676][1] = 1.2279;
 x[676][2] = 4.0309;
 x[676][3] = -4.6435;
 x[676][4] = -3.9125;
 y[676][0] = b;

 x[677][0] = 0.001;
 x[677][1] = -4.2249;
 x[677][2] = 6.2699;
 x[677][3] = 0.15822;
 x[677][4] = -5.5457;
 y[677][0] = b;

 x[678][0] = 0.001;
 x[678][1] = -2.5346;
 x[678][2] = -0.77392;
 x[678][3] = 3.3602;
 x[678][4] = 0.00171;
 y[678][0] = b;

 x[679][0] = 0.001;
 x[679][1] = -1.749;
 x[679][2] = -6.332;
 x[679][3] = 6.0987;
 x[679][4] = 0.14266;
 y[679][0] = b;

 x[680][0] = 0.001;
 x[680][1] = -0.539;
 x[680][2] = -5.167;
 x[680][3] = 3.4399;
 x[680][4] = 0.052141;
 y[680][0] = b;

 x[681][0] = 0.001;
 x[681][1] = 1.5631;
 x[681][2] = 0.89599;
 x[681][3] = -1.9702;
 x[681][4] = 0.65472;
 y[681][0] = b;

 x[682][0] = 0.001;
 x[682][1] = 2.3917;
 x[682][2] = 4.5565;
 x[682][3] = -4.9888;
 x[682][4] = -2.8987;
 y[682][0] = b;

 x[683][0] = 0.001;
 x[683][1] = 0.89512;
 x[683][2] = 4.7738;
 x[683][3] = -4.8431;
 x[683][4] = -5.5909;
 y[683][0] = b;

 x[684][0] = 0.001;
 x[684][1] = -5.4808;
 x[684][2] = 8.1819;
 x[684][3] = 0.27818;
 x[684][4] = -5.0323;
 y[684][0] = b;

 x[685][0] = 0.001;
 x[685][1] = -2.8833;
 x[685][2] = 1.7713;
 x[685][3] = 0.68946;
 x[685][4] = -0.4638;
 y[685][0] = b;

 x[686][0] = 0.001;
 x[686][1] = -1.4174;
 x[686][2] = -2.2535;
 x[686][3] = 1.518;
 x[686][4] = 0.61981;
 y[686][0] = b;

 x[687][0] = 0.001;
 x[687][1] = 0.4283;
 x[687][2] = -0.94981;
 x[687][3] = -1.0731;
 x[687][4] = 0.3211;
 y[687][0] = b;

 x[688][0] = 0.001;
 x[688][1] = 1.5904;
 x[688][2] = 2.2121;
 x[688][3] = -3.1183;
 x[688][4] = -0.11725;
 y[688][0] = b;

 x[689][0] = 0.001;
 x[689][1] = 1.7425;
 x[689][2] = 3.6833;
 x[689][3] = -4.0129;
 x[689][4] = -1.7207;
 y[689][0] = b;

 x[690][0] = 0.001;
 x[690][1] = -0.23356;
 x[690][2] = 3.2405;
 x[690][3] = -3.0669;
 x[690][4] = -2.7784;
 y[690][0] = b;

 x[691][0] = 0.001;
 x[691][1] = -3.6227;
 x[691][2] = 3.9958;
 x[691][3] = -0.35845;
 x[691][4] = -3.9047;
 y[691][0] = b;

 x[692][0] = 0.001;
 x[692][1] = -6.1536;
 x[692][2] = 7.9295;
 x[692][3] = 0.61663;
 x[692][4] = -3.2646;
 y[692][0] = b;

 x[693][0] = 0.001;
 x[693][1] = -3.9172;
 x[693][2] = 2.6652;
 x[693][3] = 0.78886;
 x[693][4] = -0.7819;
 y[693][0] = b;

 x[694][0] = 0.001;
 x[694][1] = -2.2214;
 x[694][2] = -0.23798;
 x[694][3] = 0.56008;
 x[694][4] = 0.05602;
 y[694][0] = b;

 x[695][0] = 0.001;
 x[695][1] = -0.49241;
 x[695][2] = 0.89392;
 x[695][3] = -1.6283;
 x[695][4] = -0.56854;
 y[695][0] = b;

 x[696][0] = 0.001;
 x[696][1] = 0.26517;
 x[696][2] = 2.4066;
 x[696][3] = -2.8416;
 x[696][4] = -0.59958;
 y[696][0] = b;

 x[697][0] = 0.001;
 x[697][1] = -0.10234;
 x[697][2] = 1.8189;
 x[697][3] = -2.2169;
 x[697][4] = -0.56725;
 y[697][0] = b;

 x[698][0] = 0.001;
 x[698][1] = -1.6176;
 x[698][2] = 1.0926;
 x[698][3] = -0.35502;
 x[698][4] = -0.59958;
 y[698][0] = b;

 x[699][0] = 0.001;
 x[699][1] = -1.8448;
 x[699][2] = 1.254;
 x[699][3] = 0.27218;
 x[699][4] = -1.0728;
 y[699][0] = b;

 x[700][0] = 0.001;
 x[700][1] = -1.2786;
 x[700][2] = -2.4087;
 x[700][3] = 4.5735;
 x[700][4] = 0.47627;
 y[700][0] = b;

 x[701][0] = 0.001;
 x[701][1] = -2.902;
 x[701][2] = -7.6563;
 x[701][3] = 11.8318;
 x[701][4] = -0.84268;
 y[701][0] = b;

 x[702][0] = 0.001;
 x[702][1] = -4.3773;
 x[702][2] = -5.5167;
 x[702][3] = 10.939;
 x[702][4] = -0.4082;
 y[702][0] = b;

 x[703][0] = 0.001;
 x[703][1] = -2.0529;
 x[703][2] = 3.8385;
 x[703][3] = -0.79544;
 x[703][4] = -1.2138;
 y[703][0] = b;

 x[704][0] = 0.001;
 x[704][1] = 0.18868;
 x[704][2] = 0.70148;
 x[704][3] = -0.51182;
 x[704][4] = 0.005589;
 y[704][0] = b;

 x[705][0] = 0.001;
 x[705][1] = -1.7279;
 x[705][2] = -6.841;
 x[705][3] = 8.9494;
 x[705][4] = 0.68058;
 y[705][0] = b;

 x[706][0] = 0.001;
 x[706][1] = -3.3793;
 x[706][2] = -13.7731;
 x[706][3] = 17.9274;
 x[706][4] = -2.0323;
 y[706][0] = b;

 x[707][0] = 0.001;
 x[707][1] = -3.1273;
 x[707][2] = -7.1121;
 x[707][3] = 11.3897;
 x[707][4] = -0.083634;
 y[707][0] = b;

 x[708][0] = 0.001;
 x[708][1] = -2.121;
 x[708][2] = -0.05588;
 x[708][3] = 1.949;
 x[708][4] = 1.353;
 y[708][0] = b;

 x[709][0] = 0.001;
 x[709][1] = -1.7697;
 x[709][2] = 3.4329;
 x[709][3] = -1.2144;
 x[709][4] = -2.3789;
 y[709][0] = b;

 x[710][0] = 0.001;
 x[710][1] = -0.001285;
 x[710][2] = 0.13863;
 x[710][3] = -0.19651;
 x[710][4] = 0.008175;
 y[710][0] = b;

 x[711][0] = 0.001;
 x[711][1] = -1.682;
 x[711][2] = -6.8121;
 x[711][3] = 7.1398;
 x[711][4] = 1.3323;
 y[711][0] = b;

 x[712][0] = 0.001;
 x[712][1] = -3.4917;
 x[712][2] = -12.1736;
 x[712][3] = 14.3689;
 x[712][4] = -0.61639;
 y[712][0] = b;

 x[713][0] = 0.001;
 x[713][1] = -3.1158;
 x[713][2] = -8.6289;
 x[713][3] = 10.4403;
 x[713][4] = 0.97153;
 y[713][0] = b;

 x[714][0] = 0.001;
 x[714][1] = -2.0891;
 x[714][2] = -0.48422;
 x[714][3] = 1.704;
 x[714][4] = 1.7435;
 y[714][0] = b;

 x[715][0] = 0.001;
 x[715][1] = -1.6936;
 x[715][2] = 2.7852;
 x[715][3] = -2.1835;
 x[715][4] = -1.9276;
 y[715][0] = b;

 x[716][0] = 0.001;
 x[716][1] = -1.2846;
 x[716][2] = 3.2715;
 x[716][3] = -1.7671;
 x[716][4] = -3.2608;
 y[716][0] = b;

 x[717][0] = 0.001;
 x[717][1] = -0.092194;
 x[717][2] = 0.39315;
 x[717][3] = -0.32846;
 x[717][4] = -0.13794;
 y[717][0] = b;

 x[718][0] = 0.001;
 x[718][1] = -1.0292;
 x[718][2] = -6.3879;
 x[718][3] = 5.5255;
 x[718][4] = 0.79955;
 y[718][0] = b;

 x[719][0] = 0.001;
 x[719][1] = -2.2083;
 x[719][2] = -9.1069;
 x[719][3] = 8.9991;
 x[719][4] = -0.28406;
 y[719][0] = b;

 x[720][0] = 0.001;
 x[720][1] = -1.0744;
 x[720][2] = -6.3113;
 x[720][3] = 5.355;
 x[720][4] = 0.80472;
 y[720][0] = b;

 x[721][0] = 0.001;
 x[721][1] = -0.51003;
 x[721][2] = -0.23591;
 x[721][3] = 0.020273;
 x[721][4] = 0.76334;
 y[721][0] = b;

 x[722][0] = 0.001;
 x[722][1] = -0.36372;
 x[722][2] = 3.0439;
 x[722][3] = -3.4816;
 x[722][4] = -2.7836;
 y[722][0] = b;

 x[723][0] = 0.001;
 x[723][1] = -6.3979;
 x[723][2] = 6.4479;
 x[723][3] = 1.0836;
 x[723][4] = -6.6176;
 y[723][0] = b;

 x[724][0] = 0.001;
 x[724][1] = -2.2501;
 x[724][2] = 3.3129;
 x[724][3] = -0.88369;
 x[724][4] = -2.8974;
 y[724][0] = b;

 x[725][0] = 0.001;
 x[725][1] = -1.1859;
 x[725][2] = -1.2519;
 x[725][3] = 2.2635;
 x[725][4] = 0.77239;
 y[725][0] = b;

 x[726][0] = 0.001;
 x[726][1] = -1.8076;
 x[726][2] = -8.8131;
 x[726][3] = 8.7086;
 x[726][4] = -0.21682;
 y[726][0] = b;

 x[727][0] = 0.001;
 x[727][1] = -3.3863;
 x[727][2] = -12.9889;
 x[727][3] = 13.0545;
 x[727][4] = -2.7202;
 y[727][0] = b;

 x[728][0] = 0.001;
 x[728][1] = -1.4106;
 x[728][2] = -7.108;
 x[728][3] = 5.6454;
 x[728][4] = 0.31335;
 y[728][0] = b;

 x[729][0] = 0.001;
 x[729][1] = -0.21394;
 x[729][2] = -0.68287;
 x[729][3] = 0.096532;
 x[729][4] = 1.1965;
 y[729][0] = b;

 x[730][0] = 0.001;
 x[730][1] = 0.48797;
 x[730][2] = 3.5674;
 x[730][3] = -4.3882;
 x[730][4] = -3.8116;
 y[730][0] = b;

 x[731][0] = 0.001;
 x[731][1] = -3.8167;
 x[731][2] = 5.1401;
 x[731][3] = -0.65063;
 x[731][4] = -5.4306;
 y[731][0] = b;

 x[732][0] = 0.001;
 x[732][1] = -1.9555;
 x[732][2] = 0.20692;
 x[732][3] = 1.2473;
 x[732][4] = -0.3707;
 y[732][0] = b;

 x[733][0] = 0.001;
 x[733][1] = -2.1786;
 x[733][2] = -6.4479;
 x[733][3] = 6.0344;
 x[733][4] = -0.20777;
 y[733][0] = b;

 x[734][0] = 0.001;
 x[734][1] = -2.3299;
 x[734][2] = -9.9532;
 x[734][3] = 8.4756;
 x[734][4] = -1.8733;
 y[734][0] = b;

 x[735][0] = 0.001;
 x[735][1] = 0.00312;
 x[735][2] = -4.0061;
 x[735][3] = 1.7956;
 x[735][4] = 0.91722;
 y[735][0] = b;

 x[736][0] = 0.001;
 x[736][1] = 1.3518;
 x[736][2] = 1.0595;
 x[736][3] = -2.3437;
 x[736][4] = 0.39998;
 y[736][0] = b;

 x[737][0] = 0.001;
 x[737][1] = 1.2309;
 x[737][2] = 3.8923;
 x[737][3] = -4.8277;
 x[737][4] = -4.0069;
 y[737][0] = b;

 x[738][0] = 0.001;
 x[738][1] = -5.0301;
 x[738][2] = 7.5032;
 x[738][3] = -0.13396;
 x[738][4] = -7.5034;
 y[738][0] = b;

 x[739][0] = 0.001;
 x[739][1] = -3.0799;
 x[739][2] = 0.60836;
 x[739][3] = 2.7039;
 x[739][4] = -0.23751;
 y[739][0] = b;

 x[740][0] = 0.001;
 x[740][1] = -2.2987;
 x[740][2] = -5.227;
 x[740][3] = 5.63;
 x[740][4] = 0.91722;
 y[740][0] = b;

 x[741][0] = 0.001;
 x[741][1] = -1.239;
 x[741][2] = -6.541;
 x[741][3] = 4.8151;
 x[741][4] = -0.033204;
 y[741][0] = b;

 x[742][0] = 0.001;
 x[742][1] = 0.75896;
 x[742][2] = 0.29176;
 x[742][3] = -1.6506;
 x[742][4] = 0.83834;
 y[742][0] = b;

 x[743][0] = 0.001;
 x[743][1] = 1.6799;
 x[743][2] = 4.2068;
 x[743][3] = -4.5398;
 x[743][4] = -2.3931;
 y[743][0] = b;

 x[744][0] = 0.001;
 x[744][1] = 0.63655;
 x[744][2] = 5.2022;
 x[744][3] = -5.2159;
 x[744][4] = -6.1211;
 y[744][0] = b;

 x[745][0] = 0.001;
 x[745][1] = -6.0598;
 x[745][2] = 9.2952;
 x[745][3] = -0.43642;
 x[745][4] = -6.3694;
 y[745][0] = b;

 x[746][0] = 0.001;
 x[746][1] = -3.518;
 x[746][2] = 2.8763;
 x[746][3] = 0.1548;
 x[746][4] = -1.2086;
 y[746][0] = b;

 x[747][0] = 0.001;
 x[747][1] = -2.0336;
 x[747][2] = -1.4092;
 x[747][3] = 1.1582;
 x[747][4] = 0.36507;
 y[747][0] = b;

 x[748][0] = 0.001;
 x[748][1] = -0.69745;
 x[748][2] = -1.7672;
 x[748][3] = -0.34474;
 x[748][4] = -0.12372;
 y[748][0] = b;

 x[749][0] = 0.001;
 x[749][1] = 0.75108;
 x[749][2] = 1.9161;
 x[749][3] = -3.1098;
 x[749][4] = -0.20518;
 y[749][0] = b;

 x[750][0] = 0.001;
 x[750][1] = 0.84546;
 x[750][2] = 3.4826;
 x[750][3] = -3.6307;
 x[750][4] = -1.3961;
 y[750][0] = b;

 x[751][0] = 0.001;
 x[751][1] = -0.55648;
 x[751][2] = 3.2136;
 x[751][3] = -3.3085;
 x[751][4] = -2.7965;
 y[751][0] = b;

 x[752][0] = 0.001;
 x[752][1] = -3.6817;
 x[752][2] = 3.2239;
 x[752][3] = -0.69347;
 x[752][4] = -3.4004;
 y[752][0] = b;

 x[753][0] = 0.001;
 x[753][1] = -6.7526;
 x[753][2] = 8.8172;
 x[753][3] = -0.061983;
 x[753][4] = -3.725;
 y[753][0] = b;

 x[754][0] = 0.001;
 x[754][1] = -4.577;
 x[754][2] = 3.4515;
 x[754][3] = 0.66719;
 x[754][4] = -0.94742;
 y[754][0] = b;

 x[755][0] = 0.001;
 x[755][1] = -2.9883;
 x[755][2] = 0.31245;
 x[755][3] = 0.45041;
 x[755][4] = 0.068951;
 y[755][0] = b;

 x[756][0] = 0.001;
 x[756][1] = -1.4781;
 x[756][2] = 0.14277;
 x[756][3] = -1.1622;
 x[756][4] = -0.48579;
 y[756][0] = b;

 x[757][0] = 0.001;
 x[757][1] = -0.46651;
 x[757][2] = 2.3383;
 x[757][3] = -2.9812;
 x[757][4] = -1.0431;
 y[757][0] = b;

 x[758][0] = 0.001;
 x[758][1] = -0.8734;
 x[758][2] = 1.6533;
 x[758][3] = -2.1964;
 x[758][4] = -0.78061;
 y[758][0] = b;

 x[759][0] = 0.001;
 x[759][1] = -2.1234;
 x[759][2] = 1.1815;
 x[759][3] = -0.55552;
 x[759][4] = -0.81165;
 y[759][0] = b;

 x[760][0] = 0.001;
 x[760][1] = -2.3142;
 x[760][2] = 2.0838;
 x[760][3] = -0.46813;
 x[760][4] = -1.6767;
 y[760][0] = b;

 x[761][0] = 0.001;
 x[761][1] = -1.4233;
 x[761][2] = -0.98912;
 x[761][3] = 2.3586;
 x[761][4] = 0.39481;
 y[761][0] = b;

 x[762][0] = 0.001;
 x[762][1] = -3.0866;
 x[762][2] = -6.6362;
 x[762][3] = 10.5405;
 x[762][4] = -0.89182;
 y[762][0] = b;

 x[763][0] = 0.001;
 x[763][1] = -4.7331;
 x[763][2] = -6.1789;
 x[763][3] = 11.388;
 x[763][4] = -1.0741;
 y[763][0] = b;

 x[764][0] = 0.001;
 x[764][1] = -2.8829;
 x[764][2] = 3.8964;
 x[764][3] = -0.1888;
 x[764][4] = -1.1672;
 y[764][0] = b;

 x[765][0] = 0.001;
 x[765][1] = -0.036127;
 x[765][2] = 1.525;
 x[765][3] = -1.4089;
 x[765][4] = -0.76121;
 y[765][0] = b;

 x[766][0] = 0.001;
 x[766][1] = -1.7104;
 x[766][2] = -4.778;
 x[766][3] = 6.2109;
 x[766][4] = 0.3974;
 y[766][0] = b;

 x[767][0] = 0.001;
 x[767][1] = -3.8203;
 x[767][2] = -13.0551;
 x[767][3] = 16.9583;
 x[767][4] = -2.3052;
 y[767][0] = b;

 x[768][0] = 0.001;
 x[768][1] = -3.7181;
 x[768][2] = -8.5089;
 x[768][3] = 12.363;
 x[768][4] = -0.95518;
 y[768][0] = b;

 x[769][0] = 0.001;
 x[769][1] = -2.899;
 x[769][2] = -0.60424;
 x[769][3] = 2.6045;
 x[769][4] = 1.3776;
 y[769][0] = b;

 x[770][0] = 0.001;
 x[770][1] = -0.98193;
 x[770][2] = 2.7956;
 x[770][3] = -1.2341;
 x[770][4] = -1.5668;
 y[770][0] = b;

 x[771][0] = 0.001;
 x[771][1] = -0.17296;
 x[771][2] = -1.1816;
 x[771][3] = 1.3818;
 x[771][4] = 0.7336;
 y[771][0] = b;

 x[772][0] = 0.001;
 x[772][1] = -1.9409;
 x[772][2] = -8.6848;
 x[772][3] = 9.155;
 x[772][4] = 0.94049;
 y[772][0] = b;

 x[773][0] = 0.001;
 x[773][1] = -3.5713;
 x[773][2] = -12.4922;
 x[773][3] = 14.8881;
 x[773][4] = -0.47027;
 y[773][0] = b;

 x[774][0] = 0.001;
 x[774][1] = -2.9915;
 x[774][2] = -6.6258;
 x[774][3] = 8.6521;
 x[774][4] = 1.8198;
 y[774][0] = b;

 x[775][0] = 0.001;
 x[775][1] = -1.8483;
 x[775][2] = 0.31038;
 x[775][3] = 0.77344;
 x[775][4] = 1.4189;
 y[775][0] = b;

 x[776][0] = 0.001;
 x[776][1] = -2.2677;
 x[776][2] = 3.2964;
 x[776][3] = -2.2563;
 x[776][4] = -2.4642;
 y[776][0] = b;

 x[777][0] = 0.001;
 x[777][1] = -0.50816;
 x[777][2] = 2.868;
 x[777][3] = -1.8108;
 x[777][4] = -2.2612;
 y[777][0] = b;

 x[778][0] = 0.001;
 x[778][1] = 0.14329;
 x[778][2] = -1.0885;
 x[778][3] = 1.0039;
 x[778][4] = 0.48791;
 y[778][0] = b;

 x[779][0] = 0.001;
 x[779][1] = -0.77688;
 x[779][2] = 0.13036;
 x[779][3] = -0.031137;
 x[779][4] = -0.35389;
 y[779][0] = b;

 x[780][0] = 0.001;
 x[780][1] = -2.7083;
 x[780][2] = -6.8266;
 x[780][3] = 7.5339;
 x[780][4] = 0.59007;
 y[780][0] = b;

 x[781][0] = 0.001;
 x[781][1] = -4.5531;
 x[781][2] = -12.5854;
 x[781][3] = 15.4417;
 x[781][4] = -1.4983;
 y[781][0] = b;

 x[782][0] = 0.001;
 x[782][1] = -3.8894;
 x[782][2] = -7.8322;
 x[782][3] = 9.8208;
 x[782][4] = 0.47498;
 y[782][0] = b;

 x[783][0] = 0.001;
 x[783][1] = -2.5084;
 x[783][2] = -0.22763;
 x[783][3] = 1.488;
 x[783][4] = 1.2069;
 y[783][0] = b;

 x[784][0] = 0.001;
 x[784][1] = -2.1652;
 x[784][2] = 3.0211;
 x[784][3] = -2.4132;
 x[784][4] = -2.4241;
 y[784][0] = b;

 x[785][0] = 0.001;
 x[785][1] = -1.8974;
 x[785][2] = 3.5074;
 x[785][3] = -1.7842;
 x[785][4] = -3.8491;
 y[785][0] = b;

 x[786][0] = 0.001;
 x[786][1] = -0.62043;
 x[786][2] = 0.5587;
 x[786][3] = -0.38587;
 x[786][4] = -0.66423;
 y[786][0] = b;

 x[787][0] = 0.001;
 x[787][1] = -1.8387;
 x[787][2] = -6.301;
 x[787][3] = 5.6506;
 x[787][4] = 0.19567;
 y[787][0] = b;

 x[788][0] = 0.001;
 x[788][1] = -3;
 x[788][2] = -9.1566;
 x[788][3] = 9.5766;
 x[788][4] = -0.73018;
 y[788][0] = b;

 x[789][0] = 0.001;
 x[789][1] = -1.9116;
 x[789][2] = -6.1603;
 x[789][3] = 5.606;
 x[789][4] = 0.48533;
 y[789][0] = b;

 x[790][0] = 0.001;
 x[790][1] = -1.005;
 x[790][2] = 0.084831;
 x[790][3] = -0.2462;
 x[790][4] = 0.45688;
 y[790][0] = b;

 x[791][0] = 0.001;
 x[791][1] = -0.87834;
 x[791][2] = 3.257;
 x[791][3] = -3.6778;
 x[791][4] = -3.2944;
 y[791][0] = b;

 x[792][0] = 0.001;
 x[792][1] = -6.651;
 x[792][2] = 6.7934;
 x[792][3] = 0.68604;
 x[792][4] = -7.5887;
 y[792][0] = b;

 x[793][0] = 0.001;
 x[793][1] = -2.5463;
 x[793][2] = 3.1101;
 x[793][3] = -0.83228;
 x[793][4] = -3.0358;
 y[793][0] = b;

 x[794][0] = 0.001;
 x[794][1] = -1.4377;
 x[794][2] = -1.432;
 x[794][3] = 2.1144;
 x[794][4] = 0.42067;
 y[794][0] = b;

 x[795][0] = 0.001;
 x[795][1] = -2.4554;
 x[795][2] = -9.0407;
 x[795][3] = 8.862;
 x[795][4] = -0.86983;
 y[795][0] = b;

 x[796][0] = 0.001;
 x[796][1] = -3.9411;
 x[796][2] = -12.8792;
 x[796][3] = 13.0597;
 x[796][4] = -3.3125;
 y[796][0] = b;

 x[797][0] = 0.001;
 x[797][1] = -2.1241;
 x[797][2] = -6.8969;
 x[797][3] = 5.5992;
 x[797][4] = -0.47156;
 y[797][0] = b;

 x[798][0] = 0.001;
 x[798][1] = -0.74324;
 x[798][2] = -0.32902;
 x[798][3] = -0.42785;
 x[798][4] = 0.23317;
 y[798][0] = b;

 x[799][0] = 0.001;
 x[799][1] = -0.071503;
 x[799][2] = 3.7412;
 x[799][3] = -4.5415;
 x[799][4] = -4.2526;
 y[799][0] = b;

 x[800][0] = 0.001;
 x[800][1] = -4.2333;
 x[800][2] = 4.9166;
 x[800][3] = -0.49212;
 x[800][4] = -5.3207;
 y[800][0] = b;

 x[801][0] = 0.001;
 x[801][1] = -2.3675;
 x[801][2] = -0.43663;
 x[801][3] = 1.692;
 x[801][4] = -0.43018;
 y[801][0] = b;

 x[802][0] = 0.001;
 x[802][1] = -2.5526;
 x[802][2] = -7.3625;
 x[802][3] = 6.9255;
 x[802][4] = -0.66811;
 y[802][0] = b;

 x[803][0] = 0.001;
 x[803][1] = -3.0986;
 x[803][2] = -10.4602;
 x[803][3] = 8.9717;
 x[803][4] = -2.3427;
 y[803][0] = b;

 x[804][0] = 0.001;
 x[804][1] = -0.89809;
 x[804][2] = -4.4862;
 x[804][3] = 2.2009;
 x[804][4] = 0.50731;
 y[804][0] = b;

 x[805][0] = 0.001;
 x[805][1] = 0.56232;
 x[805][2] = 1.0015;
 x[805][3] = -2.2726;
 x[805][4] = -0.006049;
 y[805][0] = b;

 x[806][0] = 0.001;
 x[806][1] = 0.53936;
 x[806][2] = 3.8944;
 x[806][3] = -4.8166;
 x[806][4] = -4.3418;
 y[806][0] = b;

 x[807][0] = 0.001;
 x[807][1] = -5.3012;
 x[807][2] = 7.3915;
 x[807][3] = 0.029699;
 x[807][4] = -7.3987;
 y[807][0] = b;

 x[808][0] = 0.001;
 x[808][1] = -3.3553;
 x[808][2] = 0.35591;
 x[808][3] = 2.6473;
 x[808][4] = -0.37846;
 y[808][0] = b;

 x[809][0] = 0.001;
 x[809][1] = -2.7908;
 x[809][2] = -5.7133;
 x[809][3] = 5.953;
 x[809][4] = 0.45946;
 y[809][0] = b;

 x[810][0] = 0.001;
 x[810][1] = -1.9983;
 x[810][2] = -6.6072;
 x[810][3] = 4.8254;
 x[810][4] = -0.41984;
 y[810][0] = b;

 x[811][0] = 0.001;
 x[811][1] = 0.15423;
 x[811][2] = 0.11794;
 x[811][3] = -1.6823;
 x[811][4] = 0.59524;
 y[811][0] = b;

 x[812][0] = 0.001;
 x[812][1] = 1.208;
 x[812][2] = 4.0744;
 x[812][3] = -4.7635;
 x[812][4] = -2.6129;
 y[812][0] = b;

 x[813][0] = 0.001;
 x[813][1] = 0.2952;
 x[813][2] = 4.8856;
 x[813][3] = -5.149;
 x[813][4] = -6.2323;
 y[813][0] = b;

 x[814][0] = 0.001;
 x[814][1] = -6.4247;
 x[814][2] = 9.5311;
 x[814][3] = 0.022844;
 x[814][4] = -6.8517;
 y[814][0] = b;

 x[815][0] = 0.001;
 x[815][1] = -3.9933;
 x[815][2] = 2.6218;
 x[815][3] = 0.62863;
 x[815][4] = -1.1595;
 y[815][0] = b;

 x[816][0] = 0.001;
 x[816][1] = -2.659;
 x[816][2] = -1.6058;
 x[816][3] = 1.3647;
 x[816][4] = 0.16464;
 y[816][0] = b;

 x[817][0] = 0.001;
 x[817][1] = -1.4094;
 x[817][2] = -2.1252;
 x[817][3] = -0.10397;
 x[817][4] = -0.19225;
 y[817][0] = b;

 x[818][0] = 0.001;
 x[818][1] = 0.11032;
 x[818][2] = 1.9741;
 x[818][3] = -3.3668;
 x[818][4] = -0.65259;
 y[818][0] = b;

 x[819][0] = 0.001;
 x[819][1] = 0.52374;
 x[819][2] = 3.644;
 x[819][3] = -4.0746;
 x[819][4] = -1.9909;
 y[819][0] = b;

 x[820][0] = 0.001;
 x[820][1] = -0.76794;
 x[820][2] = 3.4598;
 x[820][3] = -3.4405;
 x[820][4] = -3.4276;
 y[820][0] = b;

 x[821][0] = 0.001;
 x[821][1] = -3.9698;
 x[821][2] = 3.6812;
 x[821][3] = -0.60008;
 x[821][4] = -4.0133;
 y[821][0] = b;

 x[822][0] = 0.001;
 x[822][1] = -7.0364;
 x[822][2] = 9.2931;
 x[822][3] = 0.16594;
 x[822][4] = -4.5396;
 y[822][0] = b;

 x[823][0] = 0.001;
 x[823][1] = -4.9447;
 x[823][2] = 3.3005;
 x[823][3] = 1.063;
 x[823][4] = -1.444;
 y[823][0] = b;

 x[824][0] = 0.001;
 x[824][1] = -3.5933;
 x[824][2] = 0.22968;
 x[824][3] = 0.7126;
 x[824][4] = -0.3332;
 y[824][0] = b;

 x[825][0] = 0.001;
 x[825][1] = -2.1674;
 x[825][2] = 0.12415;
 x[825][3] = -1.0465;
 x[825][4] = -0.86208;
 y[825][0] = b;

 x[826][0] = 0.001;
 x[826][1] = -0.9607;
 x[826][2] = 2.6963;
 x[826][3] = -3.1226;
 x[826][4] = -1.3121;
 y[826][0] = b;

 x[827][0] = 0.001;
 x[827][1] = -1.0802;
 x[827][2] = 2.1996;
 x[827][3] = -2.5862;
 x[827][4] = -1.2759;
 y[827][0] = b;

 x[828][0] = 0.001;
 x[828][1] = -2.3277;
 x[828][2] = 1.4381;
 x[828][3] = -0.82114;
 x[828][4] = -1.2862;
 y[828][0] = b;

 x[829][0] = 0.001;
 x[829][1] = -3.7244;
 x[829][2] = 1.9037;
 x[829][3] = -0.035421;
 x[829][4] = -2.5095;
 y[829][0] = b;

 x[830][0] = 0.001;
 x[830][1] = -2.5724;
 x[830][2] = -0.95602;
 x[830][3] = 2.7073;
 x[830][4] = -0.16639;
 y[830][0] = b;

 x[831][0] = 0.001;
 x[831][1] = -3.9297;
 x[831][2] = -6.0816;
 x[831][3] = 10.0958;
 x[831][4] = -1.0147;
 y[831][0] = b;

 x[832][0] = 0.001;
 x[832][1] = -5.2943;
 x[832][2] = -5.1463;
 x[832][3] = 10.3332;
 x[832][4] = -1.1181;
 y[832][0] = b;

 x[833][0] = 0.001;
 x[833][1] = -3.8953;
 x[833][2] = 4.0392;
 x[833][3] = -0.3019;
 x[833][4] = -2.1836;
 y[833][0] = b;

 x[834][0] = 0.001;
 x[834][1] = -1.2244;
 x[834][2] = 1.7485;
 x[834][3] = -1.4801;
 x[834][4] = -1.4181;
 y[834][0] = b;

 x[835][0] = 0.001;
 x[835][1] = -2.6406;
 x[835][2] = -4.4159;
 x[835][3] = 5.983;
 x[835][4] = -0.13924;
 y[835][0] = b;

 x[836][0] = 0.001;
 x[836][1] = -4.6338;
 x[836][2] = -12.7509;
 x[836][3] = 16.7166;
 x[836][4] = -3.2168;
 y[836][0] = b;

 x[837][0] = 0.001;
 x[837][1] = -4.2887;
 x[837][2] = -7.8633;
 x[837][3] = 11.8387;
 x[837][4] = -1.8978;
 y[837][0] = b;

 x[838][0] = 0.001;
 x[838][1] = -3.3458;
 x[838][2] = -0.50491;
 x[838][3] = 2.6328;
 x[838][4] = 0.53705;
 y[838][0] = b;

 x[839][0] = 0.001;
 x[839][1] = -1.1188;
 x[839][2] = 3.3357;
 x[839][3] = -1.3455;
 x[839][4] = -1.9573;
 y[839][0] = b;

 x[840][0] = 0.001;
 x[840][1] = 0.55939;
 x[840][2] = -0.3104;
 x[840][3] = 0.18307;
 x[840][4] = 0.44653;
 y[840][0] = b;

 x[841][0] = 0.001;
 x[841][1] = -1.5078;
 x[841][2] = -7.3191;
 x[841][3] = 7.8981;
 x[841][4] = 1.2289;
 y[841][0] = b;

 x[842][0] = 0.001;
 x[842][1] = -3.506;
 x[842][2] = -12.5667;
 x[842][3] = 15.1606;
 x[842][4] = -0.75216;
 y[842][0] = b;

 x[843][0] = 0.001;
 x[843][1] = -2.9498;
 x[843][2] = -8.273;
 x[843][3] = 10.2646;
 x[843][4] = 1.1629;
 y[843][0] = b;

 x[844][0] = 0.001;
 x[844][1] = -1.6029;
 x[844][2] = -0.38903;
 x[844][3] = 1.62;
 x[844][4] = 1.9103;
 y[844][0] = b;

 x[845][0] = 0.001;
 x[845][1] = -1.2667;
 x[845][2] = 2.8183;
 x[845][3] = -2.426;
 x[845][4] = -1.8862;
 y[845][0] = b;

 x[846][0] = 0.001;
 x[846][1] = -0.49281;
 x[846][2] = 3.0605;
 x[846][3] = -1.8356;
 x[846][4] = -2.834;
 y[846][0] = b;

 x[847][0] = 0.001;
 x[847][1] = 0.66365;
 x[847][2] = -0.045533;
 x[847][3] = -0.18794;
 x[847][4] = 0.23447;
 y[847][0] = b;

 x[848][0] = 0.001;
 x[848][1] = -0.72068;
 x[848][2] = -6.7583;
 x[848][3] = 5.8408;
 x[848][4] = 0.62369;
 y[848][0] = b;

 x[849][0] = 0.001;
 x[849][1] = -1.9966;
 x[849][2] = -9.5001;
 x[849][3] = 9.682;
 x[849][4] = -0.12889;
 y[849][0] = b;

 x[850][0] = 0.001;
 x[850][1] = -0.97325;
 x[850][2] = -6.4168;
 x[850][3] = 5.6026;
 x[850][4] = 1.0323;
 y[850][0] = b;

 x[851][0] = 0.001;
 x[851][1] = -0.025314;
 x[851][2] = -0.17383;
 x[851][3] = -0.11339;
 x[851][4] = 1.2198;
 y[851][0] = b;

 x[852][0] = 0.001;
 x[852][1] = 0.062525;
 x[852][2] = 2.9301;
 x[852][3] = -3.5467;
 x[852][4] = -2.6737;
 y[852][0] = b;

 x[853][0] = 0.001;
 x[853][1] = -5.525;
 x[853][2] = 6.3258;
 x[853][3] = 0.89768;
 x[853][4] = -6.6241;
 y[853][0] = b;

 x[854][0] = 0.001;
 x[854][1] = -1.2943;
 x[854][2] = 2.6735;
 x[854][3] = -0.84085;
 x[854][4] = -2.0323;
 y[854][0] = b;

 x[855][0] = 0.001;
 x[855][1] = -0.24037;
 x[855][2] = -1.7837;
 x[855][3] = 2.135;
 x[855][4] = 1.2418;
 y[855][0] = b;

 x[856][0] = 0.001;
 x[856][1] = -1.3968;
 x[856][2] = -9.6698;
 x[856][3] = 9.4652;
 x[856][4] = -0.34872;
 y[856][0] = b;

 x[857][0] = 0.001;
 x[857][1] = -2.9672;
 x[857][2] = -13.2869;
 x[857][3] = 13.4727;
 x[857][4] = -2.6271;
 y[857][0] = b;

 x[858][0] = 0.001;
 x[858][1] = -1.1005;
 x[858][2] = -7.2508;
 x[858][3] = 6.0139;
 x[858][4] = 0.36895;
 y[858][0] = b;

 x[859][0] = 0.001;
 x[859][1] = 0.22432;
 x[859][2] = -0.52147;
 x[859][3] = -0.40386;
 x[859][4] = 1.2017;
 y[859][0] = b;

 x[860][0] = 0.001;
 x[860][1] = 0.90407;
 x[860][2] = 3.3708;
 x[860][3] = -4.4987;
 x[860][4] = -3.6965;
 y[860][0] = b;

 x[861][0] = 0.001;
 x[861][1] = -2.8619;
 x[861][2] = 4.5193;
 x[861][3] = -0.58123;
 x[861][4] = -4.2629;
 y[861][0] = b;

 x[862][0] = 0.001;
 x[862][1] = -1.0833;
 x[862][2] = -0.31247;
 x[862][3] = 1.2815;
 x[862][4] = 0.41291;
 y[862][0] = b;

 x[863][0] = 0.001;
 x[863][1] = -1.5681;
 x[863][2] = -7.2446;
 x[863][3] = 6.5537;
 x[863][4] = -0.1276;
 y[863][0] = b;

 x[864][0] = 0.001;
 x[864][1] = -2.0545;
 x[864][2] = -10.8679;
 x[864][3] = 9.4926;
 x[864][4] = -1.4116;
 y[864][0] = b;

 x[865][0] = 0.001;
 x[865][1] = 0.2346;
 x[865][2] = -4.5152;
 x[865][3] = 2.1195;
 x[865][4] = 1.4448;
 y[865][0] = b;

 x[866][0] = 0.001;
 x[866][1] = 1.581;
 x[866][2] = 0.86909;
 x[866][3] = -2.3138;
 x[866][4] = 0.82412;
 y[866][0] = b;

 x[867][0] = 0.001;
 x[867][1] = 1.5514;
 x[867][2] = 3.8013;
 x[867][3] = -4.9143;
 x[867][4] = -3.7483;
 y[867][0] = b;

 x[868][0] = 0.001;
 x[868][1] = -4.1479;
 x[868][2] = 7.1225;
 x[868][3] = -0.083404;
 x[868][4] = -6.4172;
 y[868][0] = b;

 x[869][0] = 0.001;
 x[869][1] = -2.2625;
 x[869][2] = -0.099335;
 x[869][3] = 2.8127;
 x[869][4] = 0.48662;
 y[869][0] = b;

 x[870][0] = 0.001;
 x[870][1] = -1.7479;
 x[870][2] = -5.823;
 x[870][3] = 5.8699;
 x[870][4] = 1.212;
 y[870][0] = b;

 x[871][0] = 0.001;
 x[871][1] = -0.95923;
 x[871][2] = -6.7128;
 x[871][3] = 4.9857;
 x[871][4] = 0.32886;
 y[871][0] = b;

 x[872][0] = 0.001;
 x[872][1] = 1.3451;
 x[872][2] = 0.23589;
 x[872][3] = -1.8785;
 x[872][4] = 1.3258;
 y[872][0] = b;

 x[873][0] = 0.001;
 x[873][1] = 2.2279;
 x[873][2] = 4.0951;
 x[873][3] = -4.8037;
 x[873][4] = -2.1112;
 y[873][0] = b;

 x[874][0] = 0.001;
 x[874][1] = 1.2572;
 x[874][2] = 4.8731;
 x[874][3] = -5.2861;
 x[874][4] = -5.8741;
 y[874][0] = b;

 x[875][0] = 0.001;
 x[875][1] = -5.3857;
 x[875][2] = 9.1214;
 x[875][3] = -0.41929;
 x[875][4] = -5.9181;
 y[875][0] = b;

 x[876][0] = 0.001;
 x[876][1] = -2.9786;
 x[876][2] = 2.3445;
 x[876][3] = 0.52667;
 x[876][4] = -0.40173;
 y[876][0] = b;

 x[877][0] = 0.001;
 x[877][1] = -1.5851;
 x[877][2] = -2.1562;
 x[877][3] = 1.7082;
 x[877][4] = 0.9017;
 y[877][0] = b;

 x[878][0] = 0.001;
 x[878][1] = -0.21888;
 x[878][2] = -2.2038;
 x[878][3] = -0.0954;
 x[878][4] = 0.56421;
 y[878][0] = b;

 x[879][0] = 0.001;
 x[879][1] = 1.3183;
 x[879][2] = 1.9017;
 x[879][3] = -3.3111;
 x[879][4] = 0.065071;
 y[879][0] = b;

 x[880][0] = 0.001;
 x[880][1] = 1.4896;
 x[880][2] = 3.4288;
 x[880][3] = -4.0309;
 x[880][4] = -1.4259;
 y[880][0] = b;

 x[881][0] = 0.001;
 x[881][1] = 0.11592;
 x[881][2] = 3.2219;
 x[881][3] = -3.4302;
 x[881][4] = -2.8457;
 y[881][0] = b;

 x[882][0] = 0.001;
 x[882][1] = -3.3924;
 x[882][2] = 3.3564;
 x[882][3] = -0.72004;
 x[882][4] = -3.5233;
 y[882][0] = b;

 x[883][0] = 0.001;
 x[883][1] = -6.1632;
 x[883][2] = 8.7096;
 x[883][3] = -0.21621;
 x[883][4] = -3.6345;
 y[883][0] = b;

 x[884][0] = 0.001;
 x[884][1] = -4.0786;
 x[884][2] = 2.9239;
 x[884][3] = 0.87026;
 x[884][4] = -0.65389;
 y[884][0] = b;

 x[885][0] = 0.001;
 x[885][1] = -2.5899;
 x[885][2] = -0.3911;
 x[885][3] = 0.93452;
 x[885][4] = 0.42972;
 y[885][0] = b;

 x[886][0] = 0.001;
 x[886][1] = -1.0116;
 x[886][2] = -0.19038;
 x[886][3] = -0.90597;
 x[886][4] = 0.003003;
 y[886][0] = b;

 x[887][0] = 0.001;
 x[887][1] = 0.066129;
 x[887][2] = 2.4914;
 x[887][3] = -2.9401;
 x[887][4] = -0.62156;
 y[887][0] = b;

 x[888][0] = 0.001;
 x[888][1] = -0.24745;
 x[888][2] = 1.9368;
 x[888][3] = -2.4697;
 x[888][4] = -0.80518;
 y[888][0] = b;

 x[889][0] = 0.001;
 x[889][1] = -1.5732;
 x[889][2] = 1.0636;
 x[889][3] = -0.71232;
 x[889][4] = -0.8388;
 y[889][0] = b;

 x[890][0] = 0.001;
 x[890][1] = -2.1668;
 x[890][2] = 1.5933;
 x[890][3] = 0.045122;
 x[890][4] = -1.678;
 y[890][0] = b;

 x[891][0] = 0.001;
 x[891][1] = -1.1667;
 x[891][2] = -1.4237;
 x[891][3] = 2.9241;
 x[891][4] = 0.66119;
 y[891][0] = b;

 x[892][0] = 0.001;
 x[892][1] = -2.8391;
 x[892][2] = -6.63;
 x[892][3] = 10.4849;
 x[892][4] = -0.42113;
 y[892][0] = b;

 x[893][0] = 0.001;
 x[893][1] = -4.5046;
 x[893][2] = -5.8126;
 x[893][3] = 10.8867;
 x[893][4] = -0.52846;
 y[893][0] = b;

 x[894][0] = 0.001;
 x[894][1] = -2.41;
 x[894][2] = 3.7433;
 x[894][3] = -0.40215;
 x[894][4] = -1.2953;
 y[894][0] = b;

 x[895][0] = 0.001;
 x[895][1] = 0.40614;
 x[895][2] = 1.3492;
 x[895][3] = -1.4501;
 x[895][4] = -0.55949;
 y[895][0] = b;

 x[896][0] = 0.001;
 x[896][1] = -1.3887;
 x[896][2] = -4.8773;
 x[896][3] = 6.4774;
 x[896][4] = 0.34179;
 y[896][0] = b;

 x[897][0] = 0.001;
 x[897][1] = -3.7503;
 x[897][2] = -13.4586;
 x[897][3] = 17.5932;
 x[897][4] = -2.7771;
 y[897][0] = b;

 x[898][0] = 0.001;
 x[898][1] = -3.5637;
 x[898][2] = -8.3827;
 x[898][3] = 12.393;
 x[898][4] = -1.2823;
 y[898][0] = b;

 x[899][0] = 0.001;
 x[899][1] = -2.5419;
 x[899][2] = -0.65804;
 x[899][3] = 2.6842;
 x[899][4] = 1.1952;
 y[899][0] = b;
}

}

void inicializa_pesos(double w1[][NINTER],double w2[][NSAI])
{
   int i,j,k;
   double aleatorio;
   for(i=0;i<NENT;i++)
     for(j=1;j<NINT;j++)
     {
	    aleatorio = (rand()%10);
		w1[i][j] =  (1.0 - 2.0 * aleatorio) / 100.0;
     }
   for(j=0;j<NINT;j++)
     for(k=0;k<NSAI;k++)
       {
	     aleatorio = (rand()%10);
	     w2[j][k] = (1.0 - 2.0 * aleatorio) / 100.0;
       }
}

void intermediaria(double x[][NENT],double w1[][NINTER],double h[],int m)
  {
   int i,j;
   double somatorio;
   h[0] = 1.0;
   for(j=1;j<NINT;j++)
     {
      somatorio = 0.0;
      for(i=0;i<NENT;i++)
	     somatorio = somatorio + x[m][i] * w1[i][j];
         
      
     h[j] = 1.0 / (1.0 + exp(-somatorio));
       
     }
  }

void saida(double h[],double w2[][NSAI],double o[])
  {
   int j,k;
   double somatorio;
   for(k=0;k<NSAI;k++)
     {
      somatorio = 0.0;
      for(j=0;j<NINT;j++)
    	somatorio = somatorio + h[j] * w2[j][k];
       
      
      o[k] = 1.0 / (1.0 + exp(-somatorio));
      // 1.0 / (1.0 + exp(somatorio));
     }
  }

double erro_saida(double o[],double y[][NSAI],int m)
{
   int k;
   double somatorio,erro;
   somatorio = 0.0;
   for(k=0;k<NSAI;k++)
     somatorio = somatorio + (o[k] - y[m][k]) * (o[k] - y[m][k]);
   erro = 0.5 * somatorio;

   return erro;
}

void erro2(double o[],double y[][NSAI],int m,double delta2[])
  {
   int k;
   for(k=0;k<NSAI;k++)
     delta2[k] = o[k] * (1.0 - o[k]) * (y[m][k] - o[k]);
  }

void erro1(double h[],double delta2[],double w2[][NSAI],double delta1[])
  {
   int j,k;
   double somatorio;
   for(j=1;j<NINT;j++)
     {
      somatorio = 0.0;
      for(k=0;k<NSAI;k++)
    	somatorio = somatorio + delta2[k] * w2[j][k];
      delta1[j] = h[j] * (1 - h[j]) * somatorio;
     }
  }

void ajusta2(double w2[][NSAI],double delta2[],double h[])
  {
   int j,k;
   for(j=0;j<NINT;j++)
     for(k=0;k<NSAI;k++)
       w2[j][k] = w2[j][k] + TAPR3 * delta2[k] * h[j];
   // w2[2][1]=0;
  }

void ajusta1(double w1[][NINTER],double delta1[],double x[][NENT],int m)
  {
   int i,j;
   for(i=0;i<NENT;i++)
     for(j=1;j<NINT;j++)
       w1[i][j] = w1[i][j] + TAPR3 * delta1[j] * x[m][i];
   //w1[1][1]=0;
  }
void verifica1(double x[][NENT],double w1[][NINTER],double w2[][NSAI],double y[][NSAI],FILE *pFile)
{
   int m,i,j;
   double err,h[NINTER],o[NSAI];
   
   float a = 1, b = 0;
   {
   	
//TESTE

 x[0][0] = 0.001;
 x[0][1] = 2.9719;
 x[0][2] = 6.8369;
 x[0][3] = -0.2702;
 x[0][4] = 0.71291;
 y[0][0] = a;

 x[1][0] = 0.001;
 x[1][1] = 1.6849;
 x[1][2] = 8.7489;
 x[1][3] = -1.2641;
 x[1][4] = -1.3858;
 y[1][0] = a;

 x[2][0] = 0.001;
 x[2][1] = -1.9177;
 x[2][2] = 11.6894;
 x[2][3] = 2.5454;
 x[2][4] = -3.2763;
 y[2][0] = a;

 x[3][0] = 0.001;
 x[3][1] = 2.3729;
 x[3][2] = 10.4726;
 x[3][3] = -3.0087;
 x[3][4] = -3.2013;
 y[3][0] = a;

 x[4][0] = 0.001;
 x[4][1] = 1.0284;
 x[4][2] = 9.767;
 x[4][3] = -1.3687;
 x[4][4] = -1.7853;
 y[4][0] = a;

 x[5][0] = 0.001;
 x[5][1] = 0.27451;
 x[5][2] = 9.2186;
 x[5][3] = -3.2863;
 x[5][4] = -4.8448;
 y[5][0] = a;

 x[6][0] = 0.001;
 x[6][1] = 1.6032;
 x[6][2] = -4.7863;
 x[6][3] = 8.5193;
 x[6][4] = -2.1203;
 y[6][0] = a;

 x[7][0] = 0.001;
 x[7][1] = 4.616;
 x[7][2] = 10.1788;
 x[7][3] = -4.2185;
 x[7][4] = -4.4245;
 y[7][0] = a;

 x[8][0] = 0.001;
 x[8][1] = 4.2478;
 x[8][2] = 7.6956;
 x[8][3] = -2.7696;
 x[8][4] = -1.0767;
 y[8][0] = a;

 x[9][0] = 0.001;
 x[9][1] = 4.0215;
 x[9][2] = -2.7004;
 x[9][3] = 2.4957;
 x[9][4] = 0.36636;
 y[9][0] = a;

 x[10][0] = 0.001;
 x[10][1] = 5.0297;
 x[10][2] = -4.9704;
 x[10][3] = 3.5025;
 x[10][4] = -0.23751;
 y[10][0] = a;

 x[11][0] = 0.001;
 x[11][1] = 1.5902;
 x[11][2] = 2.2948;
 x[11][3] = 3.2403;
 x[11][4] = 0.18404;
 y[11][0] = a;

 x[12][0] = 0.001;
 x[12][1] = 2.1274;
 x[12][2] = 5.1939;
 x[12][3] = -1.7971;
 x[12][4] = -1.1763;
 y[12][0] = a;

 x[13][0] = 0.001;
 x[13][1] = 1.1811;
 x[13][2] = 8.3847;
 x[13][3] = -2.0567;
 x[13][4] = -0.90345;
 y[13][0] = a;

 x[14][0] = 0.001;
 x[14][1] = 0.3292;
 x[14][2] = -4.4552;
 x[14][3] = 4.5718;
 x[14][4] = -0.9888;
 y[14][0] = a;

 x[15][0] = 0.001;
 x[15][1] = 5.7353;
 x[15][2] = 5.2808;
 x[15][3] = -2.2598;
 x[15][4] = 0.075416;
 y[15][0] = a;

 x[16][0] = 0.001;
 x[16][1] = 2.6718;
 x[16][2] = 5.6574;
 x[16][3] = 0.72974;
 x[16][4] = -1.4892;
 y[16][0] = a;

 x[17][0] = 0.001;
 x[17][1] = 1.5799;
 x[17][2] = -4.7076;
 x[17][3] = 7.9186;
 x[17][4] = -1.5487;
 y[17][0] = a;

 x[18][0] = 0.001;
 x[18][1] = 2.9499;
 x[18][2] = 2.2493;
 x[18][3] = 1.3458;
 x[18][4] = -0.037083;
 y[18][0] = a;

 x[19][0] = 0.001;
 x[19][1] = 0.5195;
 x[19][2] = -3.2633;
 x[19][3] = 3.0895;
 x[19][4] = -0.9849;
 y[19][0] = a;

 x[20][0] = 0.001;
 x[20][1] = 3.7352;
 x[20][2] = 9.5911;
 x[20][3] = -3.9032;
 x[20][4] = -3.3487;
 y[20][0] = a;

 x[21][0] = 0.001;
 x[21][1] = -1.7344;
 x[21][2] = 2.0175;
 x[21][3] = 7.7618;
 x[21][4] = 0.93532;
 y[21][0] = a;

 x[22][0] = 0.001;
 x[22][1] = 3.884;
 x[22][2] = 10.0277;
 x[22][3] = -3.9298;
 x[22][4] = -4.0819;
 y[22][0] = a;

 x[23][0] = 0.001;
 x[23][1] = 3.5257;
 x[23][2] = 1.2829;
 x[23][3] = 1.9276;
 x[23][4] = 1.7991;
 y[23][0] = a;

 x[24][0] = 0.001;
 x[24][1] = 4.4549;
 x[24][2] = 2.4976;
 x[24][3] = 1.0313;
 x[24][4] = 0.96894;
 y[24][0] = a;

 x[25][0] = 0.001;
 x[25][1] = -0.16108;
 x[25][2] = -6.4624;
 x[25][3] = 8.3573;
 x[25][4] = -1.5216;
 y[25][0] = a;

 x[26][0] = 0.001;
 x[26][1] = 4.2164;
 x[26][2] = 9.4607;
 x[26][3] = -4.9288;
 x[26][4] = -5.2366;
 y[26][0] = a;

 x[27][0] = 0.001;
 x[27][1] = 3.5152;
 x[27][2] = 6.8224;
 x[27][3] = -0.67377;
 x[27][4] = -0.46898;
 y[27][0] = a;

 x[28][0] = 0.001;
 x[28][1] = 1.6988;
 x[28][2] = 2.9094;
 x[28][3] = 2.9044;
 x[28][4] = 0.11033;
 y[28][0] = a;

 x[29][0] = 0.001;
 x[29][1] = 1.0607;
 x[29][2] = 2.4542;
 x[29][3] = 2.5188;
 x[29][4] = -0.17027;
 y[29][0] = a;

 x[30][0] = 0.001;
 x[30][1] = 2.0421;
 x[30][2] = 1.2436;
 x[30][3] = 4.2171;
 x[30][4] = 0.90429;
 y[30][0] = a;

 x[31][0] = 0.001;
 x[31][1] = 3.5594;
 x[31][2] = 1.3078;
 x[31][3] = 1.291;
 x[31][4] = 1.6556;
 y[31][0] = a;

 x[32][0] = 0.001;
 x[32][1] = 3.0009;
 x[32][2] = 5.8126;
 x[32][3] = -2.2306;
 x[32][4] = -0.66553;
 y[32][0] = a;

 x[33][0] = 0.001;
 x[33][1] = 3.9294;
 x[33][2] = 1.4112;
 x[33][3] = 1.8076;
 x[33][4] = 0.89782;
 y[33][0] = a;

 x[34][0] = 0.001;
 x[34][1] = 3.4667;
 x[34][2] = -4.0724;
 x[34][3] = 4.2882;
 x[34][4] = 1.5418;
 y[34][0] = a;

 x[35][0] = 0.001;
 x[35][1] = 3.966;
 x[35][2] = 3.9213;
 x[35][3] = 0.70574;
 x[35][4] = 0.33662;
 y[35][0] = a;

 x[36][0] = 0.001;
 x[36][1] = 1.0191;
 x[36][2] = 2.33;
 x[36][3] = 4.9334;
 x[36][4] = 0.82929;
 y[36][0] = a;

 x[37][0] = 0.001;
 x[37][1] = 0.96414;
 x[37][2] = 5.616;
 x[37][3] = 2.2138;
 x[37][4] = -0.12501;
 y[37][0] = a;

 x[38][0] = 0.001;
 x[38][1] = 1.8205;
 x[38][2] = 6.7562;
 x[38][3] = 0.009991;
 x[38][4] = 0.39481;
 y[38][0] = a;

 x[39][0] = 0.001;
 x[39][1] = 4.9923;
 x[39][2] = 7.8653;
 x[39][3] = -2.3515;
 x[39][4] = -0.71984;
 y[39][0] = a;

 x[40][0] = 0.001;
 x[40][1] = -1.1804;
 x[40][2] = 11.5093;
 x[40][3] = 0.15565;
 x[40][4] = -6.8194;
 y[40][0] = a;

 x[41][0] = 0.001;
 x[41][1] = 4.0329;
 x[41][2] = 0.23175;
 x[41][3] = 0.89082;
 x[41][4] = 1.1823;
 y[41][0] = a;

 x[42][0] = 0.001;
 x[42][1] = 0.66018;
 x[42][2] = 10.3878;
 x[42][3] = -1.4029;
 x[42][4] = -3.9151;
 y[42][0] = a;

 x[43][0] = 0.001;
 x[43][1] = 3.5982;
 x[43][2] = 7.1307;
 x[43][3] = -1.3035;
 x[43][4] = 0.21248;
 y[43][0] = a;

 x[44][0] = 0.001;
 x[44][1] = -1.8584;
 x[44][2] = 7.886;
 x[44][3] = -1.6643;
 x[44][4] = -1.8384;
 y[44][0] = a;

 x[45][0] = 0.001;
 x[45][1] = 4.0972;
 x[45][2] = 0.46972;
 x[45][3] = 1.6671;
 x[45][4] = 0.91593;
 y[45][0] = a;

 x[46][0] = 0.001;
 x[46][1] = 3.3299;
 x[46][2] = 0.91254;
 x[46][3] = 1.5806;
 x[46][4] = 0.39352;
 y[46][0] = a;

 x[47][0] = 0.001;
 x[47][1] = 3.1088;
 x[47][2] = 3.1122;
 x[47][3] = 0.80857;
 x[47][4] = 0.4336;
 y[47][0] = a;

 x[48][0] = 0.001;
 x[48][1] = -4.2859;
 x[48][2] = 8.5234;
 x[48][3] = 3.1392;
 x[48][4] = -0.91639;
 y[48][0] = a;

 x[49][0] = 0.001;
 x[49][1] = -1.2528;
 x[49][2] = 10.2036;
 x[49][3] = 2.1787;
 x[49][4] = -5.6038;
 y[49][0] = a;

 x[50][0] = 0.001;
 x[50][1] = 0.5195;
 x[50][2] = -3.2633;
 x[50][3] = 3.0895;
 x[50][4] = -0.9849;
 y[50][0] = a;

 x[51][0] = 0.001;
 x[51][1] = 0.3292;
 x[51][2] = -4.4552;
 x[51][3] = 4.5718;
 x[51][4] = -0.9888;
 y[51][0] = a;

 x[52][0] = 0.001;
 x[52][1] = 0.88872;
 x[52][2] = 5.3449;
 x[52][3] = 2.045;
 x[52][4] = -0.19355;
 y[52][0] = a;

 x[53][0] = 0.001;
 x[53][1] = 3.5458;
 x[53][2] = 9.3718;
 x[53][3] = -4.0351;
 x[53][4] = -3.9564;
 y[53][0] = a;

 x[54][0] = 0.001;
 x[54][1] = -0.21661;
 x[54][2] = 8.0329;
 x[54][3] = 1.8848;
 x[54][4] = -3.8853;
 y[54][0] = a;

 x[55][0] = 0.001;
 x[55][1] = 2.7206;
 x[55][2] = 9.0821;
 x[55][3] = -3.3111;
 x[55][4] = -0.96811;
 y[55][0] = a;

 x[56][0] = 0.001;
 x[56][1] = 3.2051;
 x[56][2] = 8.6889;
 x[56][3] = -2.9033;
 x[56][4] = -0.7819;
 y[56][0] = a;

 x[57][0] = 0.001;
 x[57][1] = 2.6917;
 x[57][2] = 10.8161;
 x[57][3] = -3.3;
 x[57][4] = -4.2888;
 y[57][0] = a;

 x[58][0] = 0.001;
 x[58][1] = -2.3242;
 x[58][2] = 11.5176;
 x[58][3] = 1.8231;
 x[58][4] = -5.375;
 y[58][0] = a;

 x[59][0] = 0.001;
 x[59][1] = 2.7161;
 x[59][2] = -4.2006;
 x[59][3] = 4.1914;
 x[59][4] = 0.16981;
 y[59][0] = a;

 x[60][0] = 0.001;
 x[60][1] = 3.3848;
 x[60][2] = 3.2674;
 x[60][3] = 0.90967;
 x[60][4] = 0.25128;
 y[60][0] = a;

 x[61][0] = 0.001;
 x[61][1] = 1.7452;
 x[61][2] = 4.8028;
 x[61][3] = 2.0878;
 x[61][4] = 0.62627;
 y[61][0] = a;

 x[62][0] = 0.001;
 x[62][1] = 2.805;
 x[62][2] = 0.57732;
 x[62][3] = 1.3424;
 x[62][4] = 1.2133;
 y[62][0] = a;

 x[63][0] = 0.001;
 x[63][1] = 5.7823;
 x[63][2] = 5.5788;
 x[63][3] = -2.4089;
 x[63][4] = -0.056479;
 y[63][0] = a;

 x[64][0] = 0.001;
 x[64][1] = 3.8999;
 x[64][2] = 1.734;
 x[64][3] = 1.6011;
 x[64][4] = 0.96765;
 y[64][0] = a;

 x[65][0] = 0.001;
 x[65][1] = 3.5189;
 x[65][2] = 6.332;
 x[65][3] = -1.7791;
 x[65][4] = -0.020273;
 y[65][0] = a;

 x[66][0] = 0.001;
 x[66][1] = 3.2294;
 x[66][2] = 7.7391;
 x[66][3] = -0.37816;
 x[66][4] = -2.5405;
 y[66][0] = a;

 x[67][0] = 0.001;
 x[67][1] = 3.4985;
 x[67][2] = 3.1639;
 x[67][3] = 0.22677;
 x[67][4] = -0.1651;
 y[67][0] = a;

 x[68][0] = 0.001;
 x[68][1] = 2.1948;
 x[68][2] = 1.3781;
 x[68][3] = 1.1582;
 x[68][4] = 0.85774;
 y[68][0] = a;

 x[69][0] = 0.001;
 x[69][1] = 2.2526;
 x[69][2] = 9.9636;
 x[69][3] = -3.1749;
 x[69][4] = -2.9944;
 y[69][0] = a;

 x[70][0] = 0.001;
 x[70][1] = 4.1529;
 x[70][2] = -3.9358;
 x[70][3] = 2.8633;
 x[70][4] = -0.017686;
 y[70][0] = a;

 x[71][0] = 0.001;
 x[71][1] = 0.74307;
 x[71][2] = 11.17;
 x[71][3] = -1.3824;
 x[71][4] = -4.0728;
 y[71][0] = a;

 x[72][0] = 0.001;
 x[72][1] = 1.9105;
 x[72][2] = 8.871;
 x[72][3] = -2.3386;
 x[72][4] = -0.75604;
 y[72][0] = a;

 x[73][0] = 0.001;
 x[73][1] = -1.5055;
 x[73][2] = 0.070346;
 x[73][3] = 6.8681;
 x[73][4] = -0.50648;
 y[73][0] = a;

 x[74][0] = 0.001;
 x[74][1] = 0.58836;
 x[74][2] = 10.7727;
 x[74][3] = -1.3884;
 x[74][4] = -4.3276;
 y[74][0] = a;

 x[75][0] = 0.001;
 x[75][1] = 3.2303;
 x[75][2] = 7.8384;
 x[75][3] = -3.5348;
 x[75][4] = -1.2151;
 y[75][0] = a;

 x[76][0] = 0.001;
 x[76][1] = -1.9922;
 x[76][2] = 11.6542;
 x[76][3] = 2.6542;
 x[76][4] = -5.2107;
 y[76][0] = a;

 x[77][0] = 0.001;
 x[77][1] = 2.8523;
 x[77][2] = 9.0096;
 x[77][3] = -3.761;
 x[77][4] = -3.3371;
 y[77][0] = a;

 x[78][0] = 0.001;
 x[78][1] = 4.2772;
 x[78][2] = 2.4955;
 x[78][3] = 0.48554;
 x[78][4] = 0.36119;
 y[78][0] = a;

 x[79][0] = 0.001;
 x[79][1] = 1.5099;
 x[79][2] = 0.039307;
 x[79][3] = 6.2332;
 x[79][4] = -0.30346;
 y[79][0] = a;

 x[80][0] = 0.001;
 x[80][1] = 5.4188;
 x[80][2] = 10.1457;
 x[80][3] = -4.084;
 x[80][4] = -3.6991;
 y[80][0] = a;

 x[81][0] = 0.001;
 x[81][1] = 0.86202;
 x[81][2] = 2.6963;
 x[81][3] = 4.2908;
 x[81][4] = 0.54739;
 y[81][0] = a;

 x[82][0] = 0.001;
 x[82][1] = 3.8117;
 x[82][2] = 10.1457;
 x[82][3] = -4.0463;
 x[82][4] = -4.5629;
 y[82][0] = a;

 x[83][0] = 0.001;
 x[83][1] = 0.54777;
 x[83][2] = 10.3754;
 x[83][3] = -1.5435;
 x[83][4] = -4.1633;
 y[83][0] = a;

 x[84][0] = 0.001;
 x[84][1] = 2.3718;
 x[84][2] = 7.4908;
 x[84][3] = 0.015989;
 x[84][4] = -1.7414;
 y[84][0] = a;

 x[85][0] = 0.001;
 x[85][1] = -2.4953;
 x[85][2] = 11.1472;
 x[85][3] = 1.9353;
 x[85][4] = -3.4638;
 y[85][0] = a;

 x[86][0] = 0.001;
 x[86][1] = 4.6361;
 x[86][2] = -2.6611;
 x[86][3] = 2.8358;
 x[86][4] = 1.1991;
 y[86][0] = a;

 x[87][0] = 0.001;
 x[87][1] = -2.2527;
 x[87][2] = 11.5321;
 x[87][3] = 2.5899;
 x[87][4] = -3.2737;
 y[87][0] = a;

 x[88][0] = 0.001;
 x[88][1] = 3.7982;
 x[88][2] = 10.423;
 x[88][3] = -4.1602;
 x[88][4] = -4.9728;
 y[88][0] = a;

 x[89][0] = 0.001;
 x[89][1] = -0.36279;
 x[89][2] = 8.2895;
 x[89][3] = -1.9213;
 x[89][4] = -3.3332;
 y[89][0] = a;

 x[90][0] = 0.001;
 x[90][1] = 2.1265;
 x[90][2] = 6.8783;
 x[90][3] = 0.44784;
 x[90][4] = -2.2224;
 y[90][0] = a;

 x[91][0] = 0.001;
 x[91][1] = 0.86736;
 x[91][2] = 5.5643;
 x[91][3] = 1.6765;
 x[91][4] = -0.16769;
 y[91][0] = a;

 x[92][0] = 0.001;
 x[92][1] = 3.7831;
 x[92][2] = 10.0526;
 x[92][3] = -3.8869;
 x[92][4] = -3.7366;
 y[92][0] = a;

 x[93][0] = 0.001;
 x[93][1] = -2.2623;
 x[93][2] = 12.1177;
 x[93][3] = 0.28846;
 x[93][4] = -7.7581;
 y[93][0] = a;

 x[94][0] = 0.001;
 x[94][1] = 1.2616;
 x[94][2] = 4.4303;
 x[94][3] = -1.3335;
 x[94][4] = -1.7517;
 y[94][0] = a;

 x[95][0] = 0.001;
 x[95][1] = 2.6799;
 x[95][2] = 3.1349;
 x[95][3] = 0.34073;
 x[95][4] = 0.58489;
 y[95][0] = a;

 x[96][0] = 0.001;
 x[96][1] = -0.39816;
 x[96][2] = 5.9781;
 x[96][3] = 1.3912;
 x[96][4] = -1.1621;
 y[96][0] = a;

 x[97][0] = 0.001;
 x[97][1] = 4.3937;
 x[97][2] = 0.35798;
 x[97][3] = 2.0416;
 x[97][4] = 1.2004;
 y[97][0] = a;

 x[98][0] = 0.001;
 x[98][1] = 2.9695;
 x[98][2] = 5.6222;
 x[98][3] = 0.27561;
 x[98][4] = -1.1556;
 y[98][0] = a;

 x[99][0] = 0.001;
 x[99][1] = 1.3049;
 x[99][2] = -0.15521;
 x[99][3] = 6.4911;
 x[99][4] = -0.75346;
 y[99][0] = a;

 x[100][0] = 0.001;
 x[100][1] = 2.2123;
 x[100][2] = -5.8395;
 x[100][3] = 7.7687;
 x[100][4] = -0.85302;
 y[100][0] = a;

 x[101][0] = 0.001;
 x[101][1] = 1.9647;
 x[101][2] = 6.9383;
 x[101][3] = 0.57722;
 x[101][4] = 0.66377;
 y[101][0] = a;

 x[102][0] = 0.001;
 x[102][1] = 3.0864;
 x[102][2] = -2.5845;
 x[102][3] = 2.2309;
 x[102][4] = 0.30947;
 y[102][0] = a;

 x[103][0] = 0.001;
 x[103][1] = 0.3798;
 x[103][2] = 0.7098;
 x[103][3] = 0.7572;
 x[103][4] = -0.4444;
 y[103][0] = a;

 x[104][0] = 0.001;
 x[104][1] = 0.58982;
 x[104][2] = 7.4266;
 x[104][3] = 1.2353;
 x[104][4] = -2.9595;
 y[104][0] = a;

 x[105][0] = 0.001;
 x[105][1] = 0.14783;
 x[105][2] = 7.946;
 x[105][3] = 1.0742;
 x[105][4] = -3.3409;
 y[105][0] = a;

 x[106][0] = 0.001;
 x[106][1] = -0.062025;
 x[106][2] = 6.1975;
 x[106][3] = 1.099;
 x[106][4] = -1.131;
 y[106][0] = a;

 x[107][0] = 0.001;
 x[107][1] = 4.223;
 x[107][2] = 1.1319;
 x[107][3] = 0.72202;
 x[107][4] = 0.96118;
 y[107][0] = a;

 x[108][0] = 0.001;
 x[108][1] = 0.64295;
 x[108][2] = 7.1018;
 x[108][3] = 0.3493;
 x[108][4] = -0.41337;
 y[108][0] = a;

 x[109][0] = 0.001;
 x[109][1] = 1.941;
 x[109][2] = 0.46351;
 x[109][3] = 4.6472;
 x[109][4] = 1.0879;
 y[109][0] = a;

 x[110][0] = 0.001;
 x[110][1] = 4.0047;
 x[110][2] = 0.45937;
 x[110][3] = 1.3621;
 x[110][4] = 1.6181;
 y[110][0] = a;

 x[111][0] = 0.001;
 x[111][1] = 3.7767;
 x[111][2] = 9.7794;
 x[111][3] = -3.9075;
 x[111][4] = -3.5323;
 y[111][0] = a;

 x[112][0] = 0.001;
 x[112][1] = 3.4769;
 x[112][2] = -0.15314;
 x[112][3] = 2.53;
 x[112][4] = 2.4495;
 y[112][0] = a;

 x[113][0] = 0.001;
 x[113][1] = 1.9818;
 x[113][2] = 9.2621;
 x[113][3] = -3.521;
 x[113][4] = -1.872;
 y[113][0] = a;

 x[114][0] = 0.001;
 x[114][1] = 3.8023;
 x[114][2] = -3.8696;
 x[114][3] = 4.044;
 x[114][4] = 0.95343;
 y[114][0] = a;

 x[115][0] = 0.001;
 x[115][1] = 4.3483;
 x[115][2] = 11.1079;
 x[115][3] = -4.0857;
 x[115][4] = -4.2539;
 y[115][0] = a;

 x[116][0] = 0.001;
 x[116][1] = 1.1518;
 x[116][2] = 1.3864;
 x[116][3] = 5.2727;
 x[116][4] = -0.43536;
 y[116][0] = a;

 x[117][0] = 0.001;
 x[117][1] = -1.2576;
 x[117][2] = 1.5892;
 x[117][3] = 7.0078;
 x[117][4] = 0.42455;
 y[117][0] = a;

 x[118][0] = 0.001;
 x[118][1] = 1.9572;
 x[118][2] = -5.1153;
 x[118][3] = 8.6127;
 x[118][4] = -1.4297;
 y[118][0] = a;

 x[119][0] = 0.001;
 x[119][1] = -2.484;
 x[119][2] = 12.1611;
 x[119][3] = 2.8204;
 x[119][4] = -3.7418;
 y[119][0] = a;

 x[120][0] = 0.001;
 x[120][1] = -1.1497;
 x[120][2] = 1.2954;
 x[120][3] = 7.701;
 x[120][4] = 0.62627;
 y[120][0] = a;

 x[121][0] = 0.001;
 x[121][1] = 4.8368;
 x[121][2] = 10.0132;
 x[121][3] = -4.3239;
 x[121][4] = -4.3276;
 y[121][0] = a;

 x[122][0] = 0.001;
 x[122][1] = -0.12196;
 x[122][2] = 8.8068;
 x[122][3] = 0.94566;
 x[122][4] = -4.2267;
 y[122][0] = a;

 x[123][0] = 0.001;
 x[123][1] = 1.9429;
 x[123][2] = 6.3961;
 x[123][3] = 0.092248;
 x[123][4] = 0.58102;
 y[123][0] = a;

 x[124][0] = 0.001;
 x[124][1] = 1.742;
 x[124][2] = -4.809;
 x[124][3] = 8.2142;
 x[124][4] = -2.0659;
 y[124][0] = a;

 x[125][0] = 0.001;
 x[125][1] = -1.5222;
 x[125][2] = 10.8409;
 x[125][3] = 2.7827;
 x[125][4] = -4.0974;
 y[125][0] = a;

 x[126][0] = 0.001;
 x[126][1] = -1.3;
 x[126][2] = 10.2678;
 x[126][3] = -2.953;
 x[126][4] = -5.8638;
 y[126][0] = a;

 x[127][0] = 0.001;
 x[127][1] = 3.4246;
 x[127][2] = -0.14693;
 x[127][3] = 0.80342;
 x[127][4] = 0.29136;
 y[127][0] = a;

 x[128][0] = 0.001;
 x[128][1] = 2.5503;
 x[128][2] = -4.9518;
 x[128][3] = 6.3729;
 x[128][4] = -0.41596;
 y[128][0] = a;

 x[129][0] = 0.001;
 x[129][1] = 1.5691;
 x[129][2] = 6.3465;
 x[129][3] = -0.1828;
 x[129][4] = -2.4099;
 y[129][0] = a;

 x[130][0] = 0.001;
 x[130][1] = 1.3087;
 x[130][2] = 4.9228;
 x[130][3] = 2.0013;
 x[130][4] = 0.22024;
 y[130][0] = a;

 x[131][0] = 0.001;
 x[131][1] = 5.1776;
 x[131][2] = 8.2316;
 x[131][3] = -3.2511;
 x[131][4] = -1.5694;
 y[131][0] = a;

 x[132][0] = 0.001;
 x[132][1] = 2.229;
 x[132][2] = 9.6325;
 x[132][3] = -3.1123;
 x[132][4] = -2.7164;
 y[132][0] = a;

 x[133][0] = 0.001;
 x[133][1] = 5.6272;
 x[133][2] = 10.0857;
 x[133][3] = -4.2931;
 x[133][4] = -3.8142;
 y[133][0] = a;

 x[134][0] = 0.001;
 x[134][1] = 1.2138;
 x[134][2] = 8.7986;
 x[134][3] = -2.1672;
 x[134][4] = -0.74182;
 y[134][0] = a;

 x[135][0] = 0.001;
 x[135][1] = 0.3798;
 x[135][2] = 0.7098;
 x[135][3] = 0.7572;
 x[135][4] = -0.4444;
 y[135][0] = a;

 x[136][0] = 0.001;
 x[136][1] = 0.5415;
 x[136][2] = 6.0319;
 x[136][3] = 1.6825;
 x[136][4] = -0.46122;
 y[136][0] = a;

 x[137][0] = 0.001;
 x[137][1] = 4.0524;
 x[137][2] = 5.6802;
 x[137][3] = -1.9693;
 x[137][4] = 0.026279;
 y[137][0] = a;

 x[138][0] = 0.001;
 x[138][1] = 4.7285;
 x[138][2] = 2.1065;
 x[138][3] = -0.28305;
 x[138][4] = 1.5625;
 y[138][0] = a;

 x[139][0] = 0.001;
 x[139][1] = 3.4359;
 x[139][2] = 0.66216;
 x[139][3] = 2.1041;
 x[139][4] = 1.8922;
 y[139][0] = a;

 x[140][0] = 0.001;
 x[140][1] = 0.86816;
 x[140][2] = 10.2429;
 x[140][3] = -1.4912;
 x[140][4] = -4.0082;
 y[140][0] = a;

 x[141][0] = 0.001;
 x[141][1] = 3.359;
 x[141][2] = 9.8022;
 x[141][3] = -3.8209;
 x[141][4] = -3.7133;
 y[141][0] = a;

 x[142][0] = 0.001;
 x[142][1] = 3.6702;
 x[142][2] = 2.9942;
 x[142][3] = 0.85141;
 x[142][4] = 0.30688;
 y[142][0] = a;

 x[143][0] = 0.001;
 x[143][1] = 1.3349;
 x[143][2] = 6.1189;
 x[143][3] = 0.46497;
 x[143][4] = 0.49826;
 y[143][0] = a;

 x[144][0] = 0.001;
 x[144][1] = 3.1887;
 x[144][2] = -3.4143;
 x[144][3] = 2.7742;
 x[144][4] = -0.2026;
 y[144][0] = a;

 x[145][0] = 0.001;
 x[145][1] = 2.4527;
 x[145][2] = 2.9653;
 x[145][3] = 0.20021;
 x[145][4] = -0.056479;
 y[145][0] = a;

 x[146][0] = 0.001;
 x[146][1] = 3.9121;
 x[146][2] = 2.9735;
 x[146][3] = 0.92852;
 x[146][4] = 0.60558;
 y[146][0] = a;

 x[147][0] = 0.001;
 x[147][1] = 3.9364;
 x[147][2] = 10.5885;
 x[147][3] = -3.725;
 x[147][4] = -4.3133;
 y[147][0] = a;

 x[148][0] = 0.001;
 x[148][1] = 3.9414;
 x[148][2] = -3.2902;
 x[148][3] = 3.1674;
 x[148][4] = 1.0866;
 y[148][0] = a;

 x[149][0] = 0.001;
 x[149][1] = 3.6922;
 x[149][2] = -3.9585;
 x[149][3] = 4.3439;
 x[149][4] = 1.3517;
 y[149][0] = a;

 x[150][0] = 0.001;
 x[150][1] = 5.681;
 x[150][2] = 7.795;
 x[150][3] = -2.6848;
 x[150][4] = -0.92544;
 y[150][0] = a;

 x[151][0] = 0.001;
 x[151][1] = 0.77124;
 x[151][2] = 9.0862;
 x[151][3] = -1.2281;
 x[151][4] = -1.4996;
 y[151][0] = a;

 x[152][0] = 0.001;
 x[152][1] = 3.5761;
 x[152][2] = 9.7753;
 x[152][3] = -3.9795;
 x[152][4] = -3.4638;
 y[152][0] = a;

 x[153][0] = 0.001;
 x[153][1] = 1.602;
 x[153][2] = 6.1251;
 x[153][3] = 0.52924;
 x[153][4] = 0.47886;
 y[153][0] = a;

 x[154][0] = 0.001;
 x[154][1] = 2.6682;
 x[154][2] = 10.216;
 x[154][3] = -3.4414;
 x[154][4] = -4.0069;
 y[154][0] = a;

 x[155][0] = 0.001;
 x[155][1] = 2.0007;
 x[155][2] = 1.8644;
 x[155][3] = 2.6491;
 x[155][4] = 0.47369;
 y[155][0] = a;

 x[156][0] = 0.001;
 x[156][1] = 0.64215;
 x[156][2] = 3.1287;
 x[156][3] = 4.2933;
 x[156][4] = 0.64696;
 y[156][0] = a;

 x[157][0] = 0.001;
 x[157][1] = 4.3848;
 x[157][2] = -3.0729;
 x[157][3] = 3.0423;
 x[157][4] = 1.2741;
 y[157][0] = a;

 x[158][0] = 0.001;
 x[158][1] = 0.77445;
 x[158][2] = 9.0552;
 x[158][3] = -2.4089;
 x[158][4] = -1.3884;
 y[158][0] = a;

 x[159][0] = 0.001;
 x[159][1] = 0.96574;
 x[159][2] = 8.393;
 x[159][3] = -1.361;
 x[159][4] = -1.4659;
 y[159][0] = a;

 x[160][0] = 0.001;
 x[160][1] = 3.0948;
 x[160][2] = 8.7324;
 x[160][3] = -2.9007;
 x[160][4] = -0.96682;
 y[160][0] = a;

 x[161][0] = 0.001;
 x[161][1] = 4.9362;
 x[161][2] = 7.6046;
 x[161][3] = -2.3429;
 x[161][4] = -0.85302;
 y[161][0] = a;

 x[162][0] = 0.001;
 x[162][1] = -1.9458;
 x[162][2] = 11.2217;
 x[162][3] = 1.9079;
 x[162][4] = -3.4405;
 y[162][0] = a;

 x[163][0] = 0.001;
 x[163][1] = 5.7403;
 x[163][2] = -0.44284;
 x[163][3] = 0.38015;
 x[163][4] = 1.3763;
 y[163][0] = a;

 x[164][0] = 0.001;
 x[164][1] = -2.6989;
 x[164][2] = 12.1984;
 x[164][3] = 0.67661;
 x[164][4] = -8.5482;
 y[164][0] = a;

 x[165][0] = 0.001;
 x[165][1] = 1.1472;
 x[165][2] = 3.5985;
 x[165][3] = 1.9387;
 x[165][4] = -0.43406;
 y[165][0] = a;

 x[166][0] = 0.001;
 x[166][1] = 2.9742;
 x[166][2] = 8.96;
 x[166][3] = -2.9024;
 x[166][4] = -1.0379;
 y[166][0] = a;

 x[167][0] = 0.001;
 x[167][1] = 4.5707;
 x[167][2] = 7.2094;
 x[167][3] = -3.2794;
 x[167][4] = -1.4944;
 y[167][0] = a;

 x[168][0] = 0.001;
 x[168][1] = 0.1848;
 x[168][2] = 6.5079;
 x[168][3] = 2.0133;
 x[168][4] = -0.87242;
 y[168][0] = a;

 x[169][0] = 0.001;
 x[169][1] = 0.87256;
 x[169][2] = 9.2931;
 x[169][3] = -0.7843;
 x[169][4] = -2.1978;
 y[169][0] = a;

 x[170][0] = 0.001;
 x[170][1] = 0.39559;
 x[170][2] = 6.8866;
 x[170][3] = 1.0588;
 x[170][4] = -0.67587;
 y[170][0] = a;

 x[171][0] = 0.001;
 x[171][1] = 3.8384;
 x[171][2] = 6.1851;
 x[171][3] = -2.0439;
 x[171][4] = -0.033204;
 y[171][0] = a;

 x[172][0] = 0.001;
 x[172][1] = 2.8209;
 x[172][2] = 7.3108;
 x[172][3] = -0.81857;
 x[172][4] = -1.8784;
 y[172][0] = a;

 x[173][0] = 0.001;
 x[173][1] = 2.5817;
 x[173][2] = 9.7546;
 x[173][3] = -3.1749;
 x[173][4] = -2.9957;
 y[173][0] = a;

 x[174][0] = 0.001;
 x[174][1] = 3.8213;
 x[174][2] = 0.23175;
 x[174][3] = 2.0133;
 x[174][4] = 2.0564;
 y[174][0] = a;

 x[175][0] = 0.001;
 x[175][1] = 0.3798;
 x[175][2] = 0.7098;
 x[175][3] = 0.7572;
 x[175][4] = -0.4444;
 y[175][0] = a;

 x[176][0] = 0.001;
 x[176][1] = 3.4893;
 x[176][2] = 6.69;
 x[176][3] = -1.2042;
 x[176][4] = -0.38751;
 y[176][0] = a;

 x[177][0] = 0.001;
 x[177][1] = -1.7781;
 x[177][2] = 0.8546;
 x[177][3] = 7.1303;
 x[177][4] = 0.027572;
 y[177][0] = a;

 x[178][0] = 0.001;
 x[178][1] = 2.0962;
 x[178][2] = 2.4769;
 x[178][3] = 1.9379;
 x[178][4] = -0.040962;
 y[178][0] = a;

 x[179][0] = 0.001;
 x[179][1] = 0.94732;
 x[179][2] = -0.57113;
 x[179][3] = 7.1903;
 x[179][4] = -0.67587;
 y[179][0] = a;

 x[180][0] = 0.001;
 x[180][1] = 2.8261;
 x[180][2] = 9.4007;
 x[180][3] = -3.3034;
 x[180][4] = -1.0509;
 y[180][0] = a;

 x[181][0] = 0.001;
 x[181][1] = 0.007125;
 x[181][2] = 8.3661;
 x[181][3] = 0.50781;
 x[181][4] = -3.8155;
 y[181][0] = a;

 x[182][0] = 0.001;
 x[182][1] = 0.96788;
 x[182][2] = 7.1907;
 x[182][3] = 1.2798;
 x[182][4] = -2.4565;
 y[182][0] = a;

 x[183][0] = 0.001;
 x[183][1] = 4.7432;
 x[183][2] = 2.1086;
 x[183][3] = 0.1368;
 x[183][4] = 1.6543;
 y[183][0] = a;

 x[184][0] = 0.001;
 x[184][1] = 3.6575;
 x[184][2] = 7.2797;
 x[184][3] = -2.2692;
 x[184][4] = -1.144;
 y[184][0] = a;

 x[185][0] = 0.001;
 x[185][1] = 3.8832;
 x[185][2] = 6.4023;
 x[185][3] = -2.432;
 x[185][4] = -0.98363;
 y[185][0] = a;

 x[186][0] = 0.001;
 x[186][1] = 3.4776;
 x[186][2] = 8.811;
 x[186][3] = -3.1886;
 x[186][4] = -0.92285;
 y[186][0] = a;

 x[187][0] = 0.001;
 x[187][1] = 1.1315;
 x[187][2] = 7.9212;
 x[187][3] = 1.093;
 x[187][4] = -2.8444;
 y[187][0] = a;

 x[188][0] = 0.001;
 x[188][1] = 2.8237;
 x[188][2] = 2.8597;
 x[188][3] = 0.19678;
 x[188][4] = 0.57196;
 y[188][0] = a;

 x[189][0] = 0.001;
 x[189][1] = 1.9321;
 x[189][2] = 6.0423;
 x[189][3] = 0.26019;
 x[189][4] = -2.053;
 y[189][0] = a;

 x[190][0] = 0.001;
 x[190][1] = 3.0632;
 x[190][2] = -3.3315;
 x[190][3] = 5.1305;
 x[190][4] = 0.8267;
 y[190][0] = a;

 x[191][0] = 0.001;
 x[191][1] = -1.8411;
 x[191][2] = 10.8306;
 x[191][3] = 2.769;
 x[191][4] = -3.0901;
 y[191][0] = a;

 x[192][0] = 0.001;
 x[192][1] = 2.8084;
 x[192][2] = 11.3045;
 x[192][3] = -3.3394;
 x[192][4] = -4.4194;
 y[192][0] = a;

 x[193][0] = 0.001;
 x[193][1] = 2.5698;
 x[193][2] = -4.4076;
 x[193][3] = 5.9856;
 x[193][4] = 0.078002;
 y[193][0] = a;

 x[194][0] = 0.001;
 x[194][1] = -0.12624;
 x[194][2] = 10.3216;
 x[194][3] = -3.7121;
 x[194][4] = -6.1185;
 y[194][0] = a;

 x[195][0] = 0.001;
 x[195][1] = 3.3756;
 x[195][2] = -4.0951;
 x[195][3] = 4.367;
 x[195][4] = 1.0698;
 y[195][0] = a;

 x[196][0] = 0.001;
 x[196][1] = -0.048008;
 x[196][2] = -1.6037;
 x[196][3] = 8.4756;
 x[196][4] = 0.75558;
 y[196][0] = a;

 x[197][0] = 0.001;
 x[197][1] = 0.5706;
 x[197][2] = -0.0248;
 x[197][3] = 1.2421;
 x[197][4] = -0.5621;
 y[197][0] = a;

 x[198][0] = 0.001;
 x[198][1] = 0.88444;
 x[198][2] = 6.5906;
 x[198][3] = 0.55837;
 x[198][4] = -0.44182;
 y[198][0] = a;

 x[199][0] = 0.001;
 x[199][1] = 3.8644;
 x[199][2] = 3.7061;
 x[199][3] = 0.70403;
 x[199][4] = 0.35214;
 y[199][0] = a;

 x[200][0] = 0.001;
 x[200][1] = 1.2999;
 x[200][2] = 2.5762;
 x[200][3] = 2.0107;
 x[200][4] = -0.18967;
 y[200][0] = a;

 x[201][0] = 0.001;
 x[201][1] = 2.0051;
 x[201][2] = -6.8638;
 x[201][3] = 8.132;
 x[201][4] = -0.2401;
 y[201][0] = a;

 x[202][0] = 0.001;
 x[202][1] = 4.9294;
 x[202][2] = 0.27727;
 x[202][3] = 0.20792;
 x[202][4] = 0.33662;
 y[202][0] = a;

 x[203][0] = 0.001;
 x[203][1] = 2.8297;
 x[203][2] = 6.3485;
 x[203][3] = -0.73546;
 x[203][4] = -0.58665;
 y[203][0] = a;

 x[204][0] = 0.001;
 x[204][1] = 2.565;
 x[204][2] = 8.633;
 x[204][3] = -2.9941;
 x[204][4] = -1.3082;
 y[204][0] = a;

 x[205][0] = 0.001;
 x[205][1] = 2.093;
 x[205][2] = 8.3061;
 x[205][3] = 0.022844;
 x[205][4] = -3.2724;
 y[205][0] = a;

 x[206][0] = 0.001;
 x[206][1] = 4.6014;
 x[206][2] = 5.6264;
 x[206][3] = -2.1235;
 x[206][4] = 0.19309;
 y[206][0] = a;

 x[207][0] = 0.001;
 x[207][1] = 5.0617;
 x[207][2] = -0.35799;
 x[207][3] = 0.44698;
 x[207][4] = 0.99868;
 y[207][0] = a;

 x[208][0] = 0.001;
 x[208][1] = -0.2951;
 x[208][2] = 9.0489;
 x[208][3] = -0.52725;
 x[208][4] = -2.0789;
 y[208][0] = a;

 x[209][0] = 0.001;
 x[209][1] = 3.577;
 x[209][2] = 2.4004;
 x[209][3] = 1.8908;
 x[209][4] = 0.73231;
 y[209][0] = a;

 x[210][0] = 0.001;
 x[210][1] = 3.9433;
 x[210][2] = 2.5017;
 x[210][3] = 1.5215;
 x[210][4] = 0.903;
 y[210][0] = a;

 x[211][0] = 0.001;
 x[211][1] = 2.6648;
 x[211][2] = 10.754;
 x[211][3] = -3.3994;
 x[211][4] = -4.1685;
 y[211][0] = a;

 x[212][0] = 0.001;
 x[212][1] = 5.9374;
 x[212][2] = 6.1664;
 x[212][3] = -2.5905;
 x[212][4] = -0.36553;
 y[212][0] = a;

 x[213][0] = 0.001;
 x[213][1] = 2.0153;
 x[213][2] = 1.8479;
 x[213][3] = 3.1375;
 x[213][4] = 0.42843;
 y[213][0] = a;

 x[214][0] = 0.001;
 x[214][1] = 5.8782;
 x[214][2] = 5.9409;
 x[214][3] = -2.8544;
 x[214][4] = -0.60863;
 y[214][0] = a;

 x[215][0] = 0.001;
 x[215][1] = -2.3983;
 x[215][2] = 12.606;
 x[215][3] = 2.9464;
 x[215][4] = -5.7888;
 y[215][0] = a;

 x[216][0] = 0.001;
 x[216][1] = 1.762;
 x[216][2] = 4.3682;
 x[216][3] = 2.1384;
 x[216][4] = 0.75429;
 y[216][0] = a;

 x[217][0] = 0.001;
 x[217][1] = 4.2406;
 x[217][2] = -2.4852;
 x[217][3] = 1.608;
 x[217][4] = 0.7155;
 y[217][0] = a;

 x[218][0] = 0.001;
 x[218][1] = 3.4669;
 x[218][2] = 6.87;
 x[218][3] = -1.0568;
 x[218][4] = -0.73147;
 y[218][0] = a;

 x[219][0] = 0.001;
 x[219][1] = 3.1896;
 x[219][2] = 5.7526;
 x[219][3] = -0.18537;
 x[219][4] = -0.30087;
 y[219][0] = a;

 x[220][0] = 0.001;
 x[220][1] = 0.81356;
 x[220][2] = 9.1566;
 x[220][3] = -2.1492;
 x[220][4] = -4.1814;
 y[220][0] = a;

 x[221][0] = 0.001;
 x[221][1] = 0.52855;
 x[221][2] = 0.96427;
 x[221][3] = 4.0243;
 x[221][4] = -1.0483;
 y[221][0] = a;

 x[222][0] = 0.001;
 x[222][1] = 2.1319;
 x[222][2] = -2.0403;
 x[222][3] = 2.5574;
 x[222][4] = -0.061652;
 y[222][0] = a;

 x[223][0] = 0.001;
 x[223][1] = 0.33111;
 x[223][2] = 4.5731;
 x[223][3] = 2.057;
 x[223][4] = -0.18967;
 y[223][0] = a;

 x[224][0] = 0.001;
 x[224][1] = 1.2746;
 x[224][2] = 8.8172;
 x[224][3] = -1.5323;
 x[224][4] = -1.7957;
 y[224][0] = a;

 x[225][0] = 0.001;
 x[225][1] = 2.2091;
 x[225][2] = 7.4556;
 x[225][3] = -1.3284;
 x[225][4] = -3.3021;
 y[225][0] = a;

 x[226][0] = 0.001;
 x[226][1] = 2.5328;
 x[226][2] = 7.528;
 x[226][3] = -0.41929;
 x[226][4] = -2.6478;
 y[226][0] = a;

 x[227][0] = 0.001;
 x[227][1] = 3.6244;
 x[227][2] = 1.4609;
 x[227][3] = 1.3501;
 x[227][4] = 1.9284;
 y[227][0] = a;

 x[228][0] = 0.001;
 x[228][1] = -1.3885;
 x[228][2] = 12.5026;
 x[228][3] = 0.69118;
 x[228][4] = -7.5487;
 y[228][0] = a;

 x[229][0] = 0.001;
 x[229][1] = 5.7227;
 x[229][2] = 5.8312;
 x[229][3] = -2.4097;
 x[229][4] = -0.24527;
 y[229][0] = a;

 x[230][0] = 0.001;
 x[230][1] = 3.3583;
 x[230][2] = 10.3567;
 x[230][3] = -3.7301;
 x[230][4] = -3.6991;
 y[230][0] = a;

 x[231][0] = 0.001;
 x[231][1] = 2.5227;
 x[231][2] = 2.2369;
 x[231][3] = 2.7236;
 x[231][4] = 0.79438;
 y[231][0] = a;

 x[232][0] = 0.001;
 x[232][1] = 0.045304;
 x[232][2] = 6.7334;
 x[232][3] = 1.0708;
 x[232][4] = -0.9332;
 y[232][0] = a;

 x[233][0] = 0.001;
 x[233][1] = 4.8278;
 x[233][2] = 7.7598;
 x[233][3] = -2.4491;
 x[233][4] = -1.2216;
 y[233][0] = a;

 x[234][0] = 0.001;
 x[234][1] = 1.9476;
 x[234][2] = -4.7738;
 x[234][3] = 8.527;
 x[234][4] = -1.8668;
 y[234][0] = a;

 x[235][0] = 0.001;
 x[235][1] = 2.7659;
 x[235][2] = 0.66216;
 x[235][3] = 4.1494;
 x[235][4] = -0.28406;
 y[235][0] = a;

 x[236][0] = 0.001;
 x[236][1] = -0.90784;
 x[236][2] = -7.9026;
 x[236][3] = 6.7807;
 x[236][4] = 0.34179;
 y[236][0] = b;

 x[237][0] = 0.001;
 x[237][1] = -2.0042;
 x[237][2] = -9.3676;
 x[237][3] = 9.3333;
 x[237][4] = -0.10303;
 y[237][0] = b;

 x[238][0] = 0.001;
 x[238][1] = -0.93587;
 x[238][2] = -5.1008;
 x[238][3] = 4.5367;
 x[238][4] = 1.3866;
 y[238][0] = b;

 x[239][0] = 0.001;
 x[239][1] = -0.40804;
 x[239][2] = 0.54214;
 x[239][3] = -0.52725;
 x[239][4] = 0.6586;
 y[239][0] = b;

 x[240][0] = 0.001;
 x[240][1] = -0.8172;
 x[240][2] = 3.3812;
 x[240][3] = -3.6684;
 x[240][4] = -3.456;
 y[240][0] = b;

 x[241][0] = 0.001;
 x[241][1] = -4.8392;
 x[241][2] = 6.6755;
 x[241][3] = -0.24278;
 x[241][4] = -6.5775;
 y[241][0] = b;

 x[242][0] = 0.001;
 x[242][1] = -1.2792;
 x[242][2] = 2.1376;
 x[242][3] = -0.47584;
 x[242][4] = -1.3974;
 y[242][0] = b;

 x[243][0] = 0.001;
 x[243][1] = -0.66008;
 x[243][2] = -3.226;
 x[243][3] = 3.8058;
 x[243][4] = 1.1836;
 y[243][0] = b;

 x[244][0] = 0.001;
 x[244][1] = -1.7713;
 x[244][2] = -10.7665;
 x[244][3] = 10.2184;
 x[244][4] = -1.0043;
 y[244][0] = b;

 x[245][0] = 0.001;
 x[245][1] = -3.0061;
 x[245][2] = -12.2377;
 x[245][3] = 11.9552;
 x[245][4] = -2.1603;
 y[245][0] = b;

 x[246][0] = 0.001;
 x[246][1] = -1.1022;
 x[246][2] = -5.8395;
 x[246][3] = 4.5641;
 x[246][4] = 0.68705;
 y[246][0] = b;

 x[247][0] = 0.001;
 x[247][1] = 0.11806;
 x[247][2] = 0.39108;
 x[247][3] = -0.98223;
 x[247][4] = 0.42843;
 y[247][0] = b;

 x[248][0] = 0.001;
 x[248][1] = 0.11686;
 x[248][2] = 3.735;
 x[248][3] = -4.4379;
 x[248][4] = -4.3741;
 y[248][0] = b;

 x[249][0] = 0.001;
 x[249][1] = -2.7264;
 x[249][2] = 3.9213;
 x[249][3] = -0.49212;
 x[249][4] = -3.6371;
 y[249][0] = b;

 x[250][0] = 0.001;
 x[250][1] = -1.2369;
 x[250][2] = -1.6906;
 x[250][3] = 2.518;
 x[250][4] = 0.51636;
 y[250][0] = b;

 x[251][0] = 0.001;
 x[251][1] = -1.8439;
 x[251][2] = -8.6475;
 x[251][3] = 7.6796;
 x[251][4] = -0.66682;
 y[251][0] = b;

 x[252][0] = 0.001;
 x[252][1] = -1.8554;
 x[252][2] = -9.6035;
 x[252][3] = 7.7764;
 x[252][4] = -0.97716;
 y[252][0] = b;

 x[253][0] = 0.001;
 x[253][1] = 0.16358;
 x[253][2] = -3.3584;
 x[253][3] = 1.3749;
 x[253][4] = 1.3569;
 y[253][0] = b;

 x[254][0] = 0.001;
 x[254][1] = 1.5077;
 x[254][2] = 1.9596;
 x[254][3] = -3.0584;
 x[254][4] = -0.12243;
 y[254][0] = b;

 x[255][0] = 0.001;
 x[255][1] = 0.67886;
 x[255][2] = 4.1199;
 x[255][3] = -4.569;
 x[255][4] = -4.1414;
 y[255][0] = b;

 x[256][0] = 0.001;
 x[256][1] = -3.9934;
 x[256][2] = 5.8333;
 x[256][3] = 0.54723;
 x[256][4] = -4.9379;
 y[256][0] = b;

 x[257][0] = 0.001;
 x[257][1] = -2.3898;
 x[257][2] = -0.78427;
 x[257][3] = 3.0141;
 x[257][4] = 0.76205;
 y[257][0] = b;

 x[258][0] = 0.001;
 x[258][1] = -1.7976;
 x[258][2] = -6.7686;
 x[258][3] = 6.6753;
 x[258][4] = 0.89912;
 y[258][0] = b;

 x[259][0] = 0.001;
 x[259][1] = -0.70867;
 x[259][2] = -5.5602;
 x[259][3] = 4.0483;
 x[259][4] = 0.903;
 y[259][0] = b;

 x[260][0] = 0.001;
 x[260][1] = 1.0194;
 x[260][2] = 1.1029;
 x[260][3] = -2.3;
 x[260][4] = 0.59395;
 y[260][0] = b;

 x[261][0] = 0.001;
 x[261][1] = 1.7875;
 x[261][2] = 4.78;
 x[261][3] = -5.1362;
 x[261][4] = -3.2362;
 y[261][0] = b;

 x[262][0] = 0.001;
 x[262][1] = 0.27331;
 x[262][2] = 4.8773;
 x[262][3] = -4.9194;
 x[262][4] = -5.8198;
 y[262][0] = b;

 x[263][0] = 0.001;
 x[263][1] = -5.1661;
 x[263][2] = 8.0433;
 x[263][3] = 0.044265;
 x[263][4] = -4.4983;
 y[263][0] = b;

 x[264][0] = 0.001;
 x[264][1] = -2.7028;
 x[264][2] = 1.6327;
 x[264][3] = 0.83598;
 x[264][4] = -0.091393;
 y[264][0] = b;

 x[265][0] = 0.001;
 x[265][1] = -1.4904;
 x[265][2] = -2.2183;
 x[265][3] = 1.6054;
 x[265][4] = 0.89394;
 y[265][0] = b;

 x[266][0] = 0.001;
 x[266][1] = -0.014902;
 x[266][2] = -1.0243;
 x[266][3] = -0.94024;
 x[266][4] = 0.64955;
 y[266][0] = b;

 x[267][0] = 0.001;
 x[267][1] = 0.88992;
 x[267][2] = 2.2638;
 x[267][3] = -3.1046;
 x[267][4] = -0.11855;
 y[267][0] = b;

 x[268][0] = 0.001;
 x[268][1] = 1.0637;
 x[268][2] = 3.6957;
 x[268][3] = -4.1594;
 x[268][4] = -1.9379;
 y[268][0] = b;

 x[269][0] = 0.001;
 x[269][1] = -0.8471;
 x[269][2] = 3.1329;
 x[269][3] = -3.0112;
 x[269][4] = -2.9388;
 y[269][0] = b;

 x[270][0] = 0.001;
 x[270][1] = -3.9594;
 x[270][2] = 4.0289;
 x[270][3] = -0.35845;
 x[270][4] = -3.8957;
 y[270][0] = b;

 x[271][0] = 0.001;
 x[271][1] = -5.8818;
 x[271][2] = 7.6584;
 x[271][3] = 0.5558;
 x[271][4] = -2.9155;
 y[271][0] = b;

 x[272][0] = 0.001;
 x[272][1] = -3.7747;
 x[272][2] = 2.5162;
 x[272][3] = 0.83341;
 x[272][4] = -0.30993;
 y[272][0] = b;

 x[273][0] = 0.001;
 x[273][1] = -2.4198;
 x[273][2] = -0.24418;
 x[273][3] = 0.70146;
 x[273][4] = 0.41809;
 y[273][0] = b;

 x[274][0] = 0.001;
 x[274][1] = -0.83535;
 x[274][2] = 0.80494;
 x[274][3] = -1.6411;
 x[274][4] = -0.19225;
 y[274][0] = b;

 x[275][0] = 0.001;
 x[275][1] = -0.30432;
 x[275][2] = 2.6528;
 x[275][3] = -2.7756;
 x[275][4] = -0.65647;
 y[275][0] = b;

 x[276][0] = 0.001;
 x[276][1] = -0.60254;
 x[276][2] = 1.7237;
 x[276][3] = -2.1501;
 x[276][4] = -0.77027;
 y[276][0] = b;

 x[277][0] = 0.001;
 x[277][1] = -2.1059;
 x[277][2] = 1.1815;
 x[277][3] = -0.53324;
 x[277][4] = -0.82716;
 y[277][0] = b;

 x[278][0] = 0.001;
 x[278][1] = -2.0441;
 x[278][2] = 1.2271;
 x[278][3] = 0.18564;
 x[278][4] = -1.091;
 y[278][0] = b;

 x[279][0] = 0.001;
 x[279][1] = -1.5621;
 x[279][2] = -2.2121;
 x[279][3] = 4.2591;
 x[279][4] = 0.27972;
 y[279][0] = b;

 x[280][0] = 0.001;
 x[280][1] = -3.2305;
 x[280][2] = -7.2135;
 x[280][3] = 11.6433;
 x[280][4] = -0.94613;
 y[280][0] = b;

 x[281][0] = 0.001;
 x[281][1] = -4.8426;
 x[281][2] = -4.9932;
 x[281][3] = 10.4052;
 x[281][4] = -0.53104;
 y[281][0] = b;

 x[282][0] = 0.001;
 x[282][1] = -2.3147;
 x[282][2] = 3.6668;
 x[282][3] = -0.6969;
 x[282][4] = -1.2474;
 y[282][0] = b;

 x[283][0] = 0.001;
 x[283][1] = -0.11716;
 x[283][2] = 0.60422;
 x[283][3] = -0.38587;
 x[283][4] = -0.059065;
 y[283][0] = b;

 x[284][0] = 0.001;
 x[284][1] = -2.0066;
 x[284][2] = -6.719;
 x[284][3] = 9.0162;
 x[284][4] = 0.099985;
 y[284][0] = b;

 x[285][0] = 0.001;
 x[285][1] = -3.6961;
 x[285][2] = -13.6779;
 x[285][3] = 17.5795;
 x[285][4] = -2.6181;
 y[285][0] = b;

 x[286][0] = 0.001;
 x[286][1] = -3.6012;
 x[286][2] = -6.5389;
 x[286][3] = 10.5234;
 x[286][4] = -0.48967;
 y[286][0] = b;

 x[287][0] = 0.001;
 x[287][1] = -2.6286;
 x[287][2] = 0.18002;
 x[287][3] = 1.7956;
 x[287][4] = 0.97282;
 y[287][0] = b;

 x[288][0] = 0.001;
 x[288][1] = -0.82601;
 x[288][2] = 2.9611;
 x[288][3] = -1.2864;
 x[288][4] = -1.4647;
 y[288][0] = b;

 x[289][0] = 0.001;
 x[289][1] = 0.31803;
 x[289][2] = -0.99326;
 x[289][3] = 1.0947;
 x[289][4] = 0.88619;
 y[289][0] = b;

 x[290][0] = 0.001;
 x[290][1] = -1.4454;
 x[290][2] = -8.4385;
 x[290][3] = 8.8483;
 x[290][4] = 0.96894;
 y[290][0] = b;

 x[291][0] = 0.001;
 x[291][1] = -3.1423;
 x[291][2] = -13.0365;
 x[291][3] = 15.6773;
 x[291][4] = -0.66165;
 y[291][0] = b;

 x[292][0] = 0.001;
 x[292][1] = -2.5373;
 x[292][2] = -6.959;
 x[292][3] = 8.8054;
 x[292][4] = 1.5289;
 y[292][0] = b;

 x[293][0] = 0.001;
 x[293][1] = -1.366;
 x[293][2] = 0.18416;
 x[293][3] = 0.90539;
 x[293][4] = 1.5806;
 y[293][0] = b;

 x[294][0] = 0.001;
 x[294][1] = -1.7064;
 x[294][2] = 3.3088;
 x[294][3] = -2.2829;
 x[294][4] = -2.1978;
 y[294][0] = b;

 x[295][0] = 0.001;
 x[295][1] = -0.41965;
 x[295][2] = 2.9094;
 x[295][3] = -1.7859;
 x[295][4] = -2.2069;
 y[295][0] = b;

 x[296][0] = 0.001;
 x[296][1] = 0.37637;
 x[296][2] = -0.82358;
 x[296][3] = 0.78543;
 x[296][4] = 0.74524;
 y[296][0] = b;

 x[297][0] = 0.001;
 x[297][1] = -0.55355;
 x[297][2] = -7.9233;
 x[297][3] = 6.7156;
 x[297][4] = 0.74394;
 y[297][0] = b;

 x[298][0] = 0.001;
 x[298][1] = -1.6001;
 x[298][2] = -9.5828;
 x[298][3] = 9.4044;
 x[298][4] = 0.081882;
 y[298][0] = b;

 x[299][0] = 0.001;
 x[299][1] = -0.37013;
 x[299][2] = -5.554;
 x[299][3] = 4.7749;
 x[299][4] = 1.547;
 y[299][0] = b;

 x[300][0] = 0.001;
 x[300][1] = 0.12126;
 x[300][2] = 0.22347;
 x[300][3] = -0.47327;
 x[300][4] = 0.97024;
 y[300][0] = b;

 x[301][0] = 0.001;
 x[301][1] = -0.27068;
 x[301][2] = 3.2674;
 x[301][3] = -3.5562;
 x[301][4] = -3.0888;
 y[301][0] = b;

 x[302][0] = 0.001;
 x[302][1] = -5.119;
 x[302][2] = 6.6486;
 x[302][3] = -0.049987;
 x[302][4] = -6.5206;
 y[302][0] = b;

 x[303][0] = 0.001;
 x[303][1] = -1.3946;
 x[303][2] = 2.3134;
 x[303][3] = -0.44499;
 x[303][4] = -1.4905;
 y[303][0] = b;

 x[304][0] = 0.001;
 x[304][1] = -0.69879;
 x[304][2] = -3.3771;
 x[304][3] = 4.1211;
 x[304][4] = 1.5043;
 y[304][0] = b;

 x[305][0] = 0.001;
 x[305][1] = -1.48;
 x[305][2] = -10.5244;
 x[305][3] = 9.9176;
 x[305][4] = -0.5026;
 y[305][0] = b;

 x[306][0] = 0.001;
 x[306][1] = -2.6649;
 x[306][2] = -12.813;
 x[306][3] = 12.6689;
 x[306][4] = -1.9082;
 y[306][0] = b;

 x[307][0] = 0.001;
 x[307][1] = -0.62684;
 x[307][2] = -6.301;
 x[307][3] = 4.7843;
 x[307][4] = 1.106;
 y[307][0] = b;

 x[308][0] = 0.001;
 x[308][1] = 0.518;
 x[308][2] = 0.25865;
 x[308][3] = -0.84085;
 x[308][4] = 0.96118;
 y[308][0] = b;

 x[309][0] = 0.001;
 x[309][1] = 0.64376;
 x[309][2] = 3.764;
 x[309][3] = -4.4738;
 x[309][4] = -4.0483;
 y[309][0] = b;

 x[310][0] = 0.001;
 x[310][1] = -2.9821;
 x[310][2] = 4.1986;
 x[310][3] = -0.5898;
 x[310][4] = -3.9642;
 y[310][0] = b;

 x[311][0] = 0.001;
 x[311][1] = -1.4628;
 x[311][2] = -1.5706;
 x[311][3] = 2.4357;
 x[311][4] = 0.49826;
 y[311][0] = b;

 x[312][0] = 0.001;
 x[312][1] = -1.7101;
 x[312][2] = -8.7903;
 x[312][3] = 7.9735;
 x[312][4] = -0.45475;
 y[312][0] = b;

 x[313][0] = 0.001;
 x[313][1] = -1.5572;
 x[313][2] = -9.8808;
 x[313][3] = 8.1088;
 x[313][4] = -1.0806;
 y[313][0] = b;

 x[314][0] = 0.001;
 x[314][1] = 0.74428;
 x[314][2] = -3.7723;
 x[314][3] = 1.6131;
 x[314][4] = 1.5754;
 y[314][0] = b;

 x[315][0] = 0.001;
 x[315][1] = 2.0177;
 x[315][2] = 1.7982;
 x[315][3] = -2.9581;
 x[315][4] = 0.2099;
 y[315][0] = b;

 x[316][0] = 0.001;
 x[316][1] = 1.164;
 x[316][2] = 3.913;
 x[316][3] = -4.5544;
 x[316][4] = -3.8672;
 y[316][0] = b;

 x[317][0] = 0.001;
 x[317][1] = -4.3667;
 x[317][2] = 6.0692;
 x[317][3] = 0.57208;
 x[317][4] = -5.4668;
 y[317][0] = b;

 x[318][0] = 0.001;
 x[318][1] = -2.5919;
 x[318][2] = -1.0553;
 x[318][3] = 3.8949;
 x[318][4] = 0.77757;
 y[318][0] = b;

 x[319][0] = 0.001;
 x[319][1] = -1.8046;
 x[319][2] = -6.8141;
 x[319][3] = 6.7019;
 x[319][4] = 1.1681;
 y[319][0] = b;

 x[320][0] = 0.001;
 x[320][1] = -0.71868;
 x[320][2] = -5.7154;
 x[320][3] = 3.8298;
 x[320][4] = 1.0233;
 y[320][0] = b;

 x[321][0] = 0.001;
 x[321][1] = 1.4378;
 x[321][2] = 0.66837;
 x[321][3] = -2.0267;
 x[321][4] = 1.0271;
 y[321][0] = b;

 x[322][0] = 0.001;
 x[322][1] = 2.1943;
 x[322][2] = 4.5503;
 x[322][3] = -4.976;
 x[322][4] = -2.7254;
 y[322][0] = b;

 x[323][0] = 0.001;
 x[323][1] = 0.7376;
 x[323][2] = 4.8525;
 x[323][3] = -4.7986;
 x[323][4] = -5.6659;
 y[323][0] = b;

 x[324][0] = 0.001;
 x[324][1] = -5.637;
 x[324][2] = 8.1261;
 x[324][3] = 0.13081;
 x[324][4] = -5.0142;
 y[324][0] = b;

 x[325][0] = 0.001;
 x[325][1] = -3.0193;
 x[325][2] = 1.7775;
 x[325][3] = 0.73745;
 x[325][4] = -0.45346;
 y[325][0] = b;

 x[326][0] = 0.001;
 x[326][1] = -1.6706;
 x[326][2] = -2.09;
 x[326][3] = 1.584;
 x[326][4] = 0.71162;
 y[326][0] = b;

 x[327][0] = 0.001;
 x[327][1] = -0.1269;
 x[327][2] = -1.1505;
 x[327][3] = -0.95138;
 x[327][4] = 0.57843;
 y[327][0] = b;

 x[328][0] = 0.001;
 x[328][1] = 1.2198;
 x[328][2] = 2.0982;
 x[328][3] = -3.1954;
 x[328][4] = 0.12843;
 y[328][0] = b;

 x[329][0] = 0.001;
 x[329][1] = 1.4501;
 x[329][2] = 3.6067;
 x[329][3] = -4.0557;
 x[329][4] = -1.5966;
 y[329][0] = b;

 x[330][0] = 0.001;
 x[330][1] = -0.40857;
 x[330][2] = 3.0977;
 x[330][3] = -2.9607;
 x[330][4] = -2.6892;
 y[330][0] = b;

 x[331][0] = 0.001;
 x[331][1] = -3.8952;
 x[331][2] = 3.8157;
 x[331][3] = -0.31304;
 x[331][4] = -3.8194;
 y[331][0] = b;

 x[332][0] = 0.001;
 x[332][1] = -6.3679;
 x[332][2] = 8.0102;
 x[332][3] = 0.4247;
 x[332][4] = -3.2207;
 y[332][0] = b;

 x[333][0] = 0.001;
 x[333][1] = -4.1429;
 x[333][2] = 2.7749;
 x[333][3] = 0.68261;
 x[333][4] = -0.71984;
 y[333][0] = b;

 x[334][0] = 0.001;
 x[334][1] = -2.6864;
 x[334][2] = -0.097265;
 x[334][3] = 0.61663;
 x[334][4] = 0.061192;
 y[334][0] = b;

 x[335][0] = 0.001;
 x[335][1] = -1.0555;
 x[335][2] = 0.79459;
 x[335][3] = -1.6968;
 x[335][4] = -0.46768;
 y[335][0] = b;

 x[336][0] = 0.001;
 x[336][1] = -0.29858;
 x[336][2] = 2.4769;
 x[336][3] = -2.9512;
 x[336][4] = -0.66165;
 y[336][0] = b;

 x[337][0] = 0.001;
 x[337][1] = -0.49948;
 x[337][2] = 1.7734;
 x[337][3] = -2.2469;
 x[337][4] = -0.68104;
 y[337][0] = b;

 x[338][0] = 0.001;
 x[338][1] = -1.9881;
 x[338][2] = 0.99945;
 x[338][3] = -0.28562;
 x[338][4] = -0.70044;
 y[338][0] = b;

 x[339][0] = 0.001;
 x[339][1] = -1.9389;
 x[339][2] = 1.5706;
 x[339][3] = 0.045979;
 x[339][4] = -1.122;
 y[339][0] = b;

 x[340][0] = 0.001;
 x[340][1] = -1.4375;
 x[340][2] = -1.8624;
 x[340][3] = 4.026;
 x[340][4] = 0.55127;
 y[340][0] = b;

 x[341][0] = 0.001;
 x[341][1] = -3.1875;
 x[341][2] = -7.5756;
 x[341][3] = 11.8678;
 x[341][4] = -0.57889;
 y[341][0] = b;

 x[342][0] = 0.001;
 x[342][1] = -4.6765;
 x[342][2] = -5.6636;
 x[342][3] = 10.969;
 x[342][4] = -0.33449;
 y[342][0] = b;

 x[343][0] = 0.001;
 x[343][1] = -2.0285;
 x[343][2] = 3.8468;
 x[343][3] = -0.63435;
 x[343][4] = -1.175;
 y[343][0] = b;

 x[344][0] = 0.001;
 x[344][1] = 0.26637;
 x[344][2] = 0.73252;
 x[344][3] = -0.67891;
 x[344][4] = 0.03533;
 y[344][0] = b;

 x[345][0] = 0.001;
 x[345][1] = -1.7589;
 x[345][2] = -6.4624;
 x[345][3] = 8.4773;
 x[345][4] = 0.31981;
 y[345][0] = b;

 x[346][0] = 0.001;
 x[346][1] = -3.5985;
 x[346][2] = -13.6593;
 x[346][3] = 17.6052;
 x[346][4] = -2.4927;
 y[346][0] = b;

 x[347][0] = 0.001;
 x[347][1] = -3.3582;
 x[347][2] = -7.2404;
 x[347][3] = 11.4419;
 x[347][4] = -0.57113;
 y[347][0] = b;

 x[348][0] = 0.001;
 x[348][1] = -2.3629;
 x[348][2] = -0.10554;
 x[348][3] = 1.9336;
 x[348][4] = 1.1358;
 y[348][0] = b;

 x[349][0] = 0.001;
 x[349][1] = -2.1802;
 x[349][2] = 3.3791;
 x[349][3] = -1.2256;
 x[349][4] = -2.6621;
 y[349][0] = b;

 x[350][0] = 0.001;
 x[350][1] = -0.40951;
 x[350][2] = -0.15521;
 x[350][3] = 0.060545;
 x[350][4] = -0.088807;
 y[350][0] = b;

 x[351][0] = 0.001;
 x[351][1] = -2.2918;
 x[351][2] = -7.257;
 x[351][3] = 7.9597;
 x[351][4] = 0.9211;
 y[351][0] = b;

 x[352][0] = 0.001;
 x[352][1] = -4.0214;
 x[352][2] = -12.8006;
 x[352][3] = 15.6199;
 x[352][4] = -0.95647;
 y[352][0] = b;

 x[353][0] = 0.001;
 x[353][1] = -3.3884;
 x[353][2] = -8.215;
 x[353][3] = 10.3315;
 x[353][4] = 0.98187;
 y[353][0] = b;

 x[354][0] = 0.001;
 x[354][1] = -2.0046;
 x[354][2] = -0.49457;
 x[354][3] = 1.333;
 x[354][4] = 1.6543;
 y[354][0] = b;

 x[355][0] = 0.001;
 x[355][1] = -1.7063;
 x[355][2] = 2.7956;
 x[355][3] = -2.378;
 x[355][4] = -2.3491;
 y[355][0] = b;

 x[356][0] = 0.001;
 x[356][1] = -1.6386;
 x[356][2] = 3.3584;
 x[356][3] = -1.7302;
 x[356][4] = -3.5646;
 y[356][0] = b;

 x[357][0] = 0.001;
 x[357][1] = -0.41645;
 x[357][2] = 0.32487;
 x[357][3] = -0.33617;
 x[357][4] = -0.36036;
 y[357][0] = b;

 x[358][0] = 0.001;
 x[358][1] = -1.5877;
 x[358][2] = -6.6072;
 x[358][3] = 5.8022;
 x[358][4] = 0.31593;
 y[358][0] = b;

 x[359][0] = 0.001;
 x[359][1] = -2.5961;
 x[359][2] = -9.349;
 x[359][3] = 9.7942;
 x[359][4] = -0.28018;
 y[359][0] = b;

 x[360][0] = 0.001;
 x[360][1] = -1.5228;
 x[360][2] = -6.4789;
 x[360][3] = 5.7568;
 x[360][4] = 0.87325;
 y[360][0] = b;

 x[361][0] = 0.001;
 x[361][1] = -0.53072;
 x[361][2] = -0.097265;
 x[361][3] = -0.21793;
 x[361][4] = 1.0426;
 y[361][0] = b;

 x[362][0] = 0.001;
 x[362][1] = -0.49081;
 x[362][2] = 2.8452;
 x[362][3] = -3.6436;
 x[362][4] = -3.1004;
 y[362][0] = b;

 x[363][0] = 0.001;
 x[363][1] = -6.5773;
 x[363][2] = 6.8017;
 x[363][3] = 0.85483;
 x[363][4] = -7.5344;
 y[363][0] = b;

 x[364][0] = 0.001;
 x[364][1] = -2.4621;
 x[364][2] = 2.7645;
 x[364][3] = -0.62578;
 x[364][4] = -2.8573;
 y[364][0] = b;

 x[365][0] = 0.001;
 x[365][1] = -1.3995;
 x[365][2] = -1.9162;
 x[365][3] = 2.5154;
 x[365][4] = 0.59912;
 y[365][0] = b;

 x[366][0] = 0.001;
 x[366][1] = -2.3221;
 x[366][2] = -9.3304;
 x[366][3] = 9.233;
 x[366][4] = -0.79871;
 y[366][0] = b;

 x[367][0] = 0.001;
 x[367][1] = -3.73;
 x[367][2] = -12.9723;
 x[367][3] = 12.9817;
 x[367][4] = -2.684;
 y[367][0] = b;

 x[368][0] = 0.001;
 x[368][1] = -1.6988;
 x[368][2] = -7.1163;
 x[368][3] = 5.7902;
 x[368][4] = 0.16723;
 y[368][0] = b;

 x[369][0] = 0.001;
 x[369][1] = -0.26654;
 x[369][2] = -0.64562;
 x[369][3] = -0.42014;
 x[369][4] = 0.89136;
 y[369][0] = b;

 x[370][0] = 0.001;
 x[370][1] = 0.33325;
 x[370][2] = 3.3108;
 x[370][3] = -4.5081;
 x[370][4] = -4.012;
 y[370][0] = b;

 x[371][0] = 0.001;
 x[371][1] = -4.2091;
 x[371][2] = 4.7283;
 x[371][3] = -0.49126;
 x[371][4] = -5.2159;
 y[371][0] = b;

 x[372][0] = 0.001;
 x[372][1] = -2.3142;
 x[372][2] = -0.68494;
 x[372][3] = 1.9833;
 x[372][4] = -0.44829;
 y[372][0] = b;

 x[373][0] = 0.001;
 x[373][1] = -2.4835;
 x[373][2] = -7.4494;
 x[373][3] = 6.8964;
 x[373][4] = -0.64484;
 y[373][0] = b;

 x[374][0] = 0.001;
 x[374][1] = -2.7611;
 x[374][2] = -10.5099;
 x[374][3] = 9.0239;
 x[374][4] = -1.9547;
 y[374][0] = b;

 x[375][0] = 0.001;
 x[375][1] = -0.36025;
 x[375][2] = -4.449;
 x[375][3] = 2.1067;
 x[375][4] = 0.94308;
 y[375][0] = b;

 x[376][0] = 0.001;
 x[376][1] = 1.0117;
 x[376][2] = 0.9022;
 x[376][3] = -2.3506;
 x[376][4] = 0.42714;
 y[376][0] = b;

 x[377][0] = 0.001;
 x[377][1] = 0.96708;
 x[377][2] = 3.8426;
 x[377][3] = -4.9314;
 x[377][4] = -4.1323;
 y[377][0] = b;

 x[378][0] = 0.001;
 x[378][1] = -5.2049;
 x[378][2] = 7.259;
 x[378][3] = 0.070827;
 x[378][4] = -7.3004;
 y[378][0] = b;

 x[379][0] = 0.001;
 x[379][1] = -3.3203;
 x[379][2] = -0.02691;
 x[379][3] = 2.9618;
 x[379][4] = -0.44958;
 y[379][0] = b;

 x[380][0] = 0.001;
 x[380][1] = -2.565;
 x[380][2] = -5.7899;
 x[380][3] = 6.0122;
 x[380][4] = 0.046968;
 y[380][0] = b;

 x[381][0] = 0.001;
 x[381][1] = -1.5951;
 x[381][2] = -6.572;
 x[381][3] = 4.7689;
 x[381][4] = -0.94354;
 y[381][0] = b;

 x[382][0] = 0.001;
 x[382][1] = 0.7049;
 x[382][2] = 0.17174;
 x[382][3] = -1.7859;
 x[382][4] = 0.36119;
 y[382][0] = b;

 x[383][0] = 0.001;
 x[383][1] = 1.7331;
 x[383][2] = 3.9544;
 x[383][3] = -4.7412;
 x[383][4] = -2.5017;
 y[383][0] = b;

 x[384][0] = 0.001;
 x[384][1] = 0.6818;
 x[384][2] = 4.8504;
 x[384][3] = -5.2133;
 x[384][4] = -6.1043;
 y[384][0] = b;

 x[385][0] = 0.001;
 x[385][1] = -6.3364;
 x[385][2] = 9.2848;
 x[385][3] = 0.014275;
 x[385][4] = -6.7844;
 y[385][0] = b;

 x[386][0] = 0.001;
 x[386][1] = -3.8053;
 x[386][2] = 2.4273;
 x[386][3] = 0.6809;
 x[386][4] = -1.0871;
 y[386][0] = b;

 x[387][0] = 0.001;
 x[387][1] = -2.1979;
 x[387][2] = -2.1252;
 x[387][3] = 1.7151;
 x[387][4] = 0.45171;
 y[387][0] = b;

 x[388][0] = 0.001;
 x[388][1] = -0.87874;
 x[388][2] = -2.2121;
 x[388][3] = -0.051701;
 x[388][4] = 0.099985;
 y[388][0] = b;

 x[389][0] = 0.001;
 x[389][1] = 0.74067;
 x[389][2] = 1.7299;
 x[389][3] = -3.1963;
 x[389][4] = -0.1457;
 y[389][0] = b;

 x[390][0] = 0.001;
 x[390][1] = 0.98296;
 x[390][2] = 3.4226;
 x[390][3] = -3.9692;
 x[390][4] = -1.7116;
 y[390][0] = b;

 x[391][0] = 0.001;
 x[391][1] = -0.3489;
 x[391][2] = 3.1929;
 x[391][3] = -3.4054;
 x[391][4] = -3.1832;
 y[391][0] = b;

 x[392][0] = 0.001;
 x[392][1] = -3.8552;
 x[392][2] = 3.5219;
 x[392][3] = -0.38415;
 x[392][4] = -3.8608;
 y[392][0] = b;

 x[393][0] = 0.001;
 x[393][1] = -6.9599;
 x[393][2] = 8.9931;
 x[393][3] = 0.2182;
 x[393][4] = -4.572;
 y[393][0] = b;

 x[394][0] = 0.001;
 x[394][1] = -4.7462;
 x[394][2] = 3.1205;
 x[394][3] = 1.075;
 x[394][4] = -1.2966;
 y[394][0] = b;

 x[395][0] = 0.001;
 x[395][1] = -3.2051;
 x[395][2] = -0.14279;
 x[395][3] = 0.97565;
 x[395][4] = 0.045675;
 y[395][0] = b;

 x[396][0] = 0.001;
 x[396][1] = -1.7549;
 x[396][2] = -0.080711;
 x[396][3] = -0.75774;
 x[396][4] = -0.3707;
 y[396][0] = b;

 x[397][0] = 0.001;
 x[397][1] = -0.59587;
 x[397][2] = 2.4811;
 x[397][3] = -2.8673;
 x[397][4] = -0.89828;
 y[397][0] = b;

 x[398][0] = 0.001;
 x[398][1] = -0.89542;
 x[398][2] = 2.0279;
 x[398][3] = -2.3652;
 x[398][4] = -1.2746;
 y[398][0] = b;

 x[399][0] = 0.001;
 x[399][1] = -2.0754;
 x[399][2] = 1.2767;
 x[399][3] = -0.64206;
 x[399][4] = -1.2642;
 y[399][0] = b;

 x[400][0] = 0.001;
 x[400][1] = -3.2778;
 x[400][2] = 1.8023;
 x[400][3] = 0.1805;
 x[400][4] = -2.3931;
 y[400][0] = b;

 x[401][0] = 0.001;
 x[401][1] = -2.2183;
 x[401][2] = -1.254;
 x[401][3] = 2.9986;
 x[401][4] = 0.36378;
 y[401][0] = b;

 x[402][0] = 0.001;
 x[402][1] = -3.5895;
 x[402][2] = -6.572;
 x[402][3] = 10.5251;
 x[402][4] = -0.16381;
 y[402][0] = b;

 x[403][0] = 0.001;
 x[403][1] = -5.0477;
 x[403][2] = -5.8023;
 x[403][3] = 11.244;
 x[403][4] = -0.3901;
 y[403][0] = b;

 x[404][0] = 0.001;
 x[404][1] = -3.5741;
 x[404][2] = 3.944;
 x[404][3] = -0.07912;
 x[404][4] = -2.1203;
 y[404][0] = b;

 x[405][0] = 0.001;
 x[405][1] = -0.7351;
 x[405][2] = 1.7361;
 x[405][3] = -1.4938;
 x[405][4] = -1.1582;
 y[405][0] = b;

 x[406][0] = 0.001;
 x[406][1] = -2.2617;
 x[406][2] = -4.7428;
 x[406][3] = 6.3489;
 x[406][4] = 0.11162;
 y[406][0] = b;

 x[407][0] = 0.001;
 x[407][1] = -4.244;
 x[407][2] = -13.0634;
 x[407][3] = 17.1116;
 x[407][4] = -2.8017;
 y[407][0] = b;

 x[408][0] = 0.001;
 x[408][1] = -4.0218;
 x[408][2] = -8.304;
 x[408][3] = 12.555;
 x[408][4] = -1.5099;
 y[408][0] = b;

 x[409][0] = 0.001;
 x[409][1] = -3.0201;
 x[409][2] = -0.67253;
 x[409][3] = 2.7056;
 x[409][4] = 0.85774;
 y[409][0] = b;

 x[410][0] = 0.001;
 x[410][1] = -2.4941;
 x[410][2] = 3.5447;
 x[410][3] = -1.3721;
 x[410][4] = -2.8483;
 y[410][0] = b;

 x[411][0] = 0.001;
 x[411][1] = -0.83121;
 x[411][2] = 0.039307;
 x[411][3] = 0.05369;
 x[411][4] = -0.23105;
 y[411][0] = b;

 x[412][0] = 0.001;
 x[412][1] = -2.5665;
 x[412][2] = -6.8824;
 x[412][3] = 7.5416;
 x[412][4] = 0.70774;
 y[412][0] = b;

 x[413][0] = 0.001;
 x[413][1] = -4.4018;
 x[413][2] = -12.9371;
 x[413][3] = 15.6559;
 x[413][4] = -1.6806;
 y[413][0] = b;

 x[414][0] = 0.001;
 x[414][1] = -3.7573;
 x[414][2] = -8.2916;
 x[414][3] = 10.3032;
 x[414][4] = 0.38059;
 y[414][0] = b;

 x[415][0] = 0.001;
 x[415][1] = -2.4725;
 x[415][2] = -0.40145;
 x[415][3] = 1.4855;
 x[415][4] = 1.1189;
 y[415][0] = b;

 x[416][0] = 0.001;
 x[416][1] = -1.9725;
 x[416][2] = 2.8825;
 x[416][3] = -2.3086;
 x[416][4] = -2.3724;
 y[416][0] = b;

 x[417][0] = 0.001;
 x[417][1] = -2.0149;
 x[417][2] = 3.6874;
 x[417][3] = -1.9385;
 x[417][4] = -3.8918;
 y[417][0] = b;

 x[418][0] = 0.001;
 x[418][1] = -0.82053;
 x[418][2] = 0.65181;
 x[418][3] = -0.48869;
 x[418][4] = -0.52716;
 y[418][0] = b;

 x[419][0] = 0.001;
 x[419][1] = -1.7886;
 x[419][2] = -6.3486;
 x[419][3] = 5.6154;
 x[419][4] = 0.42584;
 y[419][0] = b;

 x[420][0] = 0.001;
 x[420][1] = -2.9138;
 x[420][2] = -9.4711;
 x[420][3] = 9.7668;
 x[420][4] = -0.60216;
 y[420][0] = b;

 x[421][0] = 0.001;
 x[421][1] = -1.8343;
 x[421][2] = -6.5907;
 x[421][3] = 5.6429;
 x[421][4] = 0.54998;
 y[421][0] = b;

 x[422][0] = 0.001;
 x[422][1] = -0.8734;
 x[422][2] = -0.033118;
 x[422][3] = -0.20165;
 x[422][4] = 0.55774;
 y[422][0] = b;

 x[423][0] = 0.001;
 x[423][1] = -0.70346;
 x[423][2] = 2.957;
 x[423][3] = -3.5947;
 x[423][4] = -3.1457;
 y[423][0] = b;

 x[424][0] = 0.001;
 x[424][1] = -6.7387;
 x[424][2] = 6.9879;
 x[424][3] = 0.67833;
 x[424][4] = -7.5887;
 y[424][0] = b;

 x[425][0] = 0.001;
 x[425][1] = -2.7723;
 x[425][2] = 3.2777;
 x[425][3] = -0.9351;
 x[425][4] = -3.1457;
 y[425][0] = b;

 x[426][0] = 0.001;
 x[426][1] = -1.6641;
 x[426][2] = -1.3678;
 x[426][3] = 1.997;
 x[426][4] = 0.52283;
 y[426][0] = b;

 x[427][0] = 0.001;
 x[427][1] = -2.4349;
 x[427][2] = -9.2497;
 x[427][3] = 8.9922;
 x[427][4] = -0.50001;
 y[427][0] = b;

 x[428][0] = 0.001;
 x[428][1] = -3.793;
 x[428][2] = -12.7095;
 x[428][3] = 12.7957;
 x[428][4] = -2.825;
 y[428][0] = b;

 x[429][0] = 0.001;
 x[429][1] = -1.9551;
 x[429][2] = -6.9756;
 x[429][3] = 5.5383;
 x[429][4] = -0.12889;
 y[429][0] = b;

 x[430][0] = 0.001;
 x[430][1] = -0.69078;
 x[430][2] = -0.50077;
 x[430][3] = -0.35417;
 x[430][4] = 0.47498;
 y[430][0] = b;

 x[431][0] = 0.001;
 x[431][1] = 0.025013;
 x[431][2] = 3.3998;
 x[431][3] = -4.4327;
 x[431][4] = -4.2655;
 y[431][0] = b;

 x[432][0] = 0.001;
 x[432][1] = -4.3967;
 x[432][2] = 4.9601;
 x[432][3] = -0.64892;
 x[432][4] = -5.4719;
 y[432][0] = b;

 x[433][0] = 0.001;
 x[433][1] = -2.456;
 x[433][2] = -0.24418;
 x[433][3] = 1.4041;
 x[433][4] = -0.45863;
 y[433][0] = b;

 x[434][0] = 0.001;
 x[434][1] = -2.62;
 x[434][2] = -6.8555;
 x[434][3] = 6.2169;
 x[434][4] = -0.62285;
 y[434][0] = b;

 x[435][0] = 0.001;
 x[435][1] = -2.9662;
 x[435][2] = -10.3257;
 x[435][3] = 8.784;
 x[435][4] = -2.1138;
 y[435][0] = b;

 x[436][0] = 0.001;
 x[436][1] = -0.71494;
 x[436][2] = -4.4448;
 x[436][3] = 2.2241;
 x[436][4] = 0.49826;
 y[436][0] = b;

 x[437][0] = 0.001;
 x[437][1] = 0.6005;
 x[437][2] = 0.99945;
 x[437][3] = -2.2126;
 x[437][4] = 0.097399;
 y[437][0] = b;

 x[438][0] = 0.001;
 x[438][1] = 0.61652;
 x[438][2] = 3.8944;
 x[438][3] = -4.7275;
 x[438][4] = -4.3948;
 y[438][0] = b;

 x[439][0] = 0.001;
 x[439][1] = -5.4414;
 x[439][2] = 7.2363;
 x[439][3] = 0.10938;
 x[439][4] = -7.5642;
 y[439][0] = b;

 x[440][0] = 0.001;
 x[440][1] = -3.5798;
 x[440][2] = 0.45937;
 x[440][3] = 2.3457;
 x[440][4] = -0.45734;
 y[440][0] = b;

 x[441][0] = 0.001;
 x[441][1] = -2.7769;
 x[441][2] = -5.6967;
 x[441][3] = 5.9179;
 x[441][4] = 0.37671;
 y[441][0] = b;

 x[442][0] = 0.001;
 x[442][1] = -1.8356;
 x[442][2] = -6.7562;
 x[442][3] = 5.0585;
 x[442][4] = -0.55044;
 y[442][0] = b;

 x[443][0] = 0.001;
 x[443][1] = 0.30081;
 x[443][2] = 0.17381;
 x[443][3] = -1.7542;
 x[443][4] = 0.48921;
 y[443][0] = b;

 x[444][0] = 0.001;
 x[444][1] = 1.3403;
 x[444][2] = 4.1323;
 x[444][3] = -4.7018;
 x[444][4] = -2.5987;
 y[444][0] = b;

 x[445][0] = 0.001;
 x[445][1] = 0.26877;
 x[445][2] = 4.987;
 x[445][3] = -5.1508;
 x[445][4] = -6.3913;
 y[445][0] = b;

 x[446][0] = 0.001;
 x[446][1] = -6.5235;
 x[446][2] = 9.6014;
 x[446][3] = -0.25392;
 x[446][4] = -6.9642;
 y[446][0] = b;

 x[447][0] = 0.001;
 x[447][1] = -4.0679;
 x[447][2] = 2.4955;
 x[447][3] = 0.79571;
 x[447][4] = -1.1039;
 y[447][0] = b;

 x[448][0] = 0.001;
 x[448][1] = -2.564;
 x[448][2] = -1.7051;
 x[448][3] = 1.5026;
 x[448][4] = 0.32757;
 y[448][0] = b;

 x[449][0] = 0.001;
 x[449][1] = -1.3414;
 x[449][2] = -1.9162;
 x[449][3] = -0.15538;
 x[449][4] = -0.11984;
 y[449][0] = b;

 x[450][0] = 0.001;
 x[450][1] = 0.23874;
 x[450][2] = 2.0879;
 x[450][3] = -3.3522;
 x[450][4] = -0.66553;
 y[450][0] = b;

 x[451][0] = 0.001;
 x[451][1] = 0.6212;
 x[451][2] = 3.6771;
 x[451][3] = -4.0771;
 x[451][4] = -2.0711;
 y[451][0] = b;

 x[452][0] = 0.001;
 x[452][1] = -0.77848;
 x[452][2] = 3.4019;
 x[452][3] = -3.4859;
 x[452][4] = -3.5569;
 y[452][0] = b;

 x[453][0] = 0.001;
 x[453][1] = -4.1244;
 x[453][2] = 3.7909;
 x[453][3] = -0.6532;
 x[453][4] = -4.1802;
 y[453][0] = b;

 x[454][0] = 0.001;
 x[454][1] = -7.0421;
 x[454][2] = 9.2;
 x[454][3] = 0.25933;
 x[454][4] = -4.6832;
 y[454][0] = b;

 x[455][0] = 0.001;
 x[455][1] = -4.9462;
 x[455][2] = 3.5716;
 x[455][3] = 0.82742;
 x[455][4] = -1.4957;
 y[455][0] = b;

 x[456][0] = 0.001;
 x[456][1] = -3.5359;
 x[456][2] = 0.30417;
 x[456][3] = 0.6569;
 x[456][4] = -0.2957;
 y[456][0] = b;

 x[457][0] = 0.001;
 x[457][1] = -2.0662;
 x[457][2] = 0.16967;
 x[457][3] = -1.0054;
 x[457][4] = -0.82975;
 y[457][0] = b;

 x[458][0] = 0.001;
 x[458][1] = -0.88728;
 x[458][2] = 2.808;
 x[458][3] = -3.1432;
 x[458][4] = -1.2035;
 y[458][0] = b;

 x[459][0] = 0.001;
 x[459][1] = -1.0941;
 x[459][2] = 2.3072;
 x[459][3] = -2.5237;
 x[459][4] = -1.4453;
 y[459][0] = b;

 x[460][0] = 0.001;
 x[460][1] = -2.4458;
 x[460][2] = 1.6285;
 x[460][3] = -0.88541;
 x[460][4] = -1.4802;
 y[460][0] = b;

 x[461][0] = 0.001;
 x[461][1] = -3.551;
 x[461][2] = 1.8955;
 x[461][3] = 0.1865;
 x[461][4] = -2.4409;
 y[461][0] = b;

 x[462][0] = 0.001;
 x[462][1] = -2.2811;
 x[462][2] = -0.85669;
 x[462][3] = 2.7185;
 x[462][4] = 0.044382;
 y[462][0] = b;

 x[463][0] = 0.001;
 x[463][1] = -3.6053;
 x[463][2] = -5.974;
 x[463][3] = 10.0916;
 x[463][4] = -0.82846;
 y[463][0] = b;

 x[464][0] = 0.001;
 x[464][1] = -5.0676;
 x[464][2] = -5.1877;
 x[464][3] = 10.4266;
 x[464][4] = -0.86725;
 y[464][0] = b;

 x[465][0] = 0.001;
 x[465][1] = -3.9204;
 x[465][2] = 4.0723;
 x[465][3] = -0.23678;
 x[465][4] = -2.1151;
 y[465][0] = b;

 x[466][0] = 0.001;
 x[466][1] = -1.1306;
 x[466][2] = 1.8458;
 x[466][3] = -1.3575;
 x[466][4] = -1.3806;
 y[466][0] = b;

 x[467][0] = 0.001;
 x[467][1] = -2.4561;
 x[467][2] = -4.5566;
 x[467][3] = 6.4534;
 x[467][4] = -0.056479;
 y[467][0] = b;

 x[468][0] = 0.001;
 x[468][1] = -4.4775;
 x[468][2] = -13.0303;
 x[468][3] = 17.0834;
 x[468][4] = -3.0345;
 y[468][0] = b;

 x[469][0] = 0.001;
 x[469][1] = -4.1958;
 x[469][2] = -8.1819;
 x[469][3] = 12.1291;
 x[469][4] = -1.6017;
 y[469][0] = b;

 x[470][0] = 0.001;
 x[470][1] = -3.38;
 x[470][2] = -0.7077;
 x[470][3] = 2.5325;
 x[470][4] = 0.71808;
 y[470][0] = b;

 x[471][0] = 0.001;
 x[471][1] = -2.4365;
 x[471][2] = 3.6026;
 x[471][3] = -1.4166;
 x[471][4] = -2.8948;
 y[471][0] = b;



   }
   
   
   for(m=0;m<472;m++)
   {
      intermediaria(x,w1,h,m);
      saida(h,w2,o);
      err=erro_saida(o,y,m);
      fprintf(pFile,"\nPadrao>>%d",m);
      fprintf(pFile,"\ncalculado>>%f   \tdesejado>>%f  \tErro>>%2.9f",o[0],y[m][0],err);
      ////printf("\nPadrao>>%d",m); ///////////////// <------------
      printf("\ncalculado>>%f   \tdesejado>>%f  \tErro>>%2.9f",o[0],y[m][0],err); ///////////////// <------------
      emq=emq+err;

   }
   emq=emq/472;
   printf("\nemq>>%2.9f",emq);
   //fitness=(exp(-pow(emq,3)))+*pow(exp(-NINT),3)+(100.0/(NINT));
   fitness=1000*(exp(-emq)*exp(-NINT))+(1/(emq*NINT));
   fit[indic][0]=(indic+1);
   fit[indic][1]=fitness;
   fit[indic][2]=fit[indic][0];
   aux=indic+1;
   fit[indic][3]=aux;
   if(aux>1)
    fit[indic][3]=aux+fit[indic-1][3];
   else
     fit[indic][3]=aux;
   fit[indic][4]=emq;
   fit[indic][5]=(NINT);
   fit[indic][6]=((NINT)+(NENT)+NSAI);
   indic++;
//////   printf("\nemq>>%f",emq);

   fprintf(pFile,"\n\nErro Medio Quadratico>>%2.9f",emq);
//////   printf("\nfitness>>%2.9f",fitness);
   fprintf(pFile,"\n\nfitness>>%2.9f",fitness);
   fprintf(pFile,"\n\n<<Pesos Camada Entrada Oculta>>");

   for(i=0;i<NENT;i++)
   { fprintf(pFile,"\n");
     for(j=0;j<NINT;j++)
     {
      fprintf(pFile," w1[%d][%d]=%f",i,j,w1[i][j]);
     }
   }
   fprintf(pFile,"\n\n<<Pesos Camada Oculta Saida>>");
   for(i=0;i<NINT;i++)
   {
     fprintf(pFile,"\n");
     for(j=0;j<NSAI;j++)
     {
      fprintf(pFile," w2[%d][%d]=%f",i,j,w2[i][j]);
     }
   }
   fprintf(pFile,"\n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void treina_rede_(int individuos,FILE *pFile,int NINT1, int NINT2)
 {
//////   printf("\nNENT=(%d Entradas + 1 Bias)",NENT);
   fprintf(pFile,"\nNENT=(%d Entradas + 1 Bias)",NENT-1);
//////   printf("\nNINT1=(%d Int + 1 Bias)",NINT1-1);
   fprintf(pFile,"\nNINT1=(%d Int + 1 Bias)",NINT1-1);
//////   printf("\nNINT2=(%d Int + 1 Bias)",NINT2-1);
   fprintf(pFile,"\nNINT2=(%d Int + 1 Bias)",NINT2-1);
//////   printf("\nNSAI=%d",NSAI);
   fprintf(pFile,"\nNSAI=%d",NSAI);
//////   printf("\n\nTreinamento do Individuo=%d",individuos);
   fprintf(pFile,"\n\nTreinamento do Individuo=%d",individuos);

   double x[NPAD][NENT],h1[NINTER],h2[NINTER],o[NSAI],y[NPAD][NSAI],delta[NINTER],delta1[NINTER],delta2[NSAI],w1[NENT][NINTER],w2[NINTER][NSAI],w3[NINT1][NINT2],erro,erromax;
   int m;
   int l;
   /////////////////////////////////////////system("cls");
   fprintf(pFile,"\n<<Aprendizagem da Rede>>");
//////   printf("\n<<Aprendizagem da Rede\n>>");
   
   conjunto_treinamento1(x,y);
   inicializa_pesos_(w1,w2,w3);
   
   erromax = ACEITAVEL*2.0;
   
   m = 0;
   l = 0;
   
   while(erromax > ACEITAVEL && l < MAXITER)
   {
      if(m == NPAD) m = 0;
      intermediaria_(x,w1,h1,m);
      intermediaria2_(w3,h1,h2,m);
      saida_(h2,w2,o);
      
      err=erro_saida_(o,y,m);
      
      if(err > erromax) erromax = err;
      
	  erro2_(o,y,m,delta2);
      erro1_(h2,delta2,w2,delta1);
      erroN_(h1,delta1,w3,delta);
	  ajusta2_(w2,delta2,h2);
      ajusta1_(w1,delta,x,m);
      ajusta_(w3,delta1,h1,m);
      
      l = l + 1;
      
//////      printf("\nPadrao>>%d",m);
      
	  m = m + 1;
      
	//  printf("\t\t Epoca:%d \t\t Erro:%2.9f",l,err);
    
	}
   // getch();
//////    printf("\a\a\a\a\a\a\a\a\a");
	fprintf(pFile,"\n<<Verificacao da Aprendizagem>>\n");
    
	//printf("\n<<Verificacao da Aprendizagem>>\n");
    
	verifica1_(x,w1,w2,w3,y,pFile);
   // getch();
 }
 
void inicializa_pesos_(double w1[][NINTER],double w2[][NSAI], double w3[NINT1][NINT2])
{
   int i,j,k;
   double aleatorio;
   for(i=0;i<NENT;i++)
     for(j=1;j<NINT1;j++)
     {
	    aleatorio = (rand()%101);
		w1[i][j] =  (1.0 - 2.0 * aleatorio) / 100.0;
     }
   for(j=0;j<NINT2;j++)
     for(k=0;k<NSAI;k++)
       {
	     aleatorio = (rand()%101);
	     w2[j][k] = (1.0 - 2.0 * aleatorio) / 100.0;
       }
    for(j=0;j<NINT1;j++)
     for(k=0;k<NINT2;k++)
       {
	     aleatorio = (rand()%101);
	     w3[j][k] = (1.0 - 2.0 * aleatorio) / 100.0;
       }
}

void intermediaria_(double x[][NENT],double w1[][NINTER],double h1[],int m)
  {
   int i,j;
   double somatorio;
   h1[0] = 1.0;
   for(j=1;j<NINT1;j++)
     {
      somatorio = 0.0;
      for(i=0;i<NENT;i++)
	     somatorio = somatorio + x[m][i] * w1[i][j];
	     ////printf("int1 %lf\n", somatorio);
        
      
      h1[j] =1.0 / (1.0 + exp(-somatorio));
    
     }
  }

void intermediaria2_(double w3[NINT1][NINT2], double h1[],double h2[] ,int m)
  {
   int i,j;
   double somatorio;
   h2[0] = 1.0;
   for(j=1;j<NINT2;j++)
     {
      somatorio = 0.0;
      for(i=0;i<NINT1;i++)
	     somatorio = somatorio + h1[i] * w3[i][j];
	     //printf("int2 %lf\n", somatorio);
      
      h2[j] = 1.0 / (1.0 + exp(-somatorio));
     }
  }
  
  void saida_(double h2[],double w2[][NSAI],double o[])
  {
   int j,k;
   double somatorio;
   for(k=0;k<NSAI;k++)
     {
      somatorio = 0.0;
      for(j=0;j<NINT2;j++)
    	somatorio = somatorio + h2[j] * w2[j][k];
    	//printf("sai %lf\n", somatorio);
      o[k] = 1.0 / (1.0 + exp(-somatorio));//1.0 / (1.0 + exp(somatorio));
     }
  }
  
  
  double erro_saida_(double o[],double y[][NSAI],int m)
{
   int k;
   double somatorio,erro;
   somatorio = 0.0;
   for(k=0;k<NSAI;k++)
     somatorio = somatorio + (o[k] - y[m][k]) * (o[k] - y[m][k]);
   erro = 0.5 * somatorio;

   return erro;
}

void erro2_(double o[],double y[][NSAI],int m,double delta2[])
  {
   int k;
   for(k=0;k<NSAI;k++)
     delta2[k] = o[k] * (1.0 - o[k]) * (y[m][k] - o[k]);
  }
  
  void erro1_(double h2[],double delta2[],double w2[][NSAI],double delta1[])
  {
   int j,k;
   double somatorio;
   for(j=1;j<NINT2;j++)
     {
      somatorio = 0.0;
      for(k=0;k<NSAI;k++)
    	somatorio = somatorio + delta2[k] * w2[j][k];
      delta1[j] = h2[j] * (1 - h2[j]) * somatorio;
     }
  }
  
    void erroN_(double h1[],double delta1[],double w3[][NSAI],double delta[])
  {
   int j,k;
   double somatorio;
   for(j=1;j<NINT1;j++)
     {
      somatorio = 0.0;
      for(k=0;k<NINT2;k++)
    	somatorio = somatorio + delta1[k] * w3[j][k];
      delta[j] = h1[j] * (1 - h1[j]) * somatorio;
     }
  }
  
  void ajusta2_(double w2[][NSAI],double delta2[],double h2[])
  {
   int j,k;
   for(j=1;j<NINT2;j++)
     for(k=0;k<NSAI;k++)
       w2[j][k] = w2[j][k] + TAPR4 * delta2[k] * h2[j];
  }
  
  void ajusta1_(double w1[][NINTER],double delta[],double x[][NENT],int m)
  {
   int i,j;
   for(i=0;i<NENT;i++)
     for(j=1;j<NINT1;j++)
       w1[i][j] = w1[i][j] + TAPR4 * delta[j] * x[m][i];
  }
  
void ajusta_(double w3[][NINTER],double delta1[],double h1[],int m)
  {
   int i,j;
   for(i=0;i<NINT1;i++)
     for(j=1;j<NINT2;j++)
       w3[i][j] = w3[i][j] + TAPR4 * delta1[j] * h1[i];
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void verifica1_(double x[][NENT],double w1[][NINTER],double w2[][NSAI], double w3[NINT1][NINT2],double y[][NSAI],FILE *pFile)
{
   int m,i,j;
   double err,h1[NINTER],o[NSAI], h2[NINTER];

   float a = 1, b = 0;
   {
   	
//TESTE

 x[0][0] = 0.001;
 x[0][1] = 2.9719;
 x[0][2] = 6.8369;
 x[0][3] = -0.2702;
 x[0][4] = 0.71291;
 y[0][0] = a;

 x[1][0] = 0.001;
 x[1][1] = 1.6849;
 x[1][2] = 8.7489;
 x[1][3] = -1.2641;
 x[1][4] = -1.3858;
 y[1][0] = a;

 x[2][0] = 0.001;
 x[2][1] = -1.9177;
 x[2][2] = 11.6894;
 x[2][3] = 2.5454;
 x[2][4] = -3.2763;
 y[2][0] = a;

 x[3][0] = 0.001;
 x[3][1] = 2.3729;
 x[3][2] = 10.4726;
 x[3][3] = -3.0087;
 x[3][4] = -3.2013;
 y[3][0] = a;

 x[4][0] = 0.001;
 x[4][1] = 1.0284;
 x[4][2] = 9.767;
 x[4][3] = -1.3687;
 x[4][4] = -1.7853;
 y[4][0] = a;

 x[5][0] = 0.001;
 x[5][1] = 0.27451;
 x[5][2] = 9.2186;
 x[5][3] = -3.2863;
 x[5][4] = -4.8448;
 y[5][0] = a;

 x[6][0] = 0.001;
 x[6][1] = 1.6032;
 x[6][2] = -4.7863;
 x[6][3] = 8.5193;
 x[6][4] = -2.1203;
 y[6][0] = a;

 x[7][0] = 0.001;
 x[7][1] = 4.616;
 x[7][2] = 10.1788;
 x[7][3] = -4.2185;
 x[7][4] = -4.4245;
 y[7][0] = a;

 x[8][0] = 0.001;
 x[8][1] = 4.2478;
 x[8][2] = 7.6956;
 x[8][3] = -2.7696;
 x[8][4] = -1.0767;
 y[8][0] = a;

 x[9][0] = 0.001;
 x[9][1] = 4.0215;
 x[9][2] = -2.7004;
 x[9][3] = 2.4957;
 x[9][4] = 0.36636;
 y[9][0] = a;

 x[10][0] = 0.001;
 x[10][1] = 5.0297;
 x[10][2] = -4.9704;
 x[10][3] = 3.5025;
 x[10][4] = -0.23751;
 y[10][0] = a;

 x[11][0] = 0.001;
 x[11][1] = 1.5902;
 x[11][2] = 2.2948;
 x[11][3] = 3.2403;
 x[11][4] = 0.18404;
 y[11][0] = a;

 x[12][0] = 0.001;
 x[12][1] = 2.1274;
 x[12][2] = 5.1939;
 x[12][3] = -1.7971;
 x[12][4] = -1.1763;
 y[12][0] = a;

 x[13][0] = 0.001;
 x[13][1] = 1.1811;
 x[13][2] = 8.3847;
 x[13][3] = -2.0567;
 x[13][4] = -0.90345;
 y[13][0] = a;

 x[14][0] = 0.001;
 x[14][1] = 0.3292;
 x[14][2] = -4.4552;
 x[14][3] = 4.5718;
 x[14][4] = -0.9888;
 y[14][0] = a;

 x[15][0] = 0.001;
 x[15][1] = 5.7353;
 x[15][2] = 5.2808;
 x[15][3] = -2.2598;
 x[15][4] = 0.075416;
 y[15][0] = a;

 x[16][0] = 0.001;
 x[16][1] = 2.6718;
 x[16][2] = 5.6574;
 x[16][3] = 0.72974;
 x[16][4] = -1.4892;
 y[16][0] = a;

 x[17][0] = 0.001;
 x[17][1] = 1.5799;
 x[17][2] = -4.7076;
 x[17][3] = 7.9186;
 x[17][4] = -1.5487;
 y[17][0] = a;

 x[18][0] = 0.001;
 x[18][1] = 2.9499;
 x[18][2] = 2.2493;
 x[18][3] = 1.3458;
 x[18][4] = -0.037083;
 y[18][0] = a;

 x[19][0] = 0.001;
 x[19][1] = 0.5195;
 x[19][2] = -3.2633;
 x[19][3] = 3.0895;
 x[19][4] = -0.9849;
 y[19][0] = a;

 x[20][0] = 0.001;
 x[20][1] = 3.7352;
 x[20][2] = 9.5911;
 x[20][3] = -3.9032;
 x[20][4] = -3.3487;
 y[20][0] = a;

 x[21][0] = 0.001;
 x[21][1] = -1.7344;
 x[21][2] = 2.0175;
 x[21][3] = 7.7618;
 x[21][4] = 0.93532;
 y[21][0] = a;

 x[22][0] = 0.001;
 x[22][1] = 3.884;
 x[22][2] = 10.0277;
 x[22][3] = -3.9298;
 x[22][4] = -4.0819;
 y[22][0] = a;

 x[23][0] = 0.001;
 x[23][1] = 3.5257;
 x[23][2] = 1.2829;
 x[23][3] = 1.9276;
 x[23][4] = 1.7991;
 y[23][0] = a;

 x[24][0] = 0.001;
 x[24][1] = 4.4549;
 x[24][2] = 2.4976;
 x[24][3] = 1.0313;
 x[24][4] = 0.96894;
 y[24][0] = a;

 x[25][0] = 0.001;
 x[25][1] = -0.16108;
 x[25][2] = -6.4624;
 x[25][3] = 8.3573;
 x[25][4] = -1.5216;
 y[25][0] = a;

 x[26][0] = 0.001;
 x[26][1] = 4.2164;
 x[26][2] = 9.4607;
 x[26][3] = -4.9288;
 x[26][4] = -5.2366;
 y[26][0] = a;

 x[27][0] = 0.001;
 x[27][1] = 3.5152;
 x[27][2] = 6.8224;
 x[27][3] = -0.67377;
 x[27][4] = -0.46898;
 y[27][0] = a;

 x[28][0] = 0.001;
 x[28][1] = 1.6988;
 x[28][2] = 2.9094;
 x[28][3] = 2.9044;
 x[28][4] = 0.11033;
 y[28][0] = a;

 x[29][0] = 0.001;
 x[29][1] = 1.0607;
 x[29][2] = 2.4542;
 x[29][3] = 2.5188;
 x[29][4] = -0.17027;
 y[29][0] = a;

 x[30][0] = 0.001;
 x[30][1] = 2.0421;
 x[30][2] = 1.2436;
 x[30][3] = 4.2171;
 x[30][4] = 0.90429;
 y[30][0] = a;

 x[31][0] = 0.001;
 x[31][1] = 3.5594;
 x[31][2] = 1.3078;
 x[31][3] = 1.291;
 x[31][4] = 1.6556;
 y[31][0] = a;

 x[32][0] = 0.001;
 x[32][1] = 3.0009;
 x[32][2] = 5.8126;
 x[32][3] = -2.2306;
 x[32][4] = -0.66553;
 y[32][0] = a;

 x[33][0] = 0.001;
 x[33][1] = 3.9294;
 x[33][2] = 1.4112;
 x[33][3] = 1.8076;
 x[33][4] = 0.89782;
 y[33][0] = a;

 x[34][0] = 0.001;
 x[34][1] = 3.4667;
 x[34][2] = -4.0724;
 x[34][3] = 4.2882;
 x[34][4] = 1.5418;
 y[34][0] = a;

 x[35][0] = 0.001;
 x[35][1] = 3.966;
 x[35][2] = 3.9213;
 x[35][3] = 0.70574;
 x[35][4] = 0.33662;
 y[35][0] = a;

 x[36][0] = 0.001;
 x[36][1] = 1.0191;
 x[36][2] = 2.33;
 x[36][3] = 4.9334;
 x[36][4] = 0.82929;
 y[36][0] = a;

 x[37][0] = 0.001;
 x[37][1] = 0.96414;
 x[37][2] = 5.616;
 x[37][3] = 2.2138;
 x[37][4] = -0.12501;
 y[37][0] = a;

 x[38][0] = 0.001;
 x[38][1] = 1.8205;
 x[38][2] = 6.7562;
 x[38][3] = 0.009991;
 x[38][4] = 0.39481;
 y[38][0] = a;

 x[39][0] = 0.001;
 x[39][1] = 4.9923;
 x[39][2] = 7.8653;
 x[39][3] = -2.3515;
 x[39][4] = -0.71984;
 y[39][0] = a;

 x[40][0] = 0.001;
 x[40][1] = -1.1804;
 x[40][2] = 11.5093;
 x[40][3] = 0.15565;
 x[40][4] = -6.8194;
 y[40][0] = a;

 x[41][0] = 0.001;
 x[41][1] = 4.0329;
 x[41][2] = 0.23175;
 x[41][3] = 0.89082;
 x[41][4] = 1.1823;
 y[41][0] = a;

 x[42][0] = 0.001;
 x[42][1] = 0.66018;
 x[42][2] = 10.3878;
 x[42][3] = -1.4029;
 x[42][4] = -3.9151;
 y[42][0] = a;

 x[43][0] = 0.001;
 x[43][1] = 3.5982;
 x[43][2] = 7.1307;
 x[43][3] = -1.3035;
 x[43][4] = 0.21248;
 y[43][0] = a;

 x[44][0] = 0.001;
 x[44][1] = -1.8584;
 x[44][2] = 7.886;
 x[44][3] = -1.6643;
 x[44][4] = -1.8384;
 y[44][0] = a;

 x[45][0] = 0.001;
 x[45][1] = 4.0972;
 x[45][2] = 0.46972;
 x[45][3] = 1.6671;
 x[45][4] = 0.91593;
 y[45][0] = a;

 x[46][0] = 0.001;
 x[46][1] = 3.3299;
 x[46][2] = 0.91254;
 x[46][3] = 1.5806;
 x[46][4] = 0.39352;
 y[46][0] = a;

 x[47][0] = 0.001;
 x[47][1] = 3.1088;
 x[47][2] = 3.1122;
 x[47][3] = 0.80857;
 x[47][4] = 0.4336;
 y[47][0] = a;

 x[48][0] = 0.001;
 x[48][1] = -4.2859;
 x[48][2] = 8.5234;
 x[48][3] = 3.1392;
 x[48][4] = -0.91639;
 y[48][0] = a;

 x[49][0] = 0.001;
 x[49][1] = -1.2528;
 x[49][2] = 10.2036;
 x[49][3] = 2.1787;
 x[49][4] = -5.6038;
 y[49][0] = a;

 x[50][0] = 0.001;
 x[50][1] = 0.5195;
 x[50][2] = -3.2633;
 x[50][3] = 3.0895;
 x[50][4] = -0.9849;
 y[50][0] = a;

 x[51][0] = 0.001;
 x[51][1] = 0.3292;
 x[51][2] = -4.4552;
 x[51][3] = 4.5718;
 x[51][4] = -0.9888;
 y[51][0] = a;

 x[52][0] = 0.001;
 x[52][1] = 0.88872;
 x[52][2] = 5.3449;
 x[52][3] = 2.045;
 x[52][4] = -0.19355;
 y[52][0] = a;

 x[53][0] = 0.001;
 x[53][1] = 3.5458;
 x[53][2] = 9.3718;
 x[53][3] = -4.0351;
 x[53][4] = -3.9564;
 y[53][0] = a;

 x[54][0] = 0.001;
 x[54][1] = -0.21661;
 x[54][2] = 8.0329;
 x[54][3] = 1.8848;
 x[54][4] = -3.8853;
 y[54][0] = a;

 x[55][0] = 0.001;
 x[55][1] = 2.7206;
 x[55][2] = 9.0821;
 x[55][3] = -3.3111;
 x[55][4] = -0.96811;
 y[55][0] = a;

 x[56][0] = 0.001;
 x[56][1] = 3.2051;
 x[56][2] = 8.6889;
 x[56][3] = -2.9033;
 x[56][4] = -0.7819;
 y[56][0] = a;

 x[57][0] = 0.001;
 x[57][1] = 2.6917;
 x[57][2] = 10.8161;
 x[57][3] = -3.3;
 x[57][4] = -4.2888;
 y[57][0] = a;

 x[58][0] = 0.001;
 x[58][1] = -2.3242;
 x[58][2] = 11.5176;
 x[58][3] = 1.8231;
 x[58][4] = -5.375;
 y[58][0] = a;

 x[59][0] = 0.001;
 x[59][1] = 2.7161;
 x[59][2] = -4.2006;
 x[59][3] = 4.1914;
 x[59][4] = 0.16981;
 y[59][0] = a;

 x[60][0] = 0.001;
 x[60][1] = 3.3848;
 x[60][2] = 3.2674;
 x[60][3] = 0.90967;
 x[60][4] = 0.25128;
 y[60][0] = a;

 x[61][0] = 0.001;
 x[61][1] = 1.7452;
 x[61][2] = 4.8028;
 x[61][3] = 2.0878;
 x[61][4] = 0.62627;
 y[61][0] = a;

 x[62][0] = 0.001;
 x[62][1] = 2.805;
 x[62][2] = 0.57732;
 x[62][3] = 1.3424;
 x[62][4] = 1.2133;
 y[62][0] = a;

 x[63][0] = 0.001;
 x[63][1] = 5.7823;
 x[63][2] = 5.5788;
 x[63][3] = -2.4089;
 x[63][4] = -0.056479;
 y[63][0] = a;

 x[64][0] = 0.001;
 x[64][1] = 3.8999;
 x[64][2] = 1.734;
 x[64][3] = 1.6011;
 x[64][4] = 0.96765;
 y[64][0] = a;

 x[65][0] = 0.001;
 x[65][1] = 3.5189;
 x[65][2] = 6.332;
 x[65][3] = -1.7791;
 x[65][4] = -0.020273;
 y[65][0] = a;

 x[66][0] = 0.001;
 x[66][1] = 3.2294;
 x[66][2] = 7.7391;
 x[66][3] = -0.37816;
 x[66][4] = -2.5405;
 y[66][0] = a;

 x[67][0] = 0.001;
 x[67][1] = 3.4985;
 x[67][2] = 3.1639;
 x[67][3] = 0.22677;
 x[67][4] = -0.1651;
 y[67][0] = a;

 x[68][0] = 0.001;
 x[68][1] = 2.1948;
 x[68][2] = 1.3781;
 x[68][3] = 1.1582;
 x[68][4] = 0.85774;
 y[68][0] = a;

 x[69][0] = 0.001;
 x[69][1] = 2.2526;
 x[69][2] = 9.9636;
 x[69][3] = -3.1749;
 x[69][4] = -2.9944;
 y[69][0] = a;

 x[70][0] = 0.001;
 x[70][1] = 4.1529;
 x[70][2] = -3.9358;
 x[70][3] = 2.8633;
 x[70][4] = -0.017686;
 y[70][0] = a;

 x[71][0] = 0.001;
 x[71][1] = 0.74307;
 x[71][2] = 11.17;
 x[71][3] = -1.3824;
 x[71][4] = -4.0728;
 y[71][0] = a;

 x[72][0] = 0.001;
 x[72][1] = 1.9105;
 x[72][2] = 8.871;
 x[72][3] = -2.3386;
 x[72][4] = -0.75604;
 y[72][0] = a;

 x[73][0] = 0.001;
 x[73][1] = -1.5055;
 x[73][2] = 0.070346;
 x[73][3] = 6.8681;
 x[73][4] = -0.50648;
 y[73][0] = a;

 x[74][0] = 0.001;
 x[74][1] = 0.58836;
 x[74][2] = 10.7727;
 x[74][3] = -1.3884;
 x[74][4] = -4.3276;
 y[74][0] = a;

 x[75][0] = 0.001;
 x[75][1] = 3.2303;
 x[75][2] = 7.8384;
 x[75][3] = -3.5348;
 x[75][4] = -1.2151;
 y[75][0] = a;

 x[76][0] = 0.001;
 x[76][1] = -1.9922;
 x[76][2] = 11.6542;
 x[76][3] = 2.6542;
 x[76][4] = -5.2107;
 y[76][0] = a;

 x[77][0] = 0.001;
 x[77][1] = 2.8523;
 x[77][2] = 9.0096;
 x[77][3] = -3.761;
 x[77][4] = -3.3371;
 y[77][0] = a;

 x[78][0] = 0.001;
 x[78][1] = 4.2772;
 x[78][2] = 2.4955;
 x[78][3] = 0.48554;
 x[78][4] = 0.36119;
 y[78][0] = a;

 x[79][0] = 0.001;
 x[79][1] = 1.5099;
 x[79][2] = 0.039307;
 x[79][3] = 6.2332;
 x[79][4] = -0.30346;
 y[79][0] = a;

 x[80][0] = 0.001;
 x[80][1] = 5.4188;
 x[80][2] = 10.1457;
 x[80][3] = -4.084;
 x[80][4] = -3.6991;
 y[80][0] = a;

 x[81][0] = 0.001;
 x[81][1] = 0.86202;
 x[81][2] = 2.6963;
 x[81][3] = 4.2908;
 x[81][4] = 0.54739;
 y[81][0] = a;

 x[82][0] = 0.001;
 x[82][1] = 3.8117;
 x[82][2] = 10.1457;
 x[82][3] = -4.0463;
 x[82][4] = -4.5629;
 y[82][0] = a;

 x[83][0] = 0.001;
 x[83][1] = 0.54777;
 x[83][2] = 10.3754;
 x[83][3] = -1.5435;
 x[83][4] = -4.1633;
 y[83][0] = a;

 x[84][0] = 0.001;
 x[84][1] = 2.3718;
 x[84][2] = 7.4908;
 x[84][3] = 0.015989;
 x[84][4] = -1.7414;
 y[84][0] = a;

 x[85][0] = 0.001;
 x[85][1] = -2.4953;
 x[85][2] = 11.1472;
 x[85][3] = 1.9353;
 x[85][4] = -3.4638;
 y[85][0] = a;

 x[86][0] = 0.001;
 x[86][1] = 4.6361;
 x[86][2] = -2.6611;
 x[86][3] = 2.8358;
 x[86][4] = 1.1991;
 y[86][0] = a;

 x[87][0] = 0.001;
 x[87][1] = -2.2527;
 x[87][2] = 11.5321;
 x[87][3] = 2.5899;
 x[87][4] = -3.2737;
 y[87][0] = a;

 x[88][0] = 0.001;
 x[88][1] = 3.7982;
 x[88][2] = 10.423;
 x[88][3] = -4.1602;
 x[88][4] = -4.9728;
 y[88][0] = a;

 x[89][0] = 0.001;
 x[89][1] = -0.36279;
 x[89][2] = 8.2895;
 x[89][3] = -1.9213;
 x[89][4] = -3.3332;
 y[89][0] = a;

 x[90][0] = 0.001;
 x[90][1] = 2.1265;
 x[90][2] = 6.8783;
 x[90][3] = 0.44784;
 x[90][4] = -2.2224;
 y[90][0] = a;

 x[91][0] = 0.001;
 x[91][1] = 0.86736;
 x[91][2] = 5.5643;
 x[91][3] = 1.6765;
 x[91][4] = -0.16769;
 y[91][0] = a;

 x[92][0] = 0.001;
 x[92][1] = 3.7831;
 x[92][2] = 10.0526;
 x[92][3] = -3.8869;
 x[92][4] = -3.7366;
 y[92][0] = a;

 x[93][0] = 0.001;
 x[93][1] = -2.2623;
 x[93][2] = 12.1177;
 x[93][3] = 0.28846;
 x[93][4] = -7.7581;
 y[93][0] = a;

 x[94][0] = 0.001;
 x[94][1] = 1.2616;
 x[94][2] = 4.4303;
 x[94][3] = -1.3335;
 x[94][4] = -1.7517;
 y[94][0] = a;

 x[95][0] = 0.001;
 x[95][1] = 2.6799;
 x[95][2] = 3.1349;
 x[95][3] = 0.34073;
 x[95][4] = 0.58489;
 y[95][0] = a;

 x[96][0] = 0.001;
 x[96][1] = -0.39816;
 x[96][2] = 5.9781;
 x[96][3] = 1.3912;
 x[96][4] = -1.1621;
 y[96][0] = a;

 x[97][0] = 0.001;
 x[97][1] = 4.3937;
 x[97][2] = 0.35798;
 x[97][3] = 2.0416;
 x[97][4] = 1.2004;
 y[97][0] = a;

 x[98][0] = 0.001;
 x[98][1] = 2.9695;
 x[98][2] = 5.6222;
 x[98][3] = 0.27561;
 x[98][4] = -1.1556;
 y[98][0] = a;

 x[99][0] = 0.001;
 x[99][1] = 1.3049;
 x[99][2] = -0.15521;
 x[99][3] = 6.4911;
 x[99][4] = -0.75346;
 y[99][0] = a;

 x[100][0] = 0.001;
 x[100][1] = 2.2123;
 x[100][2] = -5.8395;
 x[100][3] = 7.7687;
 x[100][4] = -0.85302;
 y[100][0] = a;

 x[101][0] = 0.001;
 x[101][1] = 1.9647;
 x[101][2] = 6.9383;
 x[101][3] = 0.57722;
 x[101][4] = 0.66377;
 y[101][0] = a;

 x[102][0] = 0.001;
 x[102][1] = 3.0864;
 x[102][2] = -2.5845;
 x[102][3] = 2.2309;
 x[102][4] = 0.30947;
 y[102][0] = a;

 x[103][0] = 0.001;
 x[103][1] = 0.3798;
 x[103][2] = 0.7098;
 x[103][3] = 0.7572;
 x[103][4] = -0.4444;
 y[103][0] = a;

 x[104][0] = 0.001;
 x[104][1] = 0.58982;
 x[104][2] = 7.4266;
 x[104][3] = 1.2353;
 x[104][4] = -2.9595;
 y[104][0] = a;

 x[105][0] = 0.001;
 x[105][1] = 0.14783;
 x[105][2] = 7.946;
 x[105][3] = 1.0742;
 x[105][4] = -3.3409;
 y[105][0] = a;

 x[106][0] = 0.001;
 x[106][1] = -0.062025;
 x[106][2] = 6.1975;
 x[106][3] = 1.099;
 x[106][4] = -1.131;
 y[106][0] = a;

 x[107][0] = 0.001;
 x[107][1] = 4.223;
 x[107][2] = 1.1319;
 x[107][3] = 0.72202;
 x[107][4] = 0.96118;
 y[107][0] = a;

 x[108][0] = 0.001;
 x[108][1] = 0.64295;
 x[108][2] = 7.1018;
 x[108][3] = 0.3493;
 x[108][4] = -0.41337;
 y[108][0] = a;

 x[109][0] = 0.001;
 x[109][1] = 1.941;
 x[109][2] = 0.46351;
 x[109][3] = 4.6472;
 x[109][4] = 1.0879;
 y[109][0] = a;

 x[110][0] = 0.001;
 x[110][1] = 4.0047;
 x[110][2] = 0.45937;
 x[110][3] = 1.3621;
 x[110][4] = 1.6181;
 y[110][0] = a;

 x[111][0] = 0.001;
 x[111][1] = 3.7767;
 x[111][2] = 9.7794;
 x[111][3] = -3.9075;
 x[111][4] = -3.5323;
 y[111][0] = a;

 x[112][0] = 0.001;
 x[112][1] = 3.4769;
 x[112][2] = -0.15314;
 x[112][3] = 2.53;
 x[112][4] = 2.4495;
 y[112][0] = a;

 x[113][0] = 0.001;
 x[113][1] = 1.9818;
 x[113][2] = 9.2621;
 x[113][3] = -3.521;
 x[113][4] = -1.872;
 y[113][0] = a;

 x[114][0] = 0.001;
 x[114][1] = 3.8023;
 x[114][2] = -3.8696;
 x[114][3] = 4.044;
 x[114][4] = 0.95343;
 y[114][0] = a;

 x[115][0] = 0.001;
 x[115][1] = 4.3483;
 x[115][2] = 11.1079;
 x[115][3] = -4.0857;
 x[115][4] = -4.2539;
 y[115][0] = a;

 x[116][0] = 0.001;
 x[116][1] = 1.1518;
 x[116][2] = 1.3864;
 x[116][3] = 5.2727;
 x[116][4] = -0.43536;
 y[116][0] = a;

 x[117][0] = 0.001;
 x[117][1] = -1.2576;
 x[117][2] = 1.5892;
 x[117][3] = 7.0078;
 x[117][4] = 0.42455;
 y[117][0] = a;

 x[118][0] = 0.001;
 x[118][1] = 1.9572;
 x[118][2] = -5.1153;
 x[118][3] = 8.6127;
 x[118][4] = -1.4297;
 y[118][0] = a;

 x[119][0] = 0.001;
 x[119][1] = -2.484;
 x[119][2] = 12.1611;
 x[119][3] = 2.8204;
 x[119][4] = -3.7418;
 y[119][0] = a;

 x[120][0] = 0.001;
 x[120][1] = -1.1497;
 x[120][2] = 1.2954;
 x[120][3] = 7.701;
 x[120][4] = 0.62627;
 y[120][0] = a;

 x[121][0] = 0.001;
 x[121][1] = 4.8368;
 x[121][2] = 10.0132;
 x[121][3] = -4.3239;
 x[121][4] = -4.3276;
 y[121][0] = a;

 x[122][0] = 0.001;
 x[122][1] = -0.12196;
 x[122][2] = 8.8068;
 x[122][3] = 0.94566;
 x[122][4] = -4.2267;
 y[122][0] = a;

 x[123][0] = 0.001;
 x[123][1] = 1.9429;
 x[123][2] = 6.3961;
 x[123][3] = 0.092248;
 x[123][4] = 0.58102;
 y[123][0] = a;

 x[124][0] = 0.001;
 x[124][1] = 1.742;
 x[124][2] = -4.809;
 x[124][3] = 8.2142;
 x[124][4] = -2.0659;
 y[124][0] = a;

 x[125][0] = 0.001;
 x[125][1] = -1.5222;
 x[125][2] = 10.8409;
 x[125][3] = 2.7827;
 x[125][4] = -4.0974;
 y[125][0] = a;

 x[126][0] = 0.001;
 x[126][1] = -1.3;
 x[126][2] = 10.2678;
 x[126][3] = -2.953;
 x[126][4] = -5.8638;
 y[126][0] = a;

 x[127][0] = 0.001;
 x[127][1] = 3.4246;
 x[127][2] = -0.14693;
 x[127][3] = 0.80342;
 x[127][4] = 0.29136;
 y[127][0] = a;

 x[128][0] = 0.001;
 x[128][1] = 2.5503;
 x[128][2] = -4.9518;
 x[128][3] = 6.3729;
 x[128][4] = -0.41596;
 y[128][0] = a;

 x[129][0] = 0.001;
 x[129][1] = 1.5691;
 x[129][2] = 6.3465;
 x[129][3] = -0.1828;
 x[129][4] = -2.4099;
 y[129][0] = a;

 x[130][0] = 0.001;
 x[130][1] = 1.3087;
 x[130][2] = 4.9228;
 x[130][3] = 2.0013;
 x[130][4] = 0.22024;
 y[130][0] = a;

 x[131][0] = 0.001;
 x[131][1] = 5.1776;
 x[131][2] = 8.2316;
 x[131][3] = -3.2511;
 x[131][4] = -1.5694;
 y[131][0] = a;

 x[132][0] = 0.001;
 x[132][1] = 2.229;
 x[132][2] = 9.6325;
 x[132][3] = -3.1123;
 x[132][4] = -2.7164;
 y[132][0] = a;

 x[133][0] = 0.001;
 x[133][1] = 5.6272;
 x[133][2] = 10.0857;
 x[133][3] = -4.2931;
 x[133][4] = -3.8142;
 y[133][0] = a;

 x[134][0] = 0.001;
 x[134][1] = 1.2138;
 x[134][2] = 8.7986;
 x[134][3] = -2.1672;
 x[134][4] = -0.74182;
 y[134][0] = a;

 x[135][0] = 0.001;
 x[135][1] = 0.3798;
 x[135][2] = 0.7098;
 x[135][3] = 0.7572;
 x[135][4] = -0.4444;
 y[135][0] = a;

 x[136][0] = 0.001;
 x[136][1] = 0.5415;
 x[136][2] = 6.0319;
 x[136][3] = 1.6825;
 x[136][4] = -0.46122;
 y[136][0] = a;

 x[137][0] = 0.001;
 x[137][1] = 4.0524;
 x[137][2] = 5.6802;
 x[137][3] = -1.9693;
 x[137][4] = 0.026279;
 y[137][0] = a;

 x[138][0] = 0.001;
 x[138][1] = 4.7285;
 x[138][2] = 2.1065;
 x[138][3] = -0.28305;
 x[138][4] = 1.5625;
 y[138][0] = a;

 x[139][0] = 0.001;
 x[139][1] = 3.4359;
 x[139][2] = 0.66216;
 x[139][3] = 2.1041;
 x[139][4] = 1.8922;
 y[139][0] = a;

 x[140][0] = 0.001;
 x[140][1] = 0.86816;
 x[140][2] = 10.2429;
 x[140][3] = -1.4912;
 x[140][4] = -4.0082;
 y[140][0] = a;

 x[141][0] = 0.001;
 x[141][1] = 3.359;
 x[141][2] = 9.8022;
 x[141][3] = -3.8209;
 x[141][4] = -3.7133;
 y[141][0] = a;

 x[142][0] = 0.001;
 x[142][1] = 3.6702;
 x[142][2] = 2.9942;
 x[142][3] = 0.85141;
 x[142][4] = 0.30688;
 y[142][0] = a;

 x[143][0] = 0.001;
 x[143][1] = 1.3349;
 x[143][2] = 6.1189;
 x[143][3] = 0.46497;
 x[143][4] = 0.49826;
 y[143][0] = a;

 x[144][0] = 0.001;
 x[144][1] = 3.1887;
 x[144][2] = -3.4143;
 x[144][3] = 2.7742;
 x[144][4] = -0.2026;
 y[144][0] = a;

 x[145][0] = 0.001;
 x[145][1] = 2.4527;
 x[145][2] = 2.9653;
 x[145][3] = 0.20021;
 x[145][4] = -0.056479;
 y[145][0] = a;

 x[146][0] = 0.001;
 x[146][1] = 3.9121;
 x[146][2] = 2.9735;
 x[146][3] = 0.92852;
 x[146][4] = 0.60558;
 y[146][0] = a;

 x[147][0] = 0.001;
 x[147][1] = 3.9364;
 x[147][2] = 10.5885;
 x[147][3] = -3.725;
 x[147][4] = -4.3133;
 y[147][0] = a;

 x[148][0] = 0.001;
 x[148][1] = 3.9414;
 x[148][2] = -3.2902;
 x[148][3] = 3.1674;
 x[148][4] = 1.0866;
 y[148][0] = a;

 x[149][0] = 0.001;
 x[149][1] = 3.6922;
 x[149][2] = -3.9585;
 x[149][3] = 4.3439;
 x[149][4] = 1.3517;
 y[149][0] = a;

 x[150][0] = 0.001;
 x[150][1] = 5.681;
 x[150][2] = 7.795;
 x[150][3] = -2.6848;
 x[150][4] = -0.92544;
 y[150][0] = a;

 x[151][0] = 0.001;
 x[151][1] = 0.77124;
 x[151][2] = 9.0862;
 x[151][3] = -1.2281;
 x[151][4] = -1.4996;
 y[151][0] = a;

 x[152][0] = 0.001;
 x[152][1] = 3.5761;
 x[152][2] = 9.7753;
 x[152][3] = -3.9795;
 x[152][4] = -3.4638;
 y[152][0] = a;

 x[153][0] = 0.001;
 x[153][1] = 1.602;
 x[153][2] = 6.1251;
 x[153][3] = 0.52924;
 x[153][4] = 0.47886;
 y[153][0] = a;

 x[154][0] = 0.001;
 x[154][1] = 2.6682;
 x[154][2] = 10.216;
 x[154][3] = -3.4414;
 x[154][4] = -4.0069;
 y[154][0] = a;

 x[155][0] = 0.001;
 x[155][1] = 2.0007;
 x[155][2] = 1.8644;
 x[155][3] = 2.6491;
 x[155][4] = 0.47369;
 y[155][0] = a;

 x[156][0] = 0.001;
 x[156][1] = 0.64215;
 x[156][2] = 3.1287;
 x[156][3] = 4.2933;
 x[156][4] = 0.64696;
 y[156][0] = a;

 x[157][0] = 0.001;
 x[157][1] = 4.3848;
 x[157][2] = -3.0729;
 x[157][3] = 3.0423;
 x[157][4] = 1.2741;
 y[157][0] = a;

 x[158][0] = 0.001;
 x[158][1] = 0.77445;
 x[158][2] = 9.0552;
 x[158][3] = -2.4089;
 x[158][4] = -1.3884;
 y[158][0] = a;

 x[159][0] = 0.001;
 x[159][1] = 0.96574;
 x[159][2] = 8.393;
 x[159][3] = -1.361;
 x[159][4] = -1.4659;
 y[159][0] = a;

 x[160][0] = 0.001;
 x[160][1] = 3.0948;
 x[160][2] = 8.7324;
 x[160][3] = -2.9007;
 x[160][4] = -0.96682;
 y[160][0] = a;

 x[161][0] = 0.001;
 x[161][1] = 4.9362;
 x[161][2] = 7.6046;
 x[161][3] = -2.3429;
 x[161][4] = -0.85302;
 y[161][0] = a;

 x[162][0] = 0.001;
 x[162][1] = -1.9458;
 x[162][2] = 11.2217;
 x[162][3] = 1.9079;
 x[162][4] = -3.4405;
 y[162][0] = a;

 x[163][0] = 0.001;
 x[163][1] = 5.7403;
 x[163][2] = -0.44284;
 x[163][3] = 0.38015;
 x[163][4] = 1.3763;
 y[163][0] = a;

 x[164][0] = 0.001;
 x[164][1] = -2.6989;
 x[164][2] = 12.1984;
 x[164][3] = 0.67661;
 x[164][4] = -8.5482;
 y[164][0] = a;

 x[165][0] = 0.001;
 x[165][1] = 1.1472;
 x[165][2] = 3.5985;
 x[165][3] = 1.9387;
 x[165][4] = -0.43406;
 y[165][0] = a;

 x[166][0] = 0.001;
 x[166][1] = 2.9742;
 x[166][2] = 8.96;
 x[166][3] = -2.9024;
 x[166][4] = -1.0379;
 y[166][0] = a;

 x[167][0] = 0.001;
 x[167][1] = 4.5707;
 x[167][2] = 7.2094;
 x[167][3] = -3.2794;
 x[167][4] = -1.4944;
 y[167][0] = a;

 x[168][0] = 0.001;
 x[168][1] = 0.1848;
 x[168][2] = 6.5079;
 x[168][3] = 2.0133;
 x[168][4] = -0.87242;
 y[168][0] = a;

 x[169][0] = 0.001;
 x[169][1] = 0.87256;
 x[169][2] = 9.2931;
 x[169][3] = -0.7843;
 x[169][4] = -2.1978;
 y[169][0] = a;

 x[170][0] = 0.001;
 x[170][1] = 0.39559;
 x[170][2] = 6.8866;
 x[170][3] = 1.0588;
 x[170][4] = -0.67587;
 y[170][0] = a;

 x[171][0] = 0.001;
 x[171][1] = 3.8384;
 x[171][2] = 6.1851;
 x[171][3] = -2.0439;
 x[171][4] = -0.033204;
 y[171][0] = a;

 x[172][0] = 0.001;
 x[172][1] = 2.8209;
 x[172][2] = 7.3108;
 x[172][3] = -0.81857;
 x[172][4] = -1.8784;
 y[172][0] = a;

 x[173][0] = 0.001;
 x[173][1] = 2.5817;
 x[173][2] = 9.7546;
 x[173][3] = -3.1749;
 x[173][4] = -2.9957;
 y[173][0] = a;

 x[174][0] = 0.001;
 x[174][1] = 3.8213;
 x[174][2] = 0.23175;
 x[174][3] = 2.0133;
 x[174][4] = 2.0564;
 y[174][0] = a;

 x[175][0] = 0.001;
 x[175][1] = 0.3798;
 x[175][2] = 0.7098;
 x[175][3] = 0.7572;
 x[175][4] = -0.4444;
 y[175][0] = a;

 x[176][0] = 0.001;
 x[176][1] = 3.4893;
 x[176][2] = 6.69;
 x[176][3] = -1.2042;
 x[176][4] = -0.38751;
 y[176][0] = a;

 x[177][0] = 0.001;
 x[177][1] = -1.7781;
 x[177][2] = 0.8546;
 x[177][3] = 7.1303;
 x[177][4] = 0.027572;
 y[177][0] = a;

 x[178][0] = 0.001;
 x[178][1] = 2.0962;
 x[178][2] = 2.4769;
 x[178][3] = 1.9379;
 x[178][4] = -0.040962;
 y[178][0] = a;

 x[179][0] = 0.001;
 x[179][1] = 0.94732;
 x[179][2] = -0.57113;
 x[179][3] = 7.1903;
 x[179][4] = -0.67587;
 y[179][0] = a;

 x[180][0] = 0.001;
 x[180][1] = 2.8261;
 x[180][2] = 9.4007;
 x[180][3] = -3.3034;
 x[180][4] = -1.0509;
 y[180][0] = a;

 x[181][0] = 0.001;
 x[181][1] = 0.007125;
 x[181][2] = 8.3661;
 x[181][3] = 0.50781;
 x[181][4] = -3.8155;
 y[181][0] = a;

 x[182][0] = 0.001;
 x[182][1] = 0.96788;
 x[182][2] = 7.1907;
 x[182][3] = 1.2798;
 x[182][4] = -2.4565;
 y[182][0] = a;

 x[183][0] = 0.001;
 x[183][1] = 4.7432;
 x[183][2] = 2.1086;
 x[183][3] = 0.1368;
 x[183][4] = 1.6543;
 y[183][0] = a;

 x[184][0] = 0.001;
 x[184][1] = 3.6575;
 x[184][2] = 7.2797;
 x[184][3] = -2.2692;
 x[184][4] = -1.144;
 y[184][0] = a;

 x[185][0] = 0.001;
 x[185][1] = 3.8832;
 x[185][2] = 6.4023;
 x[185][3] = -2.432;
 x[185][4] = -0.98363;
 y[185][0] = a;

 x[186][0] = 0.001;
 x[186][1] = 3.4776;
 x[186][2] = 8.811;
 x[186][3] = -3.1886;
 x[186][4] = -0.92285;
 y[186][0] = a;

 x[187][0] = 0.001;
 x[187][1] = 1.1315;
 x[187][2] = 7.9212;
 x[187][3] = 1.093;
 x[187][4] = -2.8444;
 y[187][0] = a;

 x[188][0] = 0.001;
 x[188][1] = 2.8237;
 x[188][2] = 2.8597;
 x[188][3] = 0.19678;
 x[188][4] = 0.57196;
 y[188][0] = a;

 x[189][0] = 0.001;
 x[189][1] = 1.9321;
 x[189][2] = 6.0423;
 x[189][3] = 0.26019;
 x[189][4] = -2.053;
 y[189][0] = a;

 x[190][0] = 0.001;
 x[190][1] = 3.0632;
 x[190][2] = -3.3315;
 x[190][3] = 5.1305;
 x[190][4] = 0.8267;
 y[190][0] = a;

 x[191][0] = 0.001;
 x[191][1] = -1.8411;
 x[191][2] = 10.8306;
 x[191][3] = 2.769;
 x[191][4] = -3.0901;
 y[191][0] = a;

 x[192][0] = 0.001;
 x[192][1] = 2.8084;
 x[192][2] = 11.3045;
 x[192][3] = -3.3394;
 x[192][4] = -4.4194;
 y[192][0] = a;

 x[193][0] = 0.001;
 x[193][1] = 2.5698;
 x[193][2] = -4.4076;
 x[193][3] = 5.9856;
 x[193][4] = 0.078002;
 y[193][0] = a;

 x[194][0] = 0.001;
 x[194][1] = -0.12624;
 x[194][2] = 10.3216;
 x[194][3] = -3.7121;
 x[194][4] = -6.1185;
 y[194][0] = a;

 x[195][0] = 0.001;
 x[195][1] = 3.3756;
 x[195][2] = -4.0951;
 x[195][3] = 4.367;
 x[195][4] = 1.0698;
 y[195][0] = a;

 x[196][0] = 0.001;
 x[196][1] = -0.048008;
 x[196][2] = -1.6037;
 x[196][3] = 8.4756;
 x[196][4] = 0.75558;
 y[196][0] = a;

 x[197][0] = 0.001;
 x[197][1] = 0.5706;
 x[197][2] = -0.0248;
 x[197][3] = 1.2421;
 x[197][4] = -0.5621;
 y[197][0] = a;

 x[198][0] = 0.001;
 x[198][1] = 0.88444;
 x[198][2] = 6.5906;
 x[198][3] = 0.55837;
 x[198][4] = -0.44182;
 y[198][0] = a;

 x[199][0] = 0.001;
 x[199][1] = 3.8644;
 x[199][2] = 3.7061;
 x[199][3] = 0.70403;
 x[199][4] = 0.35214;
 y[199][0] = a;

 x[200][0] = 0.001;
 x[200][1] = 1.2999;
 x[200][2] = 2.5762;
 x[200][3] = 2.0107;
 x[200][4] = -0.18967;
 y[200][0] = a;

 x[201][0] = 0.001;
 x[201][1] = 2.0051;
 x[201][2] = -6.8638;
 x[201][3] = 8.132;
 x[201][4] = -0.2401;
 y[201][0] = a;

 x[202][0] = 0.001;
 x[202][1] = 4.9294;
 x[202][2] = 0.27727;
 x[202][3] = 0.20792;
 x[202][4] = 0.33662;
 y[202][0] = a;

 x[203][0] = 0.001;
 x[203][1] = 2.8297;
 x[203][2] = 6.3485;
 x[203][3] = -0.73546;
 x[203][4] = -0.58665;
 y[203][0] = a;

 x[204][0] = 0.001;
 x[204][1] = 2.565;
 x[204][2] = 8.633;
 x[204][3] = -2.9941;
 x[204][4] = -1.3082;
 y[204][0] = a;

 x[205][0] = 0.001;
 x[205][1] = 2.093;
 x[205][2] = 8.3061;
 x[205][3] = 0.022844;
 x[205][4] = -3.2724;
 y[205][0] = a;

 x[206][0] = 0.001;
 x[206][1] = 4.6014;
 x[206][2] = 5.6264;
 x[206][3] = -2.1235;
 x[206][4] = 0.19309;
 y[206][0] = a;

 x[207][0] = 0.001;
 x[207][1] = 5.0617;
 x[207][2] = -0.35799;
 x[207][3] = 0.44698;
 x[207][4] = 0.99868;
 y[207][0] = a;

 x[208][0] = 0.001;
 x[208][1] = -0.2951;
 x[208][2] = 9.0489;
 x[208][3] = -0.52725;
 x[208][4] = -2.0789;
 y[208][0] = a;

 x[209][0] = 0.001;
 x[209][1] = 3.577;
 x[209][2] = 2.4004;
 x[209][3] = 1.8908;
 x[209][4] = 0.73231;
 y[209][0] = a;

 x[210][0] = 0.001;
 x[210][1] = 3.9433;
 x[210][2] = 2.5017;
 x[210][3] = 1.5215;
 x[210][4] = 0.903;
 y[210][0] = a;

 x[211][0] = 0.001;
 x[211][1] = 2.6648;
 x[211][2] = 10.754;
 x[211][3] = -3.3994;
 x[211][4] = -4.1685;
 y[211][0] = a;

 x[212][0] = 0.001;
 x[212][1] = 5.9374;
 x[212][2] = 6.1664;
 x[212][3] = -2.5905;
 x[212][4] = -0.36553;
 y[212][0] = a;

 x[213][0] = 0.001;
 x[213][1] = 2.0153;
 x[213][2] = 1.8479;
 x[213][3] = 3.1375;
 x[213][4] = 0.42843;
 y[213][0] = a;

 x[214][0] = 0.001;
 x[214][1] = 5.8782;
 x[214][2] = 5.9409;
 x[214][3] = -2.8544;
 x[214][4] = -0.60863;
 y[214][0] = a;

 x[215][0] = 0.001;
 x[215][1] = -2.3983;
 x[215][2] = 12.606;
 x[215][3] = 2.9464;
 x[215][4] = -5.7888;
 y[215][0] = a;

 x[216][0] = 0.001;
 x[216][1] = 1.762;
 x[216][2] = 4.3682;
 x[216][3] = 2.1384;
 x[216][4] = 0.75429;
 y[216][0] = a;

 x[217][0] = 0.001;
 x[217][1] = 4.2406;
 x[217][2] = -2.4852;
 x[217][3] = 1.608;
 x[217][4] = 0.7155;
 y[217][0] = a;

 x[218][0] = 0.001;
 x[218][1] = 3.4669;
 x[218][2] = 6.87;
 x[218][3] = -1.0568;
 x[218][4] = -0.73147;
 y[218][0] = a;

 x[219][0] = 0.001;
 x[219][1] = 3.1896;
 x[219][2] = 5.7526;
 x[219][3] = -0.18537;
 x[219][4] = -0.30087;
 y[219][0] = a;

 x[220][0] = 0.001;
 x[220][1] = 0.81356;
 x[220][2] = 9.1566;
 x[220][3] = -2.1492;
 x[220][4] = -4.1814;
 y[220][0] = a;

 x[221][0] = 0.001;
 x[221][1] = 0.52855;
 x[221][2] = 0.96427;
 x[221][3] = 4.0243;
 x[221][4] = -1.0483;
 y[221][0] = a;

 x[222][0] = 0.001;
 x[222][1] = 2.1319;
 x[222][2] = -2.0403;
 x[222][3] = 2.5574;
 x[222][4] = -0.061652;
 y[222][0] = a;

 x[223][0] = 0.001;
 x[223][1] = 0.33111;
 x[223][2] = 4.5731;
 x[223][3] = 2.057;
 x[223][4] = -0.18967;
 y[223][0] = a;

 x[224][0] = 0.001;
 x[224][1] = 1.2746;
 x[224][2] = 8.8172;
 x[224][3] = -1.5323;
 x[224][4] = -1.7957;
 y[224][0] = a;

 x[225][0] = 0.001;
 x[225][1] = 2.2091;
 x[225][2] = 7.4556;
 x[225][3] = -1.3284;
 x[225][4] = -3.3021;
 y[225][0] = a;

 x[226][0] = 0.001;
 x[226][1] = 2.5328;
 x[226][2] = 7.528;
 x[226][3] = -0.41929;
 x[226][4] = -2.6478;
 y[226][0] = a;

 x[227][0] = 0.001;
 x[227][1] = 3.6244;
 x[227][2] = 1.4609;
 x[227][3] = 1.3501;
 x[227][4] = 1.9284;
 y[227][0] = a;

 x[228][0] = 0.001;
 x[228][1] = -1.3885;
 x[228][2] = 12.5026;
 x[228][3] = 0.69118;
 x[228][4] = -7.5487;
 y[228][0] = a;

 x[229][0] = 0.001;
 x[229][1] = 5.7227;
 x[229][2] = 5.8312;
 x[229][3] = -2.4097;
 x[229][4] = -0.24527;
 y[229][0] = a;

 x[230][0] = 0.001;
 x[230][1] = 3.3583;
 x[230][2] = 10.3567;
 x[230][3] = -3.7301;
 x[230][4] = -3.6991;
 y[230][0] = a;

 x[231][0] = 0.001;
 x[231][1] = 2.5227;
 x[231][2] = 2.2369;
 x[231][3] = 2.7236;
 x[231][4] = 0.79438;
 y[231][0] = a;

 x[232][0] = 0.001;
 x[232][1] = 0.045304;
 x[232][2] = 6.7334;
 x[232][3] = 1.0708;
 x[232][4] = -0.9332;
 y[232][0] = a;

 x[233][0] = 0.001;
 x[233][1] = 4.8278;
 x[233][2] = 7.7598;
 x[233][3] = -2.4491;
 x[233][4] = -1.2216;
 y[233][0] = a;

 x[234][0] = 0.001;
 x[234][1] = 1.9476;
 x[234][2] = -4.7738;
 x[234][3] = 8.527;
 x[234][4] = -1.8668;
 y[234][0] = a;

 x[235][0] = 0.001;
 x[235][1] = 2.7659;
 x[235][2] = 0.66216;
 x[235][3] = 4.1494;
 x[235][4] = -0.28406;
 y[235][0] = a;

 x[236][0] = 0.001;
 x[236][1] = -0.90784;
 x[236][2] = -7.9026;
 x[236][3] = 6.7807;
 x[236][4] = 0.34179;
 y[236][0] = b;

 x[237][0] = 0.001;
 x[237][1] = -2.0042;
 x[237][2] = -9.3676;
 x[237][3] = 9.3333;
 x[237][4] = -0.10303;
 y[237][0] = b;

 x[238][0] = 0.001;
 x[238][1] = -0.93587;
 x[238][2] = -5.1008;
 x[238][3] = 4.5367;
 x[238][4] = 1.3866;
 y[238][0] = b;

 x[239][0] = 0.001;
 x[239][1] = -0.40804;
 x[239][2] = 0.54214;
 x[239][3] = -0.52725;
 x[239][4] = 0.6586;
 y[239][0] = b;

 x[240][0] = 0.001;
 x[240][1] = -0.8172;
 x[240][2] = 3.3812;
 x[240][3] = -3.6684;
 x[240][4] = -3.456;
 y[240][0] = b;

 x[241][0] = 0.001;
 x[241][1] = -4.8392;
 x[241][2] = 6.6755;
 x[241][3] = -0.24278;
 x[241][4] = -6.5775;
 y[241][0] = b;

 x[242][0] = 0.001;
 x[242][1] = -1.2792;
 x[242][2] = 2.1376;
 x[242][3] = -0.47584;
 x[242][4] = -1.3974;
 y[242][0] = b;

 x[243][0] = 0.001;
 x[243][1] = -0.66008;
 x[243][2] = -3.226;
 x[243][3] = 3.8058;
 x[243][4] = 1.1836;
 y[243][0] = b;

 x[244][0] = 0.001;
 x[244][1] = -1.7713;
 x[244][2] = -10.7665;
 x[244][3] = 10.2184;
 x[244][4] = -1.0043;
 y[244][0] = b;

 x[245][0] = 0.001;
 x[245][1] = -3.0061;
 x[245][2] = -12.2377;
 x[245][3] = 11.9552;
 x[245][4] = -2.1603;
 y[245][0] = b;

 x[246][0] = 0.001;
 x[246][1] = -1.1022;
 x[246][2] = -5.8395;
 x[246][3] = 4.5641;
 x[246][4] = 0.68705;
 y[246][0] = b;

 x[247][0] = 0.001;
 x[247][1] = 0.11806;
 x[247][2] = 0.39108;
 x[247][3] = -0.98223;
 x[247][4] = 0.42843;
 y[247][0] = b;

 x[248][0] = 0.001;
 x[248][1] = 0.11686;
 x[248][2] = 3.735;
 x[248][3] = -4.4379;
 x[248][4] = -4.3741;
 y[248][0] = b;

 x[249][0] = 0.001;
 x[249][1] = -2.7264;
 x[249][2] = 3.9213;
 x[249][3] = -0.49212;
 x[249][4] = -3.6371;
 y[249][0] = b;

 x[250][0] = 0.001;
 x[250][1] = -1.2369;
 x[250][2] = -1.6906;
 x[250][3] = 2.518;
 x[250][4] = 0.51636;
 y[250][0] = b;

 x[251][0] = 0.001;
 x[251][1] = -1.8439;
 x[251][2] = -8.6475;
 x[251][3] = 7.6796;
 x[251][4] = -0.66682;
 y[251][0] = b;

 x[252][0] = 0.001;
 x[252][1] = -1.8554;
 x[252][2] = -9.6035;
 x[252][3] = 7.7764;
 x[252][4] = -0.97716;
 y[252][0] = b;

 x[253][0] = 0.001;
 x[253][1] = 0.16358;
 x[253][2] = -3.3584;
 x[253][3] = 1.3749;
 x[253][4] = 1.3569;
 y[253][0] = b;

 x[254][0] = 0.001;
 x[254][1] = 1.5077;
 x[254][2] = 1.9596;
 x[254][3] = -3.0584;
 x[254][4] = -0.12243;
 y[254][0] = b;

 x[255][0] = 0.001;
 x[255][1] = 0.67886;
 x[255][2] = 4.1199;
 x[255][3] = -4.569;
 x[255][4] = -4.1414;
 y[255][0] = b;

 x[256][0] = 0.001;
 x[256][1] = -3.9934;
 x[256][2] = 5.8333;
 x[256][3] = 0.54723;
 x[256][4] = -4.9379;
 y[256][0] = b;

 x[257][0] = 0.001;
 x[257][1] = -2.3898;
 x[257][2] = -0.78427;
 x[257][3] = 3.0141;
 x[257][4] = 0.76205;
 y[257][0] = b;

 x[258][0] = 0.001;
 x[258][1] = -1.7976;
 x[258][2] = -6.7686;
 x[258][3] = 6.6753;
 x[258][4] = 0.89912;
 y[258][0] = b;

 x[259][0] = 0.001;
 x[259][1] = -0.70867;
 x[259][2] = -5.5602;
 x[259][3] = 4.0483;
 x[259][4] = 0.903;
 y[259][0] = b;

 x[260][0] = 0.001;
 x[260][1] = 1.0194;
 x[260][2] = 1.1029;
 x[260][3] = -2.3;
 x[260][4] = 0.59395;
 y[260][0] = b;

 x[261][0] = 0.001;
 x[261][1] = 1.7875;
 x[261][2] = 4.78;
 x[261][3] = -5.1362;
 x[261][4] = -3.2362;
 y[261][0] = b;

 x[262][0] = 0.001;
 x[262][1] = 0.27331;
 x[262][2] = 4.8773;
 x[262][3] = -4.9194;
 x[262][4] = -5.8198;
 y[262][0] = b;

 x[263][0] = 0.001;
 x[263][1] = -5.1661;
 x[263][2] = 8.0433;
 x[263][3] = 0.044265;
 x[263][4] = -4.4983;
 y[263][0] = b;

 x[264][0] = 0.001;
 x[264][1] = -2.7028;
 x[264][2] = 1.6327;
 x[264][3] = 0.83598;
 x[264][4] = -0.091393;
 y[264][0] = b;

 x[265][0] = 0.001;
 x[265][1] = -1.4904;
 x[265][2] = -2.2183;
 x[265][3] = 1.6054;
 x[265][4] = 0.89394;
 y[265][0] = b;

 x[266][0] = 0.001;
 x[266][1] = -0.014902;
 x[266][2] = -1.0243;
 x[266][3] = -0.94024;
 x[266][4] = 0.64955;
 y[266][0] = b;

 x[267][0] = 0.001;
 x[267][1] = 0.88992;
 x[267][2] = 2.2638;
 x[267][3] = -3.1046;
 x[267][4] = -0.11855;
 y[267][0] = b;

 x[268][0] = 0.001;
 x[268][1] = 1.0637;
 x[268][2] = 3.6957;
 x[268][3] = -4.1594;
 x[268][4] = -1.9379;
 y[268][0] = b;

 x[269][0] = 0.001;
 x[269][1] = -0.8471;
 x[269][2] = 3.1329;
 x[269][3] = -3.0112;
 x[269][4] = -2.9388;
 y[269][0] = b;

 x[270][0] = 0.001;
 x[270][1] = -3.9594;
 x[270][2] = 4.0289;
 x[270][3] = -0.35845;
 x[270][4] = -3.8957;
 y[270][0] = b;

 x[271][0] = 0.001;
 x[271][1] = -5.8818;
 x[271][2] = 7.6584;
 x[271][3] = 0.5558;
 x[271][4] = -2.9155;
 y[271][0] = b;

 x[272][0] = 0.001;
 x[272][1] = -3.7747;
 x[272][2] = 2.5162;
 x[272][3] = 0.83341;
 x[272][4] = -0.30993;
 y[272][0] = b;

 x[273][0] = 0.001;
 x[273][1] = -2.4198;
 x[273][2] = -0.24418;
 x[273][3] = 0.70146;
 x[273][4] = 0.41809;
 y[273][0] = b;

 x[274][0] = 0.001;
 x[274][1] = -0.83535;
 x[274][2] = 0.80494;
 x[274][3] = -1.6411;
 x[274][4] = -0.19225;
 y[274][0] = b;

 x[275][0] = 0.001;
 x[275][1] = -0.30432;
 x[275][2] = 2.6528;
 x[275][3] = -2.7756;
 x[275][4] = -0.65647;
 y[275][0] = b;

 x[276][0] = 0.001;
 x[276][1] = -0.60254;
 x[276][2] = 1.7237;
 x[276][3] = -2.1501;
 x[276][4] = -0.77027;
 y[276][0] = b;

 x[277][0] = 0.001;
 x[277][1] = -2.1059;
 x[277][2] = 1.1815;
 x[277][3] = -0.53324;
 x[277][4] = -0.82716;
 y[277][0] = b;

 x[278][0] = 0.001;
 x[278][1] = -2.0441;
 x[278][2] = 1.2271;
 x[278][3] = 0.18564;
 x[278][4] = -1.091;
 y[278][0] = b;

 x[279][0] = 0.001;
 x[279][1] = -1.5621;
 x[279][2] = -2.2121;
 x[279][3] = 4.2591;
 x[279][4] = 0.27972;
 y[279][0] = b;

 x[280][0] = 0.001;
 x[280][1] = -3.2305;
 x[280][2] = -7.2135;
 x[280][3] = 11.6433;
 x[280][4] = -0.94613;
 y[280][0] = b;

 x[281][0] = 0.001;
 x[281][1] = -4.8426;
 x[281][2] = -4.9932;
 x[281][3] = 10.4052;
 x[281][4] = -0.53104;
 y[281][0] = b;

 x[282][0] = 0.001;
 x[282][1] = -2.3147;
 x[282][2] = 3.6668;
 x[282][3] = -0.6969;
 x[282][4] = -1.2474;
 y[282][0] = b;

 x[283][0] = 0.001;
 x[283][1] = -0.11716;
 x[283][2] = 0.60422;
 x[283][3] = -0.38587;
 x[283][4] = -0.059065;
 y[283][0] = b;

 x[284][0] = 0.001;
 x[284][1] = -2.0066;
 x[284][2] = -6.719;
 x[284][3] = 9.0162;
 x[284][4] = 0.099985;
 y[284][0] = b;

 x[285][0] = 0.001;
 x[285][1] = -3.6961;
 x[285][2] = -13.6779;
 x[285][3] = 17.5795;
 x[285][4] = -2.6181;
 y[285][0] = b;

 x[286][0] = 0.001;
 x[286][1] = -3.6012;
 x[286][2] = -6.5389;
 x[286][3] = 10.5234;
 x[286][4] = -0.48967;
 y[286][0] = b;

 x[287][0] = 0.001;
 x[287][1] = -2.6286;
 x[287][2] = 0.18002;
 x[287][3] = 1.7956;
 x[287][4] = 0.97282;
 y[287][0] = b;

 x[288][0] = 0.001;
 x[288][1] = -0.82601;
 x[288][2] = 2.9611;
 x[288][3] = -1.2864;
 x[288][4] = -1.4647;
 y[288][0] = b;

 x[289][0] = 0.001;
 x[289][1] = 0.31803;
 x[289][2] = -0.99326;
 x[289][3] = 1.0947;
 x[289][4] = 0.88619;
 y[289][0] = b;

 x[290][0] = 0.001;
 x[290][1] = -1.4454;
 x[290][2] = -8.4385;
 x[290][3] = 8.8483;
 x[290][4] = 0.96894;
 y[290][0] = b;

 x[291][0] = 0.001;
 x[291][1] = -3.1423;
 x[291][2] = -13.0365;
 x[291][3] = 15.6773;
 x[291][4] = -0.66165;
 y[291][0] = b;

 x[292][0] = 0.001;
 x[292][1] = -2.5373;
 x[292][2] = -6.959;
 x[292][3] = 8.8054;
 x[292][4] = 1.5289;
 y[292][0] = b;

 x[293][0] = 0.001;
 x[293][1] = -1.366;
 x[293][2] = 0.18416;
 x[293][3] = 0.90539;
 x[293][4] = 1.5806;
 y[293][0] = b;

 x[294][0] = 0.001;
 x[294][1] = -1.7064;
 x[294][2] = 3.3088;
 x[294][3] = -2.2829;
 x[294][4] = -2.1978;
 y[294][0] = b;

 x[295][0] = 0.001;
 x[295][1] = -0.41965;
 x[295][2] = 2.9094;
 x[295][3] = -1.7859;
 x[295][4] = -2.2069;
 y[295][0] = b;

 x[296][0] = 0.001;
 x[296][1] = 0.37637;
 x[296][2] = -0.82358;
 x[296][3] = 0.78543;
 x[296][4] = 0.74524;
 y[296][0] = b;

 x[297][0] = 0.001;
 x[297][1] = -0.55355;
 x[297][2] = -7.9233;
 x[297][3] = 6.7156;
 x[297][4] = 0.74394;
 y[297][0] = b;

 x[298][0] = 0.001;
 x[298][1] = -1.6001;
 x[298][2] = -9.5828;
 x[298][3] = 9.4044;
 x[298][4] = 0.081882;
 y[298][0] = b;

 x[299][0] = 0.001;
 x[299][1] = -0.37013;
 x[299][2] = -5.554;
 x[299][3] = 4.7749;
 x[299][4] = 1.547;
 y[299][0] = b;

 x[300][0] = 0.001;
 x[300][1] = 0.12126;
 x[300][2] = 0.22347;
 x[300][3] = -0.47327;
 x[300][4] = 0.97024;
 y[300][0] = b;

 x[301][0] = 0.001;
 x[301][1] = -0.27068;
 x[301][2] = 3.2674;
 x[301][3] = -3.5562;
 x[301][4] = -3.0888;
 y[301][0] = b;

 x[302][0] = 0.001;
 x[302][1] = -5.119;
 x[302][2] = 6.6486;
 x[302][3] = -0.049987;
 x[302][4] = -6.5206;
 y[302][0] = b;

 x[303][0] = 0.001;
 x[303][1] = -1.3946;
 x[303][2] = 2.3134;
 x[303][3] = -0.44499;
 x[303][4] = -1.4905;
 y[303][0] = b;

 x[304][0] = 0.001;
 x[304][1] = -0.69879;
 x[304][2] = -3.3771;
 x[304][3] = 4.1211;
 x[304][4] = 1.5043;
 y[304][0] = b;

 x[305][0] = 0.001;
 x[305][1] = -1.48;
 x[305][2] = -10.5244;
 x[305][3] = 9.9176;
 x[305][4] = -0.5026;
 y[305][0] = b;

 x[306][0] = 0.001;
 x[306][1] = -2.6649;
 x[306][2] = -12.813;
 x[306][3] = 12.6689;
 x[306][4] = -1.9082;
 y[306][0] = b;

 x[307][0] = 0.001;
 x[307][1] = -0.62684;
 x[307][2] = -6.301;
 x[307][3] = 4.7843;
 x[307][4] = 1.106;
 y[307][0] = b;

 x[308][0] = 0.001;
 x[308][1] = 0.518;
 x[308][2] = 0.25865;
 x[308][3] = -0.84085;
 x[308][4] = 0.96118;
 y[308][0] = b;

 x[309][0] = 0.001;
 x[309][1] = 0.64376;
 x[309][2] = 3.764;
 x[309][3] = -4.4738;
 x[309][4] = -4.0483;
 y[309][0] = b;

 x[310][0] = 0.001;
 x[310][1] = -2.9821;
 x[310][2] = 4.1986;
 x[310][3] = -0.5898;
 x[310][4] = -3.9642;
 y[310][0] = b;

 x[311][0] = 0.001;
 x[311][1] = -1.4628;
 x[311][2] = -1.5706;
 x[311][3] = 2.4357;
 x[311][4] = 0.49826;
 y[311][0] = b;

 x[312][0] = 0.001;
 x[312][1] = -1.7101;
 x[312][2] = -8.7903;
 x[312][3] = 7.9735;
 x[312][4] = -0.45475;
 y[312][0] = b;

 x[313][0] = 0.001;
 x[313][1] = -1.5572;
 x[313][2] = -9.8808;
 x[313][3] = 8.1088;
 x[313][4] = -1.0806;
 y[313][0] = b;

 x[314][0] = 0.001;
 x[314][1] = 0.74428;
 x[314][2] = -3.7723;
 x[314][3] = 1.6131;
 x[314][4] = 1.5754;
 y[314][0] = b;

 x[315][0] = 0.001;
 x[315][1] = 2.0177;
 x[315][2] = 1.7982;
 x[315][3] = -2.9581;
 x[315][4] = 0.2099;
 y[315][0] = b;

 x[316][0] = 0.001;
 x[316][1] = 1.164;
 x[316][2] = 3.913;
 x[316][3] = -4.5544;
 x[316][4] = -3.8672;
 y[316][0] = b;

 x[317][0] = 0.001;
 x[317][1] = -4.3667;
 x[317][2] = 6.0692;
 x[317][3] = 0.57208;
 x[317][4] = -5.4668;
 y[317][0] = b;

 x[318][0] = 0.001;
 x[318][1] = -2.5919;
 x[318][2] = -1.0553;
 x[318][3] = 3.8949;
 x[318][4] = 0.77757;
 y[318][0] = b;

 x[319][0] = 0.001;
 x[319][1] = -1.8046;
 x[319][2] = -6.8141;
 x[319][3] = 6.7019;
 x[319][4] = 1.1681;
 y[319][0] = b;

 x[320][0] = 0.001;
 x[320][1] = -0.71868;
 x[320][2] = -5.7154;
 x[320][3] = 3.8298;
 x[320][4] = 1.0233;
 y[320][0] = b;

 x[321][0] = 0.001;
 x[321][1] = 1.4378;
 x[321][2] = 0.66837;
 x[321][3] = -2.0267;
 x[321][4] = 1.0271;
 y[321][0] = b;

 x[322][0] = 0.001;
 x[322][1] = 2.1943;
 x[322][2] = 4.5503;
 x[322][3] = -4.976;
 x[322][4] = -2.7254;
 y[322][0] = b;

 x[323][0] = 0.001;
 x[323][1] = 0.7376;
 x[323][2] = 4.8525;
 x[323][3] = -4.7986;
 x[323][4] = -5.6659;
 y[323][0] = b;

 x[324][0] = 0.001;
 x[324][1] = -5.637;
 x[324][2] = 8.1261;
 x[324][3] = 0.13081;
 x[324][4] = -5.0142;
 y[324][0] = b;

 x[325][0] = 0.001;
 x[325][1] = -3.0193;
 x[325][2] = 1.7775;
 x[325][3] = 0.73745;
 x[325][4] = -0.45346;
 y[325][0] = b;

 x[326][0] = 0.001;
 x[326][1] = -1.6706;
 x[326][2] = -2.09;
 x[326][3] = 1.584;
 x[326][4] = 0.71162;
 y[326][0] = b;

 x[327][0] = 0.001;
 x[327][1] = -0.1269;
 x[327][2] = -1.1505;
 x[327][3] = -0.95138;
 x[327][4] = 0.57843;
 y[327][0] = b;

 x[328][0] = 0.001;
 x[328][1] = 1.2198;
 x[328][2] = 2.0982;
 x[328][3] = -3.1954;
 x[328][4] = 0.12843;
 y[328][0] = b;

 x[329][0] = 0.001;
 x[329][1] = 1.4501;
 x[329][2] = 3.6067;
 x[329][3] = -4.0557;
 x[329][4] = -1.5966;
 y[329][0] = b;

 x[330][0] = 0.001;
 x[330][1] = -0.40857;
 x[330][2] = 3.0977;
 x[330][3] = -2.9607;
 x[330][4] = -2.6892;
 y[330][0] = b;

 x[331][0] = 0.001;
 x[331][1] = -3.8952;
 x[331][2] = 3.8157;
 x[331][3] = -0.31304;
 x[331][4] = -3.8194;
 y[331][0] = b;

 x[332][0] = 0.001;
 x[332][1] = -6.3679;
 x[332][2] = 8.0102;
 x[332][3] = 0.4247;
 x[332][4] = -3.2207;
 y[332][0] = b;

 x[333][0] = 0.001;
 x[333][1] = -4.1429;
 x[333][2] = 2.7749;
 x[333][3] = 0.68261;
 x[333][4] = -0.71984;
 y[333][0] = b;

 x[334][0] = 0.001;
 x[334][1] = -2.6864;
 x[334][2] = -0.097265;
 x[334][3] = 0.61663;
 x[334][4] = 0.061192;
 y[334][0] = b;

 x[335][0] = 0.001;
 x[335][1] = -1.0555;
 x[335][2] = 0.79459;
 x[335][3] = -1.6968;
 x[335][4] = -0.46768;
 y[335][0] = b;

 x[336][0] = 0.001;
 x[336][1] = -0.29858;
 x[336][2] = 2.4769;
 x[336][3] = -2.9512;
 x[336][4] = -0.66165;
 y[336][0] = b;

 x[337][0] = 0.001;
 x[337][1] = -0.49948;
 x[337][2] = 1.7734;
 x[337][3] = -2.2469;
 x[337][4] = -0.68104;
 y[337][0] = b;

 x[338][0] = 0.001;
 x[338][1] = -1.9881;
 x[338][2] = 0.99945;
 x[338][3] = -0.28562;
 x[338][4] = -0.70044;
 y[338][0] = b;

 x[339][0] = 0.001;
 x[339][1] = -1.9389;
 x[339][2] = 1.5706;
 x[339][3] = 0.045979;
 x[339][4] = -1.122;
 y[339][0] = b;

 x[340][0] = 0.001;
 x[340][1] = -1.4375;
 x[340][2] = -1.8624;
 x[340][3] = 4.026;
 x[340][4] = 0.55127;
 y[340][0] = b;

 x[341][0] = 0.001;
 x[341][1] = -3.1875;
 x[341][2] = -7.5756;
 x[341][3] = 11.8678;
 x[341][4] = -0.57889;
 y[341][0] = b;

 x[342][0] = 0.001;
 x[342][1] = -4.6765;
 x[342][2] = -5.6636;
 x[342][3] = 10.969;
 x[342][4] = -0.33449;
 y[342][0] = b;

 x[343][0] = 0.001;
 x[343][1] = -2.0285;
 x[343][2] = 3.8468;
 x[343][3] = -0.63435;
 x[343][4] = -1.175;
 y[343][0] = b;

 x[344][0] = 0.001;
 x[344][1] = 0.26637;
 x[344][2] = 0.73252;
 x[344][3] = -0.67891;
 x[344][4] = 0.03533;
 y[344][0] = b;

 x[345][0] = 0.001;
 x[345][1] = -1.7589;
 x[345][2] = -6.4624;
 x[345][3] = 8.4773;
 x[345][4] = 0.31981;
 y[345][0] = b;

 x[346][0] = 0.001;
 x[346][1] = -3.5985;
 x[346][2] = -13.6593;
 x[346][3] = 17.6052;
 x[346][4] = -2.4927;
 y[346][0] = b;

 x[347][0] = 0.001;
 x[347][1] = -3.3582;
 x[347][2] = -7.2404;
 x[347][3] = 11.4419;
 x[347][4] = -0.57113;
 y[347][0] = b;

 x[348][0] = 0.001;
 x[348][1] = -2.3629;
 x[348][2] = -0.10554;
 x[348][3] = 1.9336;
 x[348][4] = 1.1358;
 y[348][0] = b;

 x[349][0] = 0.001;
 x[349][1] = -2.1802;
 x[349][2] = 3.3791;
 x[349][3] = -1.2256;
 x[349][4] = -2.6621;
 y[349][0] = b;

 x[350][0] = 0.001;
 x[350][1] = -0.40951;
 x[350][2] = -0.15521;
 x[350][3] = 0.060545;
 x[350][4] = -0.088807;
 y[350][0] = b;

 x[351][0] = 0.001;
 x[351][1] = -2.2918;
 x[351][2] = -7.257;
 x[351][3] = 7.9597;
 x[351][4] = 0.9211;
 y[351][0] = b;

 x[352][0] = 0.001;
 x[352][1] = -4.0214;
 x[352][2] = -12.8006;
 x[352][3] = 15.6199;
 x[352][4] = -0.95647;
 y[352][0] = b;

 x[353][0] = 0.001;
 x[353][1] = -3.3884;
 x[353][2] = -8.215;
 x[353][3] = 10.3315;
 x[353][4] = 0.98187;
 y[353][0] = b;

 x[354][0] = 0.001;
 x[354][1] = -2.0046;
 x[354][2] = -0.49457;
 x[354][3] = 1.333;
 x[354][4] = 1.6543;
 y[354][0] = b;

 x[355][0] = 0.001;
 x[355][1] = -1.7063;
 x[355][2] = 2.7956;
 x[355][3] = -2.378;
 x[355][4] = -2.3491;
 y[355][0] = b;

 x[356][0] = 0.001;
 x[356][1] = -1.6386;
 x[356][2] = 3.3584;
 x[356][3] = -1.7302;
 x[356][4] = -3.5646;
 y[356][0] = b;

 x[357][0] = 0.001;
 x[357][1] = -0.41645;
 x[357][2] = 0.32487;
 x[357][3] = -0.33617;
 x[357][4] = -0.36036;
 y[357][0] = b;

 x[358][0] = 0.001;
 x[358][1] = -1.5877;
 x[358][2] = -6.6072;
 x[358][3] = 5.8022;
 x[358][4] = 0.31593;
 y[358][0] = b;

 x[359][0] = 0.001;
 x[359][1] = -2.5961;
 x[359][2] = -9.349;
 x[359][3] = 9.7942;
 x[359][4] = -0.28018;
 y[359][0] = b;

 x[360][0] = 0.001;
 x[360][1] = -1.5228;
 x[360][2] = -6.4789;
 x[360][3] = 5.7568;
 x[360][4] = 0.87325;
 y[360][0] = b;

 x[361][0] = 0.001;
 x[361][1] = -0.53072;
 x[361][2] = -0.097265;
 x[361][3] = -0.21793;
 x[361][4] = 1.0426;
 y[361][0] = b;

 x[362][0] = 0.001;
 x[362][1] = -0.49081;
 x[362][2] = 2.8452;
 x[362][3] = -3.6436;
 x[362][4] = -3.1004;
 y[362][0] = b;

 x[363][0] = 0.001;
 x[363][1] = -6.5773;
 x[363][2] = 6.8017;
 x[363][3] = 0.85483;
 x[363][4] = -7.5344;
 y[363][0] = b;

 x[364][0] = 0.001;
 x[364][1] = -2.4621;
 x[364][2] = 2.7645;
 x[364][3] = -0.62578;
 x[364][4] = -2.8573;
 y[364][0] = b;

 x[365][0] = 0.001;
 x[365][1] = -1.3995;
 x[365][2] = -1.9162;
 x[365][3] = 2.5154;
 x[365][4] = 0.59912;
 y[365][0] = b;

 x[366][0] = 0.001;
 x[366][1] = -2.3221;
 x[366][2] = -9.3304;
 x[366][3] = 9.233;
 x[366][4] = -0.79871;
 y[366][0] = b;

 x[367][0] = 0.001;
 x[367][1] = -3.73;
 x[367][2] = -12.9723;
 x[367][3] = 12.9817;
 x[367][4] = -2.684;
 y[367][0] = b;

 x[368][0] = 0.001;
 x[368][1] = -1.6988;
 x[368][2] = -7.1163;
 x[368][3] = 5.7902;
 x[368][4] = 0.16723;
 y[368][0] = b;

 x[369][0] = 0.001;
 x[369][1] = -0.26654;
 x[369][2] = -0.64562;
 x[369][3] = -0.42014;
 x[369][4] = 0.89136;
 y[369][0] = b;

 x[370][0] = 0.001;
 x[370][1] = 0.33325;
 x[370][2] = 3.3108;
 x[370][3] = -4.5081;
 x[370][4] = -4.012;
 y[370][0] = b;

 x[371][0] = 0.001;
 x[371][1] = -4.2091;
 x[371][2] = 4.7283;
 x[371][3] = -0.49126;
 x[371][4] = -5.2159;
 y[371][0] = b;

 x[372][0] = 0.001;
 x[372][1] = -2.3142;
 x[372][2] = -0.68494;
 x[372][3] = 1.9833;
 x[372][4] = -0.44829;
 y[372][0] = b;

 x[373][0] = 0.001;
 x[373][1] = -2.4835;
 x[373][2] = -7.4494;
 x[373][3] = 6.8964;
 x[373][4] = -0.64484;
 y[373][0] = b;

 x[374][0] = 0.001;
 x[374][1] = -2.7611;
 x[374][2] = -10.5099;
 x[374][3] = 9.0239;
 x[374][4] = -1.9547;
 y[374][0] = b;

 x[375][0] = 0.001;
 x[375][1] = -0.36025;
 x[375][2] = -4.449;
 x[375][3] = 2.1067;
 x[375][4] = 0.94308;
 y[375][0] = b;

 x[376][0] = 0.001;
 x[376][1] = 1.0117;
 x[376][2] = 0.9022;
 x[376][3] = -2.3506;
 x[376][4] = 0.42714;
 y[376][0] = b;

 x[377][0] = 0.001;
 x[377][1] = 0.96708;
 x[377][2] = 3.8426;
 x[377][3] = -4.9314;
 x[377][4] = -4.1323;
 y[377][0] = b;

 x[378][0] = 0.001;
 x[378][1] = -5.2049;
 x[378][2] = 7.259;
 x[378][3] = 0.070827;
 x[378][4] = -7.3004;
 y[378][0] = b;

 x[379][0] = 0.001;
 x[379][1] = -3.3203;
 x[379][2] = -0.02691;
 x[379][3] = 2.9618;
 x[379][4] = -0.44958;
 y[379][0] = b;

 x[380][0] = 0.001;
 x[380][1] = -2.565;
 x[380][2] = -5.7899;
 x[380][3] = 6.0122;
 x[380][4] = 0.046968;
 y[380][0] = b;

 x[381][0] = 0.001;
 x[381][1] = -1.5951;
 x[381][2] = -6.572;
 x[381][3] = 4.7689;
 x[381][4] = -0.94354;
 y[381][0] = b;

 x[382][0] = 0.001;
 x[382][1] = 0.7049;
 x[382][2] = 0.17174;
 x[382][3] = -1.7859;
 x[382][4] = 0.36119;
 y[382][0] = b;

 x[383][0] = 0.001;
 x[383][1] = 1.7331;
 x[383][2] = 3.9544;
 x[383][3] = -4.7412;
 x[383][4] = -2.5017;
 y[383][0] = b;

 x[384][0] = 0.001;
 x[384][1] = 0.6818;
 x[384][2] = 4.8504;
 x[384][3] = -5.2133;
 x[384][4] = -6.1043;
 y[384][0] = b;

 x[385][0] = 0.001;
 x[385][1] = -6.3364;
 x[385][2] = 9.2848;
 x[385][3] = 0.014275;
 x[385][4] = -6.7844;
 y[385][0] = b;

 x[386][0] = 0.001;
 x[386][1] = -3.8053;
 x[386][2] = 2.4273;
 x[386][3] = 0.6809;
 x[386][4] = -1.0871;
 y[386][0] = b;

 x[387][0] = 0.001;
 x[387][1] = -2.1979;
 x[387][2] = -2.1252;
 x[387][3] = 1.7151;
 x[387][4] = 0.45171;
 y[387][0] = b;

 x[388][0] = 0.001;
 x[388][1] = -0.87874;
 x[388][2] = -2.2121;
 x[388][3] = -0.051701;
 x[388][4] = 0.099985;
 y[388][0] = b;

 x[389][0] = 0.001;
 x[389][1] = 0.74067;
 x[389][2] = 1.7299;
 x[389][3] = -3.1963;
 x[389][4] = -0.1457;
 y[389][0] = b;

 x[390][0] = 0.001;
 x[390][1] = 0.98296;
 x[390][2] = 3.4226;
 x[390][3] = -3.9692;
 x[390][4] = -1.7116;
 y[390][0] = b;

 x[391][0] = 0.001;
 x[391][1] = -0.3489;
 x[391][2] = 3.1929;
 x[391][3] = -3.4054;
 x[391][4] = -3.1832;
 y[391][0] = b;

 x[392][0] = 0.001;
 x[392][1] = -3.8552;
 x[392][2] = 3.5219;
 x[392][3] = -0.38415;
 x[392][4] = -3.8608;
 y[392][0] = b;

 x[393][0] = 0.001;
 x[393][1] = -6.9599;
 x[393][2] = 8.9931;
 x[393][3] = 0.2182;
 x[393][4] = -4.572;
 y[393][0] = b;

 x[394][0] = 0.001;
 x[394][1] = -4.7462;
 x[394][2] = 3.1205;
 x[394][3] = 1.075;
 x[394][4] = -1.2966;
 y[394][0] = b;

 x[395][0] = 0.001;
 x[395][1] = -3.2051;
 x[395][2] = -0.14279;
 x[395][3] = 0.97565;
 x[395][4] = 0.045675;
 y[395][0] = b;

 x[396][0] = 0.001;
 x[396][1] = -1.7549;
 x[396][2] = -0.080711;
 x[396][3] = -0.75774;
 x[396][4] = -0.3707;
 y[396][0] = b;

 x[397][0] = 0.001;
 x[397][1] = -0.59587;
 x[397][2] = 2.4811;
 x[397][3] = -2.8673;
 x[397][4] = -0.89828;
 y[397][0] = b;

 x[398][0] = 0.001;
 x[398][1] = -0.89542;
 x[398][2] = 2.0279;
 x[398][3] = -2.3652;
 x[398][4] = -1.2746;
 y[398][0] = b;

 x[399][0] = 0.001;
 x[399][1] = -2.0754;
 x[399][2] = 1.2767;
 x[399][3] = -0.64206;
 x[399][4] = -1.2642;
 y[399][0] = b;

 x[400][0] = 0.001;
 x[400][1] = -3.2778;
 x[400][2] = 1.8023;
 x[400][3] = 0.1805;
 x[400][4] = -2.3931;
 y[400][0] = b;

 x[401][0] = 0.001;
 x[401][1] = -2.2183;
 x[401][2] = -1.254;
 x[401][3] = 2.9986;
 x[401][4] = 0.36378;
 y[401][0] = b;

 x[402][0] = 0.001;
 x[402][1] = -3.5895;
 x[402][2] = -6.572;
 x[402][3] = 10.5251;
 x[402][4] = -0.16381;
 y[402][0] = b;

 x[403][0] = 0.001;
 x[403][1] = -5.0477;
 x[403][2] = -5.8023;
 x[403][3] = 11.244;
 x[403][4] = -0.3901;
 y[403][0] = b;

 x[404][0] = 0.001;
 x[404][1] = -3.5741;
 x[404][2] = 3.944;
 x[404][3] = -0.07912;
 x[404][4] = -2.1203;
 y[404][0] = b;

 x[405][0] = 0.001;
 x[405][1] = -0.7351;
 x[405][2] = 1.7361;
 x[405][3] = -1.4938;
 x[405][4] = -1.1582;
 y[405][0] = b;

 x[406][0] = 0.001;
 x[406][1] = -2.2617;
 x[406][2] = -4.7428;
 x[406][3] = 6.3489;
 x[406][4] = 0.11162;
 y[406][0] = b;

 x[407][0] = 0.001;
 x[407][1] = -4.244;
 x[407][2] = -13.0634;
 x[407][3] = 17.1116;
 x[407][4] = -2.8017;
 y[407][0] = b;

 x[408][0] = 0.001;
 x[408][1] = -4.0218;
 x[408][2] = -8.304;
 x[408][3] = 12.555;
 x[408][4] = -1.5099;
 y[408][0] = b;

 x[409][0] = 0.001;
 x[409][1] = -3.0201;
 x[409][2] = -0.67253;
 x[409][3] = 2.7056;
 x[409][4] = 0.85774;
 y[409][0] = b;

 x[410][0] = 0.001;
 x[410][1] = -2.4941;
 x[410][2] = 3.5447;
 x[410][3] = -1.3721;
 x[410][4] = -2.8483;
 y[410][0] = b;

 x[411][0] = 0.001;
 x[411][1] = -0.83121;
 x[411][2] = 0.039307;
 x[411][3] = 0.05369;
 x[411][4] = -0.23105;
 y[411][0] = b;

 x[412][0] = 0.001;
 x[412][1] = -2.5665;
 x[412][2] = -6.8824;
 x[412][3] = 7.5416;
 x[412][4] = 0.70774;
 y[412][0] = b;

 x[413][0] = 0.001;
 x[413][1] = -4.4018;
 x[413][2] = -12.9371;
 x[413][3] = 15.6559;
 x[413][4] = -1.6806;
 y[413][0] = b;

 x[414][0] = 0.001;
 x[414][1] = -3.7573;
 x[414][2] = -8.2916;
 x[414][3] = 10.3032;
 x[414][4] = 0.38059;
 y[414][0] = b;

 x[415][0] = 0.001;
 x[415][1] = -2.4725;
 x[415][2] = -0.40145;
 x[415][3] = 1.4855;
 x[415][4] = 1.1189;
 y[415][0] = b;

 x[416][0] = 0.001;
 x[416][1] = -1.9725;
 x[416][2] = 2.8825;
 x[416][3] = -2.3086;
 x[416][4] = -2.3724;
 y[416][0] = b;

 x[417][0] = 0.001;
 x[417][1] = -2.0149;
 x[417][2] = 3.6874;
 x[417][3] = -1.9385;
 x[417][4] = -3.8918;
 y[417][0] = b;

 x[418][0] = 0.001;
 x[418][1] = -0.82053;
 x[418][2] = 0.65181;
 x[418][3] = -0.48869;
 x[418][4] = -0.52716;
 y[418][0] = b;

 x[419][0] = 0.001;
 x[419][1] = -1.7886;
 x[419][2] = -6.3486;
 x[419][3] = 5.6154;
 x[419][4] = 0.42584;
 y[419][0] = b;

 x[420][0] = 0.001;
 x[420][1] = -2.9138;
 x[420][2] = -9.4711;
 x[420][3] = 9.7668;
 x[420][4] = -0.60216;
 y[420][0] = b;

 x[421][0] = 0.001;
 x[421][1] = -1.8343;
 x[421][2] = -6.5907;
 x[421][3] = 5.6429;
 x[421][4] = 0.54998;
 y[421][0] = b;

 x[422][0] = 0.001;
 x[422][1] = -0.8734;
 x[422][2] = -0.033118;
 x[422][3] = -0.20165;
 x[422][4] = 0.55774;
 y[422][0] = b;

 x[423][0] = 0.001;
 x[423][1] = -0.70346;
 x[423][2] = 2.957;
 x[423][3] = -3.5947;
 x[423][4] = -3.1457;
 y[423][0] = b;

 x[424][0] = 0.001;
 x[424][1] = -6.7387;
 x[424][2] = 6.9879;
 x[424][3] = 0.67833;
 x[424][4] = -7.5887;
 y[424][0] = b;

 x[425][0] = 0.001;
 x[425][1] = -2.7723;
 x[425][2] = 3.2777;
 x[425][3] = -0.9351;
 x[425][4] = -3.1457;
 y[425][0] = b;

 x[426][0] = 0.001;
 x[426][1] = -1.6641;
 x[426][2] = -1.3678;
 x[426][3] = 1.997;
 x[426][4] = 0.52283;
 y[426][0] = b;

 x[427][0] = 0.001;
 x[427][1] = -2.4349;
 x[427][2] = -9.2497;
 x[427][3] = 8.9922;
 x[427][4] = -0.50001;
 y[427][0] = b;

 x[428][0] = 0.001;
 x[428][1] = -3.793;
 x[428][2] = -12.7095;
 x[428][3] = 12.7957;
 x[428][4] = -2.825;
 y[428][0] = b;

 x[429][0] = 0.001;
 x[429][1] = -1.9551;
 x[429][2] = -6.9756;
 x[429][3] = 5.5383;
 x[429][4] = -0.12889;
 y[429][0] = b;

 x[430][0] = 0.001;
 x[430][1] = -0.69078;
 x[430][2] = -0.50077;
 x[430][3] = -0.35417;
 x[430][4] = 0.47498;
 y[430][0] = b;

 x[431][0] = 0.001;
 x[431][1] = 0.025013;
 x[431][2] = 3.3998;
 x[431][3] = -4.4327;
 x[431][4] = -4.2655;
 y[431][0] = b;

 x[432][0] = 0.001;
 x[432][1] = -4.3967;
 x[432][2] = 4.9601;
 x[432][3] = -0.64892;
 x[432][4] = -5.4719;
 y[432][0] = b;

 x[433][0] = 0.001;
 x[433][1] = -2.456;
 x[433][2] = -0.24418;
 x[433][3] = 1.4041;
 x[433][4] = -0.45863;
 y[433][0] = b;

 x[434][0] = 0.001;
 x[434][1] = -2.62;
 x[434][2] = -6.8555;
 x[434][3] = 6.2169;
 x[434][4] = -0.62285;
 y[434][0] = b;

 x[435][0] = 0.001;
 x[435][1] = -2.9662;
 x[435][2] = -10.3257;
 x[435][3] = 8.784;
 x[435][4] = -2.1138;
 y[435][0] = b;

 x[436][0] = 0.001;
 x[436][1] = -0.71494;
 x[436][2] = -4.4448;
 x[436][3] = 2.2241;
 x[436][4] = 0.49826;
 y[436][0] = b;

 x[437][0] = 0.001;
 x[437][1] = 0.6005;
 x[437][2] = 0.99945;
 x[437][3] = -2.2126;
 x[437][4] = 0.097399;
 y[437][0] = b;

 x[438][0] = 0.001;
 x[438][1] = 0.61652;
 x[438][2] = 3.8944;
 x[438][3] = -4.7275;
 x[438][4] = -4.3948;
 y[438][0] = b;

 x[439][0] = 0.001;
 x[439][1] = -5.4414;
 x[439][2] = 7.2363;
 x[439][3] = 0.10938;
 x[439][4] = -7.5642;
 y[439][0] = b;

 x[440][0] = 0.001;
 x[440][1] = -3.5798;
 x[440][2] = 0.45937;
 x[440][3] = 2.3457;
 x[440][4] = -0.45734;
 y[440][0] = b;

 x[441][0] = 0.001;
 x[441][1] = -2.7769;
 x[441][2] = -5.6967;
 x[441][3] = 5.9179;
 x[441][4] = 0.37671;
 y[441][0] = b;

 x[442][0] = 0.001;
 x[442][1] = -1.8356;
 x[442][2] = -6.7562;
 x[442][3] = 5.0585;
 x[442][4] = -0.55044;
 y[442][0] = b;

 x[443][0] = 0.001;
 x[443][1] = 0.30081;
 x[443][2] = 0.17381;
 x[443][3] = -1.7542;
 x[443][4] = 0.48921;
 y[443][0] = b;

 x[444][0] = 0.001;
 x[444][1] = 1.3403;
 x[444][2] = 4.1323;
 x[444][3] = -4.7018;
 x[444][4] = -2.5987;
 y[444][0] = b;

 x[445][0] = 0.001;
 x[445][1] = 0.26877;
 x[445][2] = 4.987;
 x[445][3] = -5.1508;
 x[445][4] = -6.3913;
 y[445][0] = b;

 x[446][0] = 0.001;
 x[446][1] = -6.5235;
 x[446][2] = 9.6014;
 x[446][3] = -0.25392;
 x[446][4] = -6.9642;
 y[446][0] = b;

 x[447][0] = 0.001;
 x[447][1] = -4.0679;
 x[447][2] = 2.4955;
 x[447][3] = 0.79571;
 x[447][4] = -1.1039;
 y[447][0] = b;

 x[448][0] = 0.001;
 x[448][1] = -2.564;
 x[448][2] = -1.7051;
 x[448][3] = 1.5026;
 x[448][4] = 0.32757;
 y[448][0] = b;

 x[449][0] = 0.001;
 x[449][1] = -1.3414;
 x[449][2] = -1.9162;
 x[449][3] = -0.15538;
 x[449][4] = -0.11984;
 y[449][0] = b;

 x[450][0] = 0.001;
 x[450][1] = 0.23874;
 x[450][2] = 2.0879;
 x[450][3] = -3.3522;
 x[450][4] = -0.66553;
 y[450][0] = b;

 x[451][0] = 0.001;
 x[451][1] = 0.6212;
 x[451][2] = 3.6771;
 x[451][3] = -4.0771;
 x[451][4] = -2.0711;
 y[451][0] = b;

 x[452][0] = 0.001;
 x[452][1] = -0.77848;
 x[452][2] = 3.4019;
 x[452][3] = -3.4859;
 x[452][4] = -3.5569;
 y[452][0] = b;

 x[453][0] = 0.001;
 x[453][1] = -4.1244;
 x[453][2] = 3.7909;
 x[453][3] = -0.6532;
 x[453][4] = -4.1802;
 y[453][0] = b;

 x[454][0] = 0.001;
 x[454][1] = -7.0421;
 x[454][2] = 9.2;
 x[454][3] = 0.25933;
 x[454][4] = -4.6832;
 y[454][0] = b;

 x[455][0] = 0.001;
 x[455][1] = -4.9462;
 x[455][2] = 3.5716;
 x[455][3] = 0.82742;
 x[455][4] = -1.4957;
 y[455][0] = b;

 x[456][0] = 0.001;
 x[456][1] = -3.5359;
 x[456][2] = 0.30417;
 x[456][3] = 0.6569;
 x[456][4] = -0.2957;
 y[456][0] = b;

 x[457][0] = 0.001;
 x[457][1] = -2.0662;
 x[457][2] = 0.16967;
 x[457][3] = -1.0054;
 x[457][4] = -0.82975;
 y[457][0] = b;

 x[458][0] = 0.001;
 x[458][1] = -0.88728;
 x[458][2] = 2.808;
 x[458][3] = -3.1432;
 x[458][4] = -1.2035;
 y[458][0] = b;

 x[459][0] = 0.001;
 x[459][1] = -1.0941;
 x[459][2] = 2.3072;
 x[459][3] = -2.5237;
 x[459][4] = -1.4453;
 y[459][0] = b;

 x[460][0] = 0.001;
 x[460][1] = -2.4458;
 x[460][2] = 1.6285;
 x[460][3] = -0.88541;
 x[460][4] = -1.4802;
 y[460][0] = b;

 x[461][0] = 0.001;
 x[461][1] = -3.551;
 x[461][2] = 1.8955;
 x[461][3] = 0.1865;
 x[461][4] = -2.4409;
 y[461][0] = b;

 x[462][0] = 0.001;
 x[462][1] = -2.2811;
 x[462][2] = -0.85669;
 x[462][3] = 2.7185;
 x[462][4] = 0.044382;
 y[462][0] = b;

 x[463][0] = 0.001;
 x[463][1] = -3.6053;
 x[463][2] = -5.974;
 x[463][3] = 10.0916;
 x[463][4] = -0.82846;
 y[463][0] = b;

 x[464][0] = 0.001;
 x[464][1] = -5.0676;
 x[464][2] = -5.1877;
 x[464][3] = 10.4266;
 x[464][4] = -0.86725;
 y[464][0] = b;

 x[465][0] = 0.001;
 x[465][1] = -3.9204;
 x[465][2] = 4.0723;
 x[465][3] = -0.23678;
 x[465][4] = -2.1151;
 y[465][0] = b;

 x[466][0] = 0.001;
 x[466][1] = -1.1306;
 x[466][2] = 1.8458;
 x[466][3] = -1.3575;
 x[466][4] = -1.3806;
 y[466][0] = b;

 x[467][0] = 0.001;
 x[467][1] = -2.4561;
 x[467][2] = -4.5566;
 x[467][3] = 6.4534;
 x[467][4] = -0.056479;
 y[467][0] = b;

 x[468][0] = 0.001;
 x[468][1] = -4.4775;
 x[468][2] = -13.0303;
 x[468][3] = 17.0834;
 x[468][4] = -3.0345;
 y[468][0] = b;

 x[469][0] = 0.001;
 x[469][1] = -4.1958;
 x[469][2] = -8.1819;
 x[469][3] = 12.1291;
 x[469][4] = -1.6017;
 y[469][0] = b;

 x[470][0] = 0.001;
 x[470][1] = -3.38;
 x[470][2] = -0.7077;
 x[470][3] = 2.5325;
 x[470][4] = 0.71808;
 y[470][0] = b;

 x[471][0] = 0.001;
 x[471][1] = -2.4365;
 x[471][2] = 3.6026;
 x[471][3] = -1.4166;
 x[471][4] = -2.8948;
 y[471][0] = b;



   }
   
   for(m=0;m<472;m++)
   {
      intermediaria_(x,w1,h1,m);
      intermediaria2_(w3,h1,h2,m);
      saida_(h2,w2,o);
      err=erro_saida_(o,y,m);
      fprintf(pFile,"\nPadrao>>%d",m);
      fprintf(pFile,"\ncalculado>>%f   \tdesejado>>%f  \tErro>>%2.9f",o[0],y[m][0],err);
      printf("\nPadrao>>%d",m);
      printf("\ncalculado>>%f   \tdesejado>>%f  \tErro>>%2.9f",o[0],y[m][0],err);
      emq=emq+err;

   }
   emq=emq/472;
   fitness=1000*(exp(-emq)*exp(-NINT1)*exp(-NINT2))+(1/(emq*NINT1*NINT2));
   fit[indic][0]=(indic+1);
   fit[indic][1]=fitness;
   fit[indic][2]=fit[indic][0];
   aux=indic+1;
   fit[indic][3]=aux;
   if(aux>1)
    fit[indic][3]=aux+fit[indic-1][3];
   else
     fit[indic][3]=aux;
   fit[indic][4]=emq;
   fit[indic][5]=(NINT-1);
   fit[indic][6]=(NINT1-1)+(NINT2-1)+(NENT-1)+NSAI;
   if(n == 4 )fit[indic][7]= NINT2-1;
   
   indic++;
   printf("\nemq>>%f",emq);
  
//   system("pause");
   fprintf(pFile,"\n\nErro Medio Quadratico>>%2.9f",emq);
   fprintf(pFile,"\n\n<<Pesos Camada Entrada Oculta>>");
   for(i=0;i<NENT;i++)
   {
     for(j=0;j<NINT1;j++)
     {
      fprintf(pFile,"\n");
      fprintf(pFile," w1[%d][%d]=%f",i,j,w1[i][j]);
     }
     fprintf(pFile,"\n");
   }
   fprintf(pFile,"\n\n<<Pesos Camada Oculta 1 - Oculta 2>>");
   for(i=0;i<NINT1;i++)
   {
     for(j=0;j<NINT2;j++)
     {
      fprintf(pFile,"\n");
      fprintf(pFile," w3[%d][%d]=%f",i,j,w3[i][j]);
     }
     fprintf(pFile,"\n");
   }
   fprintf(pFile,"\n\n<<Pesos Camada Oculta Saida>>");
   for(i=0;i<NINT2;i++)
   {
     for(j=0;j<NSAI;j++)
     {
      fprintf(pFile,"\n");
      fprintf(pFile," w2[%d][%d]=%f",i,j,w2[i][j]);
     }
     fprintf(pFile,"\n");
   }

}
  
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imprime_fitness(double fit[][10],FILE* pFile,int contador)
{
   int i;
   double soma_fit=0.0;
   double soma_nint=0.0;
   int contador_val=0;
   double var_nint=0.0;

   fprintf(pFile,"\n\nResumo da Geracao:%d",contador);
   fprintf(pFile,"\n_______________________________________________________________");
   n==3?
   fprintf(pFile,"\nIND  |Fitness      |Posto  |Acum    |EMQ      |NINT  |NT"):
   fprintf(pFile,"\nIND  |Fitness      |Posto  |Acum    |EMQ      |NINT1 |NINT2  |NT");
   fprintf(pFile,"\n_______________________________________________________________");

   for(i=0;i<individuos;i++)
   {
    fprintf(pFile,"\n%3.0f",fit[i][0]);
    fprintf(pFile,"    %3.4f",fit[i][1]);
    fprintf(pFile,"     %3.0f",fit[i][2]);
    fprintf(pFile,"     %3.0f",fit[i][3]);
    fprintf(pFile,"      %3.6f",fit[i][4]);
    fprintf(pFile,"   %3.0f",fit[i][5]);
    if(n==4)fprintf(pFile,"      %d",fit[i][7]);
    fprintf(pFile,"   %3.0f",fit[i][6]);
    soma_fit=soma_fit+fit[i][1];
    if(fit[i][5]>0)
    {
      soma_nint=soma_nint+fit[i][5];
      contador_val++;
    }
    sup=fit[i][3];

   }
   fprintf(pFile,"\n===============================================================\n");
   if (contador_val>0)
    {
     soma_nint_total=soma_nint_total+(soma_nint/contador_val);
     fprintf(pFile,"\nMedia do Fitness=%3.4f\n",soma_fit/contador_val);
     fprintf(pFile,"Media de Neuronios na Camada Intermediaria na Geracao=%3.2f\n",soma_nint/contador_val);
    }
   hist_fit[contador-1][0]=fit[0][1];
   if(contador_val>0)
    hist_fit[contador-1][1]=soma_fit/contador_val;
   else
     hist_fit[contador-1][1]=0.0;
   for(i=0;i<individuos;i++)
   {
     if(fit[i][5]>0)
        var_nint=var_nint+pow((fit[i][5]-(soma_nint/contador_val)),2)/contador_val;
   }
   fprintf(pFile,"Variancia de Neuronios na Camada Intermediaria na Geracao=%3.2f\n",var_nint);
   fprintf(pFile,"Desvio Padrao de Neuronios na Camada Intermediaria na Geracao=%3.2f\n",sqrt(var_nint));
   if(contador_val>0)
   {
    fprintf(pFile,"Erro Padrao de Neuronios na Camada Intermediaria na Geracao=%3.2f\n",sqrt(var_nint)/sqrt(contador_val));
    fprintf(pFile,"Coeficiente de Variacao de Neuronios na Camada Intermediaria na Geracao=%3.2f\n",sqrt(var_nint)/(soma_nint/contador_val));
   }
   if(contador==geracao)
     fprintf(pFile,"\n\nMedia de Neuronios na Camada Intermediaria na Simulacao=%3.2f\n",(soma_nint_total/contador));

 }


void ordena(double fit[][10],FILE* pFile)
{
 int min,i,j;
 double ch,ch1,ch2,ch3,ch4,ch5,AUX1,AUX2;
 for(i=0;i<=(individuos-2);i++)
 {
    min=i;
    for(j=i+1;j<=(individuos-1);j++)
    {
     AUX1=fit[j][1];
     AUX2=fit[min][1];
     if(AUX1>AUX2)
     {
       min=j;
     }
    }
    ch=fit[i][0];
    ch1=fit[i][1];
    ch2=fit[i][4];
    ch3=fit[i][5];
    ch4=fit[i][6];
    ch5=fit[min][0];
    fit[i][0]=ch5;
    fit[i][1]=fit[min][1];
    fit[i][4]=fit[min][4];
    fit[i][5]=fit[min][5];
    fit[i][6]=fit[min][6];
    fit[min][0]=ch;
    fit[min][1]=ch1;
    fit[min][4]=ch2;
    fit[min][5]=ch3;
    fit[min][6]=ch4;
  }
}

void zera_fitness(double fit[][10])
{
    int i,j;
    for(i=0;i<individuos;i++)
      for(j=0;j<gene;j++)
        fit[i][j]=0;
}

//selecao por roleta ponderada
void selecao(int** gen,char** gen_string,int gene_dec)
{
  int j,e,e1,e2,e3,e4,e5,hab_pc;
  e1=0;
  e=0;
  hab_pc=0;
  j=(individuos-1);
  //printf("\nj=");
  //printf("%d",j);

  while(e==e1)
  {
   time_t t;
   srand((unsigned) time(&t));
   e=(rand()%(individuos/2));
//////   printf("ind %d",individuos);
//////   printf("e %d", e);
   e2=(int)(fit[e][0]);
//////   printf("e2 %d",e2);
   e2--;
   if(e2>=(individuos/2))
    e2=1;
   e1=(rand()%(individuos/2));
//////   printf("e1 %d",e1);
   e3=(int)(fit[e1][0]);
//////   printf("e3 %d",e3);
   e3--;
   if(e3>=(individuos/2))
    e3=0;
//////    printf("\nasddas %d",e3);
  }
  e4=rand()%(individuos/2);
  e5=(int)(fit[e4][0]);
  e5--;
  if(e5<=(individuos/2)-1)
    e5=j-e5;
  //printf("\nIndividuo 1 Selecionado para o Cruzamento =");
  //printf("%d",e2);
  //printf("\nIndividuo 2 Selecionado para o Cruzamento =");
  //printf("%d",e3);
  //printf("\nIndividuo 3 Selecionado para o Cruzamento =");
  //printf("%d",e5);
  //cruzamento pai e mae
  //printf("\ncruza1");
//  if(e2 < 0) e2 = e2*(-1);
//  if(e3 < 0) e3 = e3*(-1);
//  if(e3 > 4 && n == 4) e5 = 4;
  
  //printf("\n %d %d %d", e2, e3, e5);
  //printf("t - %d",sizeof(gen[0][0]));
  cruzamento(e2,e3,e3,gen,gen_string,gene_dec);
  //printf("\ncruza");
  

}

void cruzamento(int i,int j,int pos,int** gen,char** gen_string,int gene_dec)
{
  int aux1,aux2,aux3,aux4,aux5,aux6;
  int aux7,aux8,aux9,aux10,aux11,aux12;
  //printf("\ncruza3");
  float pc;
  time_t t;
  srand((unsigned) time(&t));
  pc=(rand()%individuos);
  pc=(pc/100);
  //printf("\npc=%1.2f",pc);
  const int ale=(gene-12);
  const int gen_ale=(rand()%ale);
  const int gen_ale1=(rand()%ale);
 // printf("\ncruza4");
  //printf("\ngen1=%d",gen_ale);
 // printf("\ngen2=%d",gen_ale1);
 // printf("\npos=%d",pos);
 // printf("\ni=%d",i);
 // printf("\nj=%d",j);
  //getch();
//printf("\n %f", pc);
  if(pc<=0.6)
  {
  	
	//printf("\ncruza5 %d", sizeof(gen[0]));
    aux1=gen[i][gen_ale];
    
//////    printf("\ncruza8 %d %d",pos, gen_ale );
	gen[pos][gen_ale]=aux1;
	//printf("\ncruza9");
    aux2=gen[j][gen_ale1];
//////    printf("\ncruza10");
    gen[pos][gen_ale1]=aux2;
//////    printf("\ncruza7");
    //
    aux3=gen[i][gen_ale+1];
    gen[pos][gen_ale+1]=aux3;
    aux4=gen[j][gen_ale1+1];
    gen[pos][gen_ale1+1]=aux4;
    //
    aux5=gen[i][gen_ale+2];
    gen[pos][gen_ale+2]=aux5;
    aux6=gen[j][gen_ale1+2];
    gen[pos][gen_ale1+2]=aux6;
    //
    aux7=gen[i][gen_ale+3];
    gen[pos][gen_ale+3]=aux7;
    aux8=gen[j][gen_ale1+3];
    gen[pos][gen_ale1+3]=aux8;
    //
    aux9=gen[i][gen_ale+4];
    gen[pos][gen_ale+4]=aux8;
    aux10=gen[j][gen_ale1+4];
    gen[pos][gen_ale1+4]=aux10;
    //
    aux11=gen[i][gen_ale+5];
    gen[pos][gen_ale+5]=aux11;
    aux12=gen[j][gen_ale1+5];
    gen[pos][gen_ale1+5]=aux12;

//printf("\ncruza6");
   // for(i=1;i<=(individuos);i++)
   // {
   //  if(pc<=tx_trans)
   //  {
      // gen_string[individuos-i][0]='.';
      // gen_string[individuos-i][1]='f';
      // gen_string[individuos-i][2]='[';
      // gen_string[individuos-i][3]='F';
      // gen_string[individuos-i][4]='f';
      // gen_string[individuos-i][5]='n';
      // gen_string[individuos-i][6]='B';
      // gen_string[individuos-i][7]=']';
      // gen_string[individuos-i][gene_dec+1]='V';
    // }
    //}

    //tx_trans=tx_trans+0.01;

  }
}

void mutacao(int** gen)
{
   int e6,e7;
   float pm;
   time_t t;
   srand((unsigned) time(&t));
   pm=(rand()%10);
   pm=(pm/100);

   e6=rand()%(individuos/2);
   e7=(int)(fit[e6][0]);
   e7--;
   if(e7<=(individuos/2)-1)
    e7=((individuos-1)-e6);
 //  printf("e7=%d",e7);
 //  printf("pm=%1.2f",pm);
   const int ale=(gene-1);
   const int gen_ale=(rand()%ale);
  // printf("\ngen=%d",gen_ale);
 //  getch();
   if((gen[e7][gen_ale]==0)&&(pm<=0.1))
     gen[e7][gen_ale]=1;
   else
     gen[e7][gen_ale]=0;
}

void transloca(char** gen_string,int gene_dec)
{
    int i;
   // printf("vez=%d",vez);
   // getch();
    int cont_trans=2;

    for(i=1;i<=(individuos/cont_trans);i++)
    {
      if((tx_trans>=0.1))
      {
       gen_string[individuos-i][0]='.';
       gen_string[individuos-i][1]='f';
       gen_string[individuos-i][2]='[';
       gen_string[individuos-i][3]='F';
       gen_string[individuos-i][4]='f';
       gen_string[individuos-i][5]='n';
       gen_string[individuos-i][6]='B';
       gen_string[individuos-i][7]=']';
       gen_string[individuos-i][gene_dec+1]='V';
      }
    }

    tx_trans=tx_trans+0.01;
}

void imprime_cabec(FILE* pFile)
{
  fprintf(pFile,"\nPercentagem de Regras Validas=%3.2f\n %",(contregrasval/(individuos*geracao))*100.0);
  fprintf(pFile,"\n\nHistorico do Fitness na Simulacao");
  fprintf(pFile,"\n________________________________________________________________");
  fprintf(pFile,"\n Geracao  |Melhor Fitness         |Fitness Medio");
  fprintf(pFile,"\n________________________________________________________________");
}

void imprime_hist_fitness(double hist_fit[][2],FILE* pFile,int contador)
{

   fprintf(pFile,"\n%d",contador);
   fprintf(pFile,"\t\t%3.4f",hist_fit[contador-1][0]);
   fprintf(pFile,"\t\t%3.4f",hist_fit[contador-1][1]);
}

// MAPEAMENTO
int Mapeamento(int NENT, int NSAI, int NR[],int SIZE, char TIPO[], FILE *pFile){
    int q = 0;
	//FILE *txt = pFile; // apaga possiveis lixos no arquivo
    //fclose(txt);

	///////////////////////////////////////////////
	////////////////// STRUCT /////////////////////
	///////////////////////////////////////////////
	struct ED{ // estrutura Lista

	   char ele;  // Dados
	   struct ED *Prox; // Proximo
	   struct ED *Ant; // Anterior

	};

	struct ED* Primeiro[1000];
	struct ED *Ultimo[1000];  // Ultimo elemento
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	int size = 0;

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////
	///////////////// METODOS ED /////////////////
	//////////////////////////////////////////////

	int pList = 0;


	void QueueNull(){ // inicaliza fila NULL

	     struct ED *aux;

	     aux = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED)); // aoca??o de memoria
	     Primeiro[q] = aux;
	     Ultimo[q] = Primeiro[q];
	     Primeiro[q]->Ant = NULL;

	}

	void InsertInFronte(char elemento){ // insere dentro da ED no inicio

	     struct ED *aux;
	     aux = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED)); // aloca??o de memoria

				 aux->ele = elemento;
	             aux -> Ant = NULL;
	             aux -> Prox = NULL;

	             if(!Primeiro[q]->ele){
	             	//Primeiro[q] = Primeiro[q]->Prox;
	                Primeiro[q] = Ultimo[q] = aux;
	                   }

	             else{
	             aux-> Prox = Primeiro[q];
	             Primeiro[q] -> Ant = aux;
	             Primeiro[q] = aux;
	                             }
	                             size++;
	}

	void ReplaceOnPosition(char elemento,int pos){

	     int p = 1;

	     struct ED *aux = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED));
	     aux = Primeiro[q];
	     int flag = 1;
	     while(aux!=NULL){
	     if(p == pos){
	                 aux->ele = elemento;
	                 flag = 0;
	                 aux = NULL;

	                            }
	   else{
	        aux = aux->Prox;
	        p++;
	        }

	                   }
	   if( flag == 1 ){
//////	     printf("O Elemento %c N?o existe na ED\n", elemento);
	   }
	}

	void InsertMiddle(char elemento, int posicao){ // insere em qualquer lugar
		char count = posicao;
		count--;
		struct ED *aux;
		struct ED *aux2;
		aux = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED)); //aloca??o de memoria para aux
		aux2 = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED)); //aloca??o de memoria para aux2

		if( count < 0){
			printf("\nn?o Existe Esta Posi??o");
		}else{
			if(count == 0){
			aux->ele = elemento;
			aux->Prox = Primeiro[q];
			aux->Ant  = NULL;
			Primeiro[q]->Ant = aux;

			aux = Primeiro[q];
			size++;
		}else{

		aux->Prox = Primeiro[q]->Prox;
		aux->Ant  = Primeiro[q]->Ant;
		aux = Primeiro[q];

			}
			while(count!=0){

				aux = aux->Prox;
				count--;
			}
			aux2->ele = elemento;

			aux->Ant->Prox =aux2;
			aux2->Prox = aux;
			aux2->Ant = aux->Ant;
			aux->Ant = aux2;
			size++;
		}
	}

	void InsertLast(char elemento){ // insere pelo final
	     struct ED *aux;

	     aux = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED));
	     aux->ele = elemento;

	      if(!Primeiro[q]->ele){
	             	Primeiro[q] = Primeiro[q]->Prox;
	                Primeiro[q] = Ultimo[q] = aux;
	                   size = 0;
					   size++;
					   }else{

	     Ultimo[q]->Prox = aux;
	     aux->Ant = Ultimo[q];
	     aux->Prox = NULL;
		 Ultimo[q] = aux;
		 size++;
		 }
	}

	void ShowStack(){ // Imprime (PILHA)
		//FILE *save;
		//save = pFile;
	  	struct ED *aux;
	  	aux = Primeiro[q];//->Prox;
	  	if(!aux->ele){
		  printf("A ED est? Vazia !");

		  }
		  else{
		  //	printf("  ||  ");
		    fprintf(pFile,"  ||  ");
		  while(aux !=NULL){

			//printf("\nElemento: %c \n", aux->ele);
		//	printf(" %c ",aux->ele);
	    	fprintf(pFile," %c ",aux->ele);
			aux = aux->Prox;

	  		   }
	  	  }
	  		  if(q == NENT-1){
	  	//	  	printf("  ||");
	  		    fprintf(pFile,"  ||");
	  	//	   	printf("\n\n");
	  		    fprintf(pFile,"\n\n");
			  }
	}
	/*
	void SaveStack(){ // Imprime (PILHA)
		FILE *save;
		save = pFile;

	  	struct ED *aux;
	  	q = 0;
	  	aux = Primeiro[q];
	  	if(!aux->ele){

		  printf("\n\nNenhum dado\nOs dados n?o podem ser salvos\n\n");

		  }
		  else

		  	while(aux !=NULL){

			fprintf(save," %c ", aux->ele);

	    	aux = aux->Prox;

	  		}

	  fprintf(save,"\n");
	 // fclose(save);
	}
	*/
	/*
	void ShowQueue(){ // imprime (Fila)

	  	struct ED *aux;
	  	aux = Ultimo[q];

		  if(!aux->ele){
		  printf("A ED est? Vazia !");
		                 }
		  else{

		  while(aux !=NULL){
	    	printf(" %c ", aux->ele);
	    	aux = aux->Ant;
	      }
	    }
	}
	*/
	int Search(char elemento){ // pesquisa por elemntos dentro da ED

	   int pos = 1;
	   struct ED *aux = (struct ED*)calloc(sizeof(struct ED),sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1;

	   while((Primeiro[q]->ele == 'f') && (Primeiro[q]->Prox) && Primeiro[q]->ele == 'f'){

						aux = Primeiro[q]->Prox;
	                    pos++;
	    }

	   while(aux!=NULL){
	     if(aux->ele == elemento || (aux->ele == 'f' && aux->Ant->Ant->ele == ']')){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                 flag = 0;
	                 aux = NULL;

	                            }
	     else{
	        aux = aux->Prox;
	        pos++;
	     }

	    }

	   if( flag == 1 ){
	    // printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     pos = 0;
	   }

	   return pos;
	}

	int SearchR36(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1;

	   while(pos_start >= c && aux!=NULL){
	   		aux = aux->Prox;
	   		pos++;
	   		c++;
	   }

	   while(aux!=NULL){
	     if(aux->ele == elemento && aux->Prox->ele == ']' && aux->Prox!=NULL){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                 flag = 0;
	                 aux = NULL;
	      }else{
	       		aux = aux->Prox;
	        	pos++;
	        	if(aux->ele == ']' && aux->Prox->ele == 'f'){
	        		flag = 1;
	        		aux = NULL;
	        		pos = 0;
				}
	       }
	    }

	   if( flag == 1 || pos >= size){
	     //printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     pos = 0;
	   }
	   return pos;
	}

	/*
	int SearchR362(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1;

	   while(aux!=NULL){
	     if(aux->ele == elemento && aux->Ant->ele == ']'){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                 flag = 0;
	                 aux = NULL;
	      }else{
	       		aux = aux->Prox;
	        	pos++;
	        	if(aux->ele == ']' && aux->Prox->ele == 'f'){
	        		flag = 1;
	        		aux = NULL;
	        		pos = 0;
				}
	       }
	    }

	   if( flag == 1 || pos >= size){
	     //printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     pos = 0;
	   }
	   return pos;
	}

	*/

	int SearchR322(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1;

	   while(pos_start >= c && aux!=NULL){
	   		aux = aux->Prox;
	   		pos++;
	   		c++;
	   }

	   while(aux!=NULL){
	     if(aux->ele == elemento && aux->Ant->Ant->ele == ']' || (aux->ele == elemento && aux->Ant->Ant->ele == '[')){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                 flag = 0;
	                 aux = NULL;
	      }else{
	       		aux = aux->Prox;
	        	pos++;
	       }
	    }

	   if( flag == 1 || pos >= size){
	  ////////   printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     pos = 0;
	   }
	   return pos;
	}

	int SearchR341(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1, flag2 = 0;


	   while(pos_start >= c && aux!=NULL){
	   		aux = aux->Prox;
	   		pos++;
	   		c++;
	   }

	   while(aux!=NULL){
	   	 if(!aux->Prox){
	            flag = 1;
	            break;
	     }
	     if((aux->ele == elemento && aux->Prox->ele == 'B' && aux->ele != 'B')){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                // printf("|%c|",aux->ele);
					 flag = 0;
	                 if(!aux->Prox){
	                 	flag2 = 1;
	                 	break;
					 }
	                 aux = NULL;

	      }else{
	      		flag = 1;
	       		aux = aux->Prox;
	        	pos++;
	       }

	    }

	   if( flag == 1 || flag2 == 1){
	     //printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     if(Ultimo[q]->ele == 'f'){
	                       Ultimo[q]->ele = 'n';
	     }
	    // printf("pas");
	     pos = 0;
	   }
	   return pos;
	}

	int SearchfF(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1, flag2 = 0;


	   while(pos_start >= c && aux!=NULL){
	   		aux = aux->Prox;
	   		pos++;
	   		c++;
	   }

	   while(aux!=NULL){
	   	 if(!aux->Prox){
	            flag = 1;
	            break;
	     }
	     if(aux->ele == elemento && (aux->Ant->Ant->ele == ']' || aux->Ant->Ant->ele == '[') && aux->Prox->ele != 'B'){
					 flag = 0;
	                 if(!aux->Prox){
	                 	break;
					 }
	                 aux = NULL;

	      }else{
	      		flag = 1;
	       		aux = aux->Prox;
	        	pos++;
	       }

	    }
		if (flag == 1){
			pos = 0;
		}
	   return pos;
	}

	/*
	int SearchR342(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1, flag2 = 0;


	   while(pos_start >= c && aux!=NULL){
	   		aux = aux->Prox;
	   		printf("%dp",pos);
	   		pos++;
	   		c++;
	   }

	   while(aux!=NULL){
	   	 if(!aux->Prox){
	            flag = 1;
	            break;
	     }
	     if((aux->ele == elemento && aux->Prox->ele == 'F')){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                 //printf("|%c|",Ultimo[q]->ele);
//////	                 printf("%dvvvvv",pos);
					 flag = 0;
	                 if(!aux->Prox){
	                 	flag2 = 1;
	                 	pos = 0;
	                 	break;
					 }
	                 aux = NULL;

	      }else{
	      		flag = 1;
	       		aux = aux->Prox;
	        	pos++;
	       }

	    }

	   if( flag == 1 || flag2 == 1){
	     //printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     if(Ultimo[q]->ele != 'B'){
	                    Ultimo[q]->ele = 'n';
	     }
	     //printf("pas");
	     //pos = 0;
	   }
	   return pos;
	}
	*/
	int SearchR33(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
	   int pos = 1;
	   int c = 1;
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1;

	   while(pos_start >= c && aux!=NULL){
	   		aux = aux->Prox;
	   		pos++;
	   		c++;
	   }

	   while(aux!=NULL){
	     if(aux->ele == elemento && aux->Ant->Ant->ele == ']'){
	                 //printf("O Elemento %c est? contido na ED \n", elemento);
	                 flag = 0;
	                 aux = NULL;
	      }else{
	       		aux = aux->Prox;
	        	pos++;
	      }
	    }

	   if( flag == 1 ){
	     //printf("O Elemento %c N?o existe na ED\n\n", elemento);
	     pos = 0;
	   }
	   return pos;
	}
	
	
int SearchAllf(char elemento, int pos_start){ // pesquisa por elemntos dentro da ED
   int pos = 1;
   int c = 1;
   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
   aux = Primeiro[q];
   int flag = 1;

   while(pos_start >= c && aux!=NULL){
   		aux = aux->Prox;
   		pos++;
   		c++;
   }

   while(aux!=NULL){
     if(aux->ele == elemento && !(aux->Ant->ele == ']')){
                 //printf("O Elemento %c est? contido na ED \n", elemento);
                 flag = 0;
                 aux = NULL;
      }else{
       		aux = aux->Prox;
        	pos++;
       }
    }

   if( flag == 1 || pos >= size){
  ////////   printf("O Elemento %c N?o existe na ED\n\n", elemento);
     pos = 0;
   }
   return pos;
}
	/*
	void Remove(char elemento){ // remove elemento especifico
	   struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	   aux = Primeiro[q];
	   int flag = 1;

	   while(aux!=NULL){
	     if(aux->ele == elemento){
//////	        printf("\n\nElemento %c Removido...\n\n",elemento);
	        if(aux->Prox==NULL){
	            Ultimo[q] = aux->Ant;
	            aux->Ant->Prox = NULL;
	            aux = NULL;
	            flag = 0;

			}else if(aux->Ant == NULL){
	              Primeiro[q] = aux->Prox;
	              Primeiro[q]->Ant = NULL;
	              aux = NULL;
	              }
	        else{
	            aux->Ant->Prox = aux->Prox;//O proximo do item de tr?s
	            aux->Prox->Ant = aux->Ant;//Item posterior de aux pega o item anterior ao aux
	            flag = 1;
	            aux = NULL;
	            size--;
	            }
		  }
	      else{
	        aux = aux->Prox;
	           }
	   }// fim while
	   free(aux); // liberar memoria

	   if(flag == 1){
//////	     printf("\n\nO Elemento %c N?o esta contido na ED\n\n",elemento);
	     size--; // tamanho-1
	    }
	}

	*/
	/*
	void RemoveFronte(){

	     struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	  	 aux = Primeiro[q];//->Prox;
	  	 if(!aux->ele){ // verifica se a exitse elementos
		 printf("A ED est? Vazia !\n");
		  }else if(!aux->Prox)
		  {
		  QueueNull();
		  }else {
	           Primeiro[q]->Prox->Ant = NULL;
	           Primeiro[q] = Primeiro[q]->Prox;
	           }
	}
	*//*
	void RemoveLast(){

	      struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	  	 aux = Primeiro[q];//->Prox;
	  	 if(!aux->ele){ // verifica se a exitse elementos
		 printf("A ED est? Vazia !\n");
		  }else if(!aux->Prox)
		  {
		  QueueNull();
		  }
		  else {
	           Ultimo[q]->Ant->Prox = NULL;
	           Ultimo[q] = Ultimo[q]->Ant;
	           }
	}
	*//*
	void RemoveAll(){ // deleta tudo
	          QueueNull();
//////	          printf("ED Vazia ! \n");
	          system("pause<<NULL");
	}
	*/
	int nf(int NENT, int NSAI, int NR, int tipo){
		int c = 0;
		struct ED *aux = (struct ED*) calloc (NENT,sizeof(struct ED));
	  	aux = Primeiro[0];
	  	if(!aux->ele){
		  printf("A ED est? Vazia !");
		}
		else
		while(aux != NULL){

			if(aux->ele == 'f' && tipo == 1) c++;
			else if(aux->ele == 'n') c++;
				aux = aux->Prox;
	  	}
	  	return (c-1-NSAI); // "1" por causa que os neuronios s? s?o gerados na 1 string e s?o confirmados nela msm
	}
	//////////////////////////////////////////////
	//////////////////////////////////////////////

void SetQ(int x){
	q = x;
}

void R1(int NENT){ // S -> *
int x = Search('S');
     if(x != 0){
          ReplaceOnPosition('.',1);
     }
}

void R2(int NENT){  // * -> (f..f)*NENT
        int x = Search('.');

		if(x != 0){
           ReplaceOnPosition('f', 1);	
		}
}

void R31(){   // f -> [f
      int x = Search('f');
      if(x != 0){
           InsertInFronte('[');
      }
}

void R32(int NR){   // f -> f F f
     int c = NR;
     int x = 1;
          if(NSAI == 1){
				  x = Search('f');
                  InsertMiddle('F',x+1);
				  InsertMiddle('f',x+2);
           }else
		    if(NSAI > 1){
           		while( x != 0 && c > 0){
				  x = SearchR322('f',x);
				  if(x == 0)break;
                  InsertMiddle('F',x+1);
				  InsertMiddle('f',x+2);
				  c--;
				}
		   }
} // fim R32

void R322(int SIZE,int MAX,int NSAI){   // f -> f F f
    int Bigger;
	if(NSAI >= MAX){Bigger = NSAI;}else{Bigger = MAX;}
	 int c = Bigger;
     int x = 1;
            
		    	int i;
		    	for(i = 0; i < SIZE ; i++){
		    		x = 1; c = Bigger;
	           		while( x != 0 && c > 0){
					  x = SearchR322('f',x);
				//	  printf("%d while\n",x);
				//	  system("pause");
					  if(x == 0)break;
	                  InsertMiddle('F',x+1);
					  InsertMiddle('f',x+2);
					  c--;
					  ShowStack();
					}
				//	printf("i>>%d for\n",i);
				//	system("pause");
				}
		   
} // fim R32

void R33(int NSAI){  // f -> f F
     int c = 1;
     int x = 1;
          while(c <= NSAI){
          		while(x != 0){
                  x = SearchR33('f',x);
                  if( x == 0 )break;
				  InsertMiddle('F',x+1);
				  //cout << endl << endl << x  <<" < - "<< endl;
                  c++;
				  }
        		}
}

//////////////////////////////////////////////////
void R34(int NENT, int NSAI, int NR){
	int nrr = NR;
		int x = 1;
		while(x != 0 ){

                  x = SearchR341('f',x); // coloca n B
                   if( x == 0 ){
				   break;
				   }else{
				  	ReplaceOnPosition('n',x);
				   }
        }
        x = 1;
        while(x != 0 && nrr > 0){

                  x = SearchfF('f',x);
                   if( x == 0 ){
                   	
				   break;
				   }else{
				  	ReplaceOnPosition('n',x);
				  	nrr--;
				   	
                    }
                   }
    	if(NENT!=1){ // verifica se existe mais de 1 string v?lida
				   		while(q < NENT){
                           Ultimo[q]->Ant->ele = 'n';
                           q++;
                		}
    }
}

void R344(int NENT, int NSAI, int N, int SIZE,int MAX, int NR[]){
	
	int Bigger;
	int Vet[SIZE+1];
	Vet[SIZE] = NSAI;
	
	int k;
	
	for(k = 0 ; k < SIZE ; k++){
		Vet[k] = NR[k];
//		printf("k%d\n",Vet[k]);
	}
//	printf("k%d\n",Vet[SIZE]);
	
	if(NSAI >= MAX){Bigger = NSAI;}else{ Bigger = MAX;}
	
	int nrr = N;
		int x=1;
		int i=0;
		int j=0;
		
//		printf("%d Size", SIZE);
	//for(i = 0 ; i < SIZE ; i++){
	//	printf("%d\n\n",NR[i]);
	//}
		
		while(x != 0 ){
                   x = SearchR341('f',x); // coloca n B
                   if( x == 0 ){
				   break;
				   }else{
				   	if(Vet[j]-i > 0){
				  		ReplaceOnPosition('n',x);
				   	}else ReplaceOnPosition('f',x); //
				   }
				   
////////				   printf("i: %d  j: %d x: %d\n",i,j,x);
				   j++;
				   	if(j >= SIZE+1){
						i++;
						j = 0;
					}
					if(i >= Bigger){
						i=0;
						return;
					}
//        system("pause");
//        ShowStack();    XX
		}
}
/////////////////////////////////////////////////////////////

void R36(int NENT,int NSAI){
	int x=1;
	while( x != 0 ){
		x = SearchR36('f',x);
			if(x != 0){
				InsertMiddle('B',x+1);
			}
		}
		x = 1;
		if(NENT > 1){
			
		InsertLast('B');			
		}
}

void R366(int NENT,int NSAI){
	int x=1;

	while( x != 0 ){
		int x = SearchAllf('f', x);
			if(x != 0){
				InsertMiddle('B',x+1);
			}else{
				return;
			}
		}
		
		x = 1;
		if(NENT > 1){
		InsertLast('B');			
		}
}

void R4(int NSAI, int MAX){  // [ -> [Ff]
          int c = 1;
          int Bigger;
if(NSAI >= MAX){Bigger = NSAI;}else{ Bigger = MAX;}
          while(c <= Bigger){
                  int x = Search('[');
                  InsertMiddle('F',x+1);
                  InsertMiddle('f',x+2);
                  InsertMiddle(']',x+3);
                  c++;
          }
}// fim R4

////////////////////////////////////////////
////////////////////////////////////////////


	/////////////////////////////////////////////////////////////
	//////////////////////// FIM REGRAS /////////////////////////
	/////////////////////////////////////////////////////////////


//////    printf("\n\n\n\nMAPEAMENTO GENOTIPO / FENOTIPO");
    fprintf(pFile,"\n\n\n\nMAPEAMENTO GENOTIPO / FENOTIPO\n\n\n");
     
    int i;
	int MAX = NR[0];
	for(i = 0 ; i < SIZE ; i++){
		if(MAX < NR[i]) MAX = NR[i];
	}

     QueueNull(); // Inicializa a Lista

     InsertLast('S'); // COLOCA 'S' NA PILHA

	 // 1to1
     //////////////////////////////////////////////////////////////
	 if(NENT == 1 && NSAI == 1 && tipo == 1){ // Type 1 to 1
	// printf("\n\nTIPO 1x1\n\n");
	 fprintf(pFile,"\n\nTIPO 1x1\n\n");
	    q = 0;
             if(NR[0] == 1){
             	// CASO ESPECIAL NENT = 1 | NSAI = 1 | NR = 1

				   QueueNull();
                   InsertLast('S'); // COLOCA 'S' NA PILHA

				   R1 (NENT);
                   ShowStack();

                   R2(NENT);
                   ShowStack();

                   R31();
                   ShowStack();

                   R4 (1,MAX);
                   ShowStack();

                   R32(MAX);
                   ShowStack();

             }else{
             	// CASO NENT = 1 | NSAI = 1 | NR = n

				   QueueNull();
                   InsertLast('S'); // COLOCA 'S' NA PILHA

				   R1 (NENT);
                   ShowStack();

                   R2(NENT);
                   ShowStack();

				   R31();
                   ShowStack();

                   R4(1,MAX);
                   ShowStack();

                   R32(NSAI);
                   ShowStack();

                   R33(NSAI);
                   ShowStack();

             }// fim else

     } // Fim Type 1 to 1
	 //1toN
     //////////////////////////////////////////////////////////////////
	 if(NENT == 1 && NSAI >= 1 && tipo == 2){
	 //	printf("\n\nTIPO 1xN\n\n");
	 	fprintf(pFile,"\n\nTIPO 1xN\n\n");

		q = 0;
		QueueNull();
        InsertLast('S'); // COLOCA 'S' NA PILHA

		R1(NENT);
		ShowStack();

	    R2(NENT);
     	ShowStack();

        R31();
        ShowStack();

        R4(NSAI,MAX);
        ShowStack();

        R32(NSAI);
        ShowStack();

        R36(NENT, NSAI);
        ShowStack();

        R34(NENT,NSAI,MAX);
        ShowStack();


	 } // fim tipo 1.2

     if(NENT >= 1 && NSAI == 1 && tipo == 3){
    //	printf("\n\nTIPO Nx1\n\n");
	    fprintf(pFile,"\n\nTIPO Nx1\n\n");

			// para o Primeiro[0]

			q = 0;
            //  QueueNull();
		    //  InsertLast('S'); // COLOCA 'S' NA PILHA
            while(q+1 <= NENT){
			 	//q++;
              	QueueNull();
              	InsertLast('S'); // COLOCA 'S' NA PILHA DAS OUTRAS STRINGS CASO EXISTAM
              	q++;
			}

        	q = 0;
            while(q+1 <= NENT){
			 	R1(NENT); // S -> .  | R1 ? APLICADO NENT VEZES
		     	ShowStack();
			 	q++;
			}

			 q = 0;
			while(q+1 <= NENT){
				R2(NENT); // . -> f  | R2 ? APLICADO NENT VEZES
		    	ShowStack();
				q++;
			}

			q = 0;
			while(q+1 <= NENT){
				R31(); // f -> [f  | R31 ? APLICADO EM CADA STRING
		    //	printf("  ||  ");
		    	ShowStack();
				q++;
			}

			q = 0;
			while(q+1 <= NENT){
				R4(NSAI,MAX); // [ -> [ F f  | NR > NSAI ENT?O ? APLICADO NR VEZES, SE N?O, ? APLICADO NSAI VEZES EM CADA STRING
		    	ShowStack();
				q++;
			}

			q = 0;
			R32(NSAI);  // f -> f F f  | APLICADO NSAI VEZES SOMENTE PARA A PRIMEIRA STRING
		    ShowStack();
			q++;
			while(q+1 <= NENT){
				ShowStack();
		    	q++;
			}

			q = 0;
			while(q+1 <= NENT){
				R36(NENT, NSAI); // f -> f B  | APLICADO NR VEZES PARA O ULTIMO F DE CADA BLOCO E APLICADO PARA O PRIMEIRO f DE
		    	ShowStack();
				q++;
			}

        	q = 0;
			while(q+1 <= NENT){
				if(q == 0){
                     R34(NENT,NSAI,MAX); // f -> n  | APLICADO PARA OS PRIMEIROS f DE CADA BLOCO E APLICADO NSAI VEZES PARA O PRIMEIRO f
                     q = 0;
                }
		     	ShowStack();
			 	q++;
			}
	 } //  fim tipo 1.3

     if(NENT >= 1 && NSAI >= 1 && tipo == 4){
     //	printf("\n\nTIPO MxN\n\n");
	    fprintf(pFile,"\n\nTIPO MxN\n\n");
     	             // para o Primeiro[0]
			q = 0;
            QueueNull();
			InsertLast('S'); // COLOCA 'S' NA PILHA

            while(q+1 <= NENT){
			  q++;
              QueueNull();
              InsertLast('S'); // COLOCA 'S' NA PILHA DAS OUTRAS STRINGS CASO EXISTAM
			}

        	q = 0;
            while(q+1 <= NENT){
				R1(NENT); // S -> .  | R1 ? APLICADO NENT VEZES
			    ShowStack();
				q++;
			}

			q = 0;
			while(q+1 <= NENT){
				R2(NENT); // . -> f  | R2 ? APLICADO NENT VEZES
		    	ShowStack();
				q++;
			}

			q = 0;
			while(q+1 <= NENT){
				R31(); // f -> [f  | R31 ? APLICADO EM CADA STRING
		    	ShowStack();
				q++;
			}

			q = 0;
			while(q+1 <= NENT){
				R4(NSAI,MAX); // [ -> [ F f  | NR > NSAI ENT?O ? APLICADO NR VEZES, SE N?O, ? APLICADO NSAI VEZES EM CADA STRING
		     	ShowStack();
			 	q++;
			}
			
			q = 0;	
			while(q+1 <= NENT){ 
				R322(SIZE,MAX,NSAI);  // f -> f F f  | APLICADO NSAI VEZES SOMENTE PARA A PRIMEIRA STRING /////// alterado
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
			while(q+1 <= NENT){ 
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
			while(q+1 <= NENT){ 
				if(q == 0){
					R344(NENT,NSAI,5,SIZE,MAX,NR);
                     //R34(NENT,NSAI,NR); // f -> n  | APLICADO PARA OS PRIMEIROS f DE CADA BLOCO E APLICADO NSAI VEZES PARA O PRIMEIRO f  
                     q = 0;
                }
		     	ShowStack();
			 	q++;
			}
//        	q = 0;
//			while(q+1 <= NENT){
//				if(q == 0){
//                     R34(NENT,NSAI,NR); // f -> n  | APLICADO PARA OS PRIMEIROS f DE CADA BLOCO E APLICADO NSAI VEZES PARA O PRIMEIRO f
//                     q = 0;
//                }
//		     	ShowStack();
//			 	q++;
//			}

	 } // fim tipo 1.4
     ///////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////
//	printf("| Numero de Neuronios Camada Intermediaria 1 -> %d \n", NR[0]);
	if(n == 4)
//		printf("| Numero de Neuronios Camada Intermediaria 2 -> %d \n", NR[1]);
//	printf("| Numero de Neuronios por Ramificacao -> %d \n", NR);
//////    printf("| NENT -> %d \n", NENT);
	//printf("| NSAI -> %d \n\n", NSAI);

	fprintf(pFile,"\n");
	fprintf(pFile,"| Numero de Neuronios Camada Intermediaria 1 -> %d \n", NR[0]);
	if(n == 4)
		fprintf(pFile,"| Numero de Neuronios Camada Intermediaria 2 -> %d \n", NR[1]);
//    fprintf(pFile,"| Numero de Neuronios por Ramificacao -> %d \n", NR);
	fprintf(pFile,"| NENT -> %d \n", NENT);
    fprintf(pFile,"| NSAI -> %d \n\n", NSAI);

	return 0;

} // Fim Mapeamento

int NR_RAND(){ // Neuronios Por Ramifica??o 1-10
	srand( (unsigned)time(NULL) );
	int x = 2 + ( rand() % 8);
	return x;
}
int NR_RAND_N(int n){ // Neuronios Por Ramifica??o 1-10
    srand( (unsigned)time(NULL) );
    int x = n + ( rand() % (10-n) );
    return x;
}

