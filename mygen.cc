#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <string.h>

#define PI 3.141592654

int Number_Items, Num_Large_Itemset, Average_LargeItemset_Size;
int Average_Transaction_Size, Num_Transactions;
double Correllation_Level;

int distinct_items, max_itemlen;

char fname[300];
int *fd=NULL;
FILE *fp=NULL;
int nproc=1;
int binary_output = 0;
#define ITEMSIZE sizeof(int)
#define OUTBUFSIZE 65536
int curr_proc=0;
int randseed=100;

int *BLKS;
int OUTBUF[OUTBUFSIZE];
int outbuf_pos=0;

void parse_args(int argc, char **argv)
{
   extern char * optarg;
   int c, i, blk;
   char *aux;
   
   if (argc < 6){
      fprintf(stderr,"USAGE: data_set -d Number of Tranasactions");
      fprintf(stderr,"-t Average Size of Transactions");
      fprintf(stderr,"-i Average size of the maximal potentially large itemsets");
      fprintf(stderr,"-l Number of maximal potentially large itemsets");
      fprintf(stderr,"-n Number of items");
      fprintf(stderr,"-c correllation_level");
      fprintf(stderr,"-b binary output flag");
   }
   else{
      while ((c=getopt(argc,argv,"d:t:i:l:n:c:o:p:br:"))!=-1){
         switch(c){
         case 'd':
            Num_Transactions = atoi(optarg);
            break;
         case 't':
            Average_Transaction_Size = atoi(optarg);
            break;
         case 'i':
            Average_LargeItemset_Size = atoi(optarg);
            break;
         case 'l':
            Num_Large_Itemset = atoi(optarg);     
            break;
         case 'n':
            Number_Items = atoi(optarg);
            break;
         case 'c':
            Correllation_Level = atof(optarg);
            break;
         case 'b':
            binary_output = 1;
            break;
         case 'p':
            nproc = atoi(optarg);
            break;
         case 'o':
            strcpy(fname, optarg);
            break;
         case 'r':
            randseed = atoi(optarg);
            break;
         }
      }
   }
   fd = (int *)malloc(nproc*ITEMSIZE);
   BLKS = (int *) malloc(nproc*ITEMSIZE);
   blk = (Num_Transactions+nproc-1)/nproc;
   for (i=0; i < nproc; i++) BLKS[i] = blk*i;
}	

double random_gen()
{
   double new_n;
   int seed;
   seed = random();
   new_n = (double) (seed)*1.0/2147483647.0;
   return new_n;
}

double exdpdef(double mean){
   
   double x;
   do {x=random_gen();}
   while (x == 0.0);
   
   return ((-log(x))*mean);
   
}

double gauss_normal(double variance, double mean){
   static int iset=0;
   static double gset;
   double fac,r,v1,v2;
   double sigma = sqrt(variance);
   if(iset==0){
      do{
         v1=2.0*random_gen()-1.0; 
         v2=2.0*random_gen()-1.0;
         r=v1*v1+v2*v2;
      } while(r>=1.0 || r==0.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac*sigma+mean;
      iset=1;
      return v2*fac*sigma+mean;
   }
   else{
      iset=0;
      return gset;
   }
}



double gammln(double xx){
   
   double x, tmp,ser;
   static double cof[6]={76.18009173, -86.50532033, 24.01409822, -1.231739516,
                         0.120858003e-2, -0.536382e-5};
   int j;
   
   x = xx-1.0;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.0;
   for (j=0;j<=5;j++){
      x +=1.0;
      ser += cof[j]/x;
   }
   return -tmp+log(2.50662827465*ser);
}


double poidev(double xm){
   
   static double sq, alxm,g,oldm=(-1.0);
   double em, t, y;
   
   if (xm < 12) {
      if (xm != oldm){
         oldm = xm;
         g = exp(-xm); 
      }
      em = -1;
      t = 1.0;
      do {
         em += 1.0;
         t *= random_gen();
      } while (t >g);
   }
   else {
      if (xm != oldm){
         oldm = xm;
         sq = sqrt(2.0*xm);
         alxm=log(xm);
         g=xm*alxm- gammln(xm+1.0);
      }
      do {
         do {
            y= tan(PI*random_gen());
            em = sq*y+xm;
         } while (em <0.0);
         em=floor(em);
         t= 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (random_gen() >t);
   }
   return em;
}

typedef struct itemsets{
   int *items;
   int size;
   double weight;
   double corrupt;
} Itemsets;

Itemsets *Tau;

void print_itemset(Itemsets *item, int id)
{
   int i;
   
   printf("%d : %d %f %f -- ", id, item->size, item->weight, item->corrupt);
   for(i=0; i < item->size; i++)
      printf("%d,", item->items[i]);
   printf("\n");
}

/* generate set Tau of set of potentially large itemsets bounded by L
   using correllation level etc. also includes normalizing them */
void Generate_Item_Sets(){
   int size, corr_items;
   int *bitvec;
   int i, item;
   int j, k;
   double sum, sum2;
   double corr;
   int rand_temp;

   bitvec = (int *)malloc(Number_Items*ITEMSIZE);
   Tau = (Itemsets*)malloc(sizeof(Itemsets)*Num_Large_Itemset);
   if(Tau== NULL) {
      perror("Insufficient Memory \n");
      exit(3);
   }
   while(!(size = (int)poidev((double)Average_LargeItemset_Size)));
   if (size > Number_Items) size =  Number_Items;
   Tau[0].items = (int*)malloc(ITEMSIZE*size);
   Tau[0].size = size;
   Tau[0].weight = exdpdef(1.0);
   Tau[0].corrupt = gauss_normal(0.1, 0.5);
   
   if(Tau[0].items == NULL) {
      perror("Insufficient Memory \n");
      exit(3);
   }
   for(i=0; i < Number_Items; i++) bitvec[i]= 0;
   for(j=0;j<size;j++){
      item = (int) (Number_Items*random_gen());
      if (bitvec[item] == 1) j--;
      else{
         Tau[0].items[j] = item;
         bitvec[item]=1;
      }
   }
   
/* Doubts about random generation.... */
   
   sum = Tau[0].weight;
   for(i=1; i< Num_Large_Itemset; i++){
      
      while(!(size = (int)poidev((double)Average_LargeItemset_Size)));
      Tau[i].items = (int*)malloc(ITEMSIZE*size);
      Tau[i].size = size;
      Tau[i].weight = exdpdef(1.0);
      Tau[i].corrupt = gauss_normal(0.1, 0.5);
      sum = sum +Tau[i].weight;
      if(Tau[i].items == NULL) {
         perror("Insufficient Memory \n");
         exit(3);
      }
      corr = exdpdef(Correllation_Level);
      corr_items = (int) (corr*Tau[i].size);
      if (corr_items > Tau[i].size) corr_items = Tau[i].size;
      if (corr_items > Tau[i-1].size) corr_items = Tau[i-1].size;
      for(k=0; k < Number_Items; k++) bitvec[k]= 0;
      for(j=0;j<corr_items;j++){
         rand_temp= (int) (Tau[i-1].size*random_gen());
         /*printf("rand_tmp %d\n", rand_temp);*/
         if (bitvec[Tau[i-1].items[rand_temp]] == 1) j--;
         else{
            Tau[i].items[j] = Tau[i-1].items[rand_temp];
            bitvec[Tau[i-1].items[rand_temp]] = 1;
         }
      }
      for(j=corr_items;j<size;j++){
         item = (int) (Number_Items*random_gen());
         if (bitvec[item] == 1) j--;
         else{
            Tau[i].items[j] = item;
            bitvec[item]=1;
         }
      }
   }

   sum2 = 0;
   for(i=0; i< Num_Large_Itemset; i++){
      Tau[i].weight = Tau[i].weight/sum;
      sum2 += Tau[i].size;
      //print_itemset(&Tau[i], i);
   }

   printf("AVERAGE ITEMSET SIZE = %f\n", sum2/Num_Large_Itemset);
   free(bitvec);
   return;
}



int get_item_set_no(){
   
   double x;
   int i;
   double sum;
   
   x= random_gen();
   /*printf("rand %f \n", x);*/
   i=0;
   sum = 0;
   while(i<Num_Large_Itemset){
      sum=sum+Tau[i].weight ;
      /*printf("sum %f \n", sum);*/
      if (sum>x){
         break;
      }
      i++;
   }
   if (i == Num_Large_Itemset) return i-1;
   return i;
}



int compare(const void *ptr1, const void *ptr2)
{
   int a, b;

   a = *(int *)ptr1;
   b = *(int *)ptr2;
   if (a < b) return -1;
   else if ( a > b) return 1;
   else return 0;
}


/* write this to a file, memory will definitely be a problem*/
/* Basically assign items to a transaction using progressive
   potentially large itemsets from Tau, including a corruption
   level c */	

int Generate_Transaction(int tid){
/* To be done */
   int size;
   double corrupt, rtmp;
   static int iset=0;
   static int gset;
   int i;
   int j;
   int num;
   int size_item; 
   int *trans;
   int *bitvec;

   trans = (int *) malloc(ITEMSIZE*(Number_Items+2));
   bitvec = (int *) malloc(ITEMSIZE*Number_Items);
   for(i=0; i < Number_Items; i++){
      trans[i+2] = Number_Items+10;
      bitvec[i] = 0;
   }

   while(!(size = (int)poidev((double)Average_Transaction_Size)));
   if (size > Number_Items) size = Number_Items;
   if (size >= distinct_items+max_itemlen)
      size = distinct_items+max_itemlen-1;
   i= 0;
   while(i<size){
      if (iset==1){
         num = gset;
         iset = 0;
      }
      else num= get_item_set_no();
      /*printf("num = %d\n", num);*/
      corrupt = Tau[num].corrupt;
      size_item = Tau[num].size;
      /*printf(" [--%d %d %d %d--] \n", i, i+size_item, size, num);*/
      if(i+size_item>size){
         if(random_gen() > 0.5 && i != 0){
            iset = 1;
            gset = num;
            break;
         }
      }
      for(j=0 ; j<size_item;j++){
         rtmp = random_gen();
         /*printf("[[%f %f %d %d]]\n", rtmp, corrupt,  Tau[num].items[j],
                bitvec[Tau[num].items[j]]);*/
         if(rtmp >= corrupt){
            if (bitvec[Tau[num].items[j]] == 0){
               trans[i+2] = Tau[num].items[j];
               bitvec[Tau[num].items[j]] = 1;
               i++;
               /*printf("***** %d\n", i);*/
            }
         }
      }
      /* printf("\n");*/
   }
   size = i;
   qsort(&trans[2], size, ITEMSIZE, compare);
   if (binary_output){
      trans[0] = tid;
      trans[1] = size;
      if (curr_proc+1 < nproc && tid == BLKS[curr_proc+1]){
         if (outbuf_pos > 0)
            write(fd[curr_proc], (char *)OUTBUF,outbuf_pos*ITEMSIZE);
         outbuf_pos = 0;
         curr_proc++;
      }
      if (outbuf_pos+size+3 >= OUTBUFSIZE){
         write(fd[curr_proc], (char *)OUTBUF,outbuf_pos*ITEMSIZE);
         outbuf_pos = 0;
      }
      
      /*printf("CURR PROC %d %d\n", curr_proc, tid);*/
      /*write(fd[curr_proc],(char *)trans, (size+3)*ITEMSIZE);*/
      /*memcpy((void *) (&OUTBUF[outbuf_pos]), (void *)trans, 
        (size+2)*ITEMSIZE);
      outbuf_pos += size+3;*/
      OUTBUF[outbuf_pos++] = tid;
      OUTBUF[outbuf_pos++] = tid;
      OUTBUF[outbuf_pos++] = size;
      for (i=0; i < size; i++)
         OUTBUF[outbuf_pos++] = trans[i+2];
   }
   else{
      fprintf(fp,"%d %d %d", tid, tid, size);
      for(i=0; i < size; i++)
         fprintf(fp," %d", trans[i+2]);
      fprintf(fp,"\n");
   }
   free(trans);
   free(bitvec);
   return size;
}



void get_distinct_items()
{
   int *bitvec;
   int i,j;
   int num = 0;
   int max = 0;
   
   bitvec = (int *) malloc(ITEMSIZE*Number_Items);
   for(i=0; i < Number_Items; i++)
      bitvec[i] = 0;
   for(i=0 ; i < Num_Large_Itemset; i++){
      if (Tau[i].size > max) max = Tau[i].size;
      for(j=0; j < Tau[i].size;j++)
         bitvec[Tau[i].items[j]] = 1;
   }
   
   for (i=0; i < Number_Items; i++)
      if (bitvec[i] == 1)
         num++;
   distinct_items = num;
   max_itemlen = max;
   printf("distinct items = %d, max_itemlen = %d\n",distinct_items, max_itemlen);
   free(bitvec);
}


int main(int argc, char **argv){
   int i,p;
   struct timeval tp;
   double sum;
   int *buf;
   char *filen;
   filen = (char *)malloc (256);      
   //gettimeofday(&tp,(struct timezone *)0);   
   //srandom(tp.tv_usec);
   //srandom(942046);
   srandom(randseed);
  
   parse_args(argc, argv);
   if (binary_output){
      sprintf(filen,"%s.conf", fname);
      fd[0] = open(filen,(O_CREAT|O_WRONLY),(S_IRWXU|S_IRWXG|S_IRWXO));
      if (fd[0] <0){
         printf("ERROR OPENING FILE CONF %s\n", filen);
         exit(3);
      }
      write(fd[0], (char *)&Num_Transactions, ITEMSIZE);
      write(fd[0], (char *)&Number_Items, ITEMSIZE);
      float tsz = (float) Average_Transaction_Size;
      write(fd[0], (char *)&tsz, sizeof(float));
      close(fd[0]);
      for (p = 0; p < nproc; p++){
         if (!fname) fd[p] = (int)stdout;
         else{
            if (nproc > 1) sprintf(filen,"%s.data.P%d", fname, p);
            else sprintf(filen,"%s.data", fname);
            fd[p] = open(filen,(O_CREAT|O_WRONLY),(S_IRWXU|S_IRWXG|S_IRWXO));
            /*fd[p] = open(filen, (O_WRONLY | O_CREAT | O_EXCL));*/
            if (fd[p] <0){
               printf("ERROR OPENING FILE %s\n", filen);
               exit(3);
            }
            else printf("OPENED %s, %d\n", filen, fd[p]);
         }
      }
   }
   else{
      if (!fname) fp = stdout;
      else{
         sprintf(filen,"%s.conf", fname);
         fp = fopen(filen, "w+");
         if (fp == NULL){
            printf("ERROR OPENING FILE %s\n", filen);
            exit(3);
         }
         fprintf(fp, "%d\n", Num_Transactions);
         fprintf(fp, "%d\n", Number_Items);
         fprintf(fp, "%f\n", (float)Average_Transaction_Size);
         fclose(fp);
         
         sprintf(filen,"%s.data", fname);
         fp = fopen(filen, "w+");
         if (fp == NULL){
            printf("ERROR OPENING FILE %s\n", fname);
            exit(3);
         }
      }
   }

   free(filen);

/* initialization function not needed as parseargs does basic initialization
   init();*/
   
   
   Generate_Item_Sets();
   printf("seed = %d\n", randseed);
   get_distinct_items();
   i= 0;
   sum = 0;
   while(i < Num_Transactions){      
      sum += Generate_Transaction(i);
      i++;
      
   }
   printf("AVERAGE TRANSACTION  SIZE = %f\n", sum/i);
   //printf("seed = %d\n", tp.tv_usec);
   if (binary_output){
      if (outbuf_pos > 0)
         write(fd[curr_proc], (char *)OUTBUF,outbuf_pos*ITEMSIZE);
      for (i=0; i < nproc; i++)
         close(fd[i]);
   }
   else fclose(fp);
   exit(0);
}





