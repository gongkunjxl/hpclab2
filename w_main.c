/**
 * 预处理文件: data.txt  第1: int nSize(矩阵中非零个数)
 *N个int值表示每行中第一个数值的偏移 nSize个列 nSize个double类型的值
 *v1_b.txt 下面N个double类型的值
 *v1_x0.txt下面N个double类型的值
 *
**/
#include <stdio.h>
#include <stdlib.h>

#define NX 360
#define NY 180
#define NZ 38
#define N NX*NY*NZ

void handleA()
{
    FILE *fp=fopen("A.txt","r");
    if(fp==NULL){
        printf("can not open file!\n");
        exit(0);
    }
    int i,j,k,num;
    double data;

    double *AV=(double *)malloc(N*19* sizeof(double));
    int *AI=(int*)malloc(N*sizeof(int));
    int *AJ=(int*)malloc(N*19*sizeof(int));
    int ind=0;  
    int m_c;    

    int m_i,m_j ;
    for(m_i=0;m_i<N;m_i++) {
        
        AI[m_i]=ind;
        for(m_j=0;m_j<19;m_j++)
        {
            fscanf(fp,"%d",&i); fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&k); fscanf(fp,"%d",&num);
            fscanf(fp,"%lf",&data);
           /* char buf[100];
            fscanf(fp,"%s",&buf);
            printf("---%s\n",buf);
            sscanf(buf,"%.30lf",&data);
	    printf("---->%.35lf\n",data);
           data=1.0;*/
            switch (m_j)
            {
                case 0:
                {
                    m_c=i*NY*NZ+j*NZ+k;
                    AV[ind]=data;    AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 1:
                {
                    int tmp_i=i-1;
                    if(i==0) tmp_i=(tmp_i+NX)%NX;
                    m_c=tmp_i*NY*NZ+j*NZ+k;
                    AV[ind]=data;    AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 2:
                {
                    int tmp_i=i+1;
                    if(i==NX-1) tmp_i=(tmp_i+NX)%NX;;
                    m_c=tmp_i*NY*NZ+j*NZ+k;
                    AV[ind]=data;   AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 3:
                {
                    int tmp_j=j-1; int tmp_i=i;
                    if(j==0){  tmp_i=(tmp_i+NX/2)%NX; tmp_j=j; }
                    m_c=tmp_i*NY*NZ+tmp_j*NZ+k;
                    AV[ind]=data;    AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 4:
                {
                    int tmp_i=i;   int tmp_j=j+1;
                    if(j==NY-1){ tmp_i=(tmp_i+NX/2)%NX; tmp_j=j;}
                    m_c=tmp_i*NY*NZ+tmp_j*NZ+k;
                    AV[ind]=data;    AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 5:
                {
                    int tmp_i=i+1;  int tmp_j=j+1;
                    if(i==NX-1) tmp_i=(tmp_i+NX)%NX;
                    if(j==NY-1){
                        tmp_i=(tmp_i+NX/2)%NX;
                        tmp_j=j;
                    }
                    m_c=tmp_i*NY*NZ+tmp_j*NZ+k;
                    AV[ind]=data;   AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 6:
                {
                    int tmp_i=i+1;  int tmp_j=j-1;
                    if(i==NX-1) tmp_i=(tmp_i+NX)%NX;
                    if(j==0){
                        tmp_i=(tmp_i+NX/2)%NX;
                        tmp_j=j;
                    }
                    m_c=tmp_i*NY*NZ+tmp_j*NZ+k;
                    AV[ind]=data;   AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 7:
                {
                    int tmp_i=i-1;  int tmp_j=j-1;
                    if(i==0) tmp_i=(tmp_i+NX)%NX;
                    if(j==0){
                        tmp_i=(tmp_i+NX/2)%NX;
                        tmp_j=j;
                    }
                    m_c=tmp_i*NY*NZ+tmp_j*NZ+k;
                    AV[ind]=data;     AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 8:
                {
                    int tmp_i=i-1;  int tmp_j=j+1;
                    if(i==0) tmp_i=(tmp_i+NX)%NX;
                    if(j==NY-1){
                        tmp_i=(tmp_i+NX/2)%NX;
                        tmp_j=j;
                    }
                    m_c=tmp_i*NY*NZ+tmp_j*NZ+k;
                    AV[ind]=data;   AJ[ind]=m_c;
                    ind++;
                    break;
                }
                case 9:
                {
                    int tmp_k=k-1;
                    if(k==0) break;
                    else{
                        m_c=i*NY*NZ+j*NZ+tmp_k;
                        AV[ind]=data;    AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 10:
                {
                    int tmp_i=i-1;  int tmp_k=k-1;
                    if(i==0) tmp_i=(tmp_i+NX)%NX;
                    if(k==0) break;
                    else{
                        m_c=tmp_i*NY*NZ+j*NZ+tmp_k;
                        AV[ind]=data;   AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 11:
                {
                    int tmp_i=i+1;  int tmp_k=k-1;
                    if(i==NX-1) tmp_i=(tmp_i+NX)%NX;
                    if(k==0) break;
                    else{
                        m_c=tmp_i*NY*NZ+j*NZ+tmp_k;
                        AV[ind]=data;   AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 12:
                {
                    int tmp_i=i;  int tmp_j=j-1; int tmp_k=k-1;
                    if(j==0){ tmp_i=(tmp_i+NX/2)%NX; tmp_j=j;}
                    if(k==0) break;
                    else{
                        m_c=tmp_i*NY*NZ+tmp_j*NZ+tmp_k;
                        AV[ind]=data;     AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 13:
                {
                    int tmp_i=i;  int tmp_j=j+1; int tmp_k=k-1;
                    if(j==NY-1){ tmp_i=(tmp_i+NX/2)%NX; tmp_j=j;}
                    if(k==0) break;
                    else{
                        m_c=tmp_i*NY*NZ+tmp_j*NZ+tmp_k;
                        AV[ind]=data;   AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 14:
                {
                    int tmp_k=k+1;
                    if(k==NZ-1) break;
                    else{
                        m_c=i*NY*NZ+j*NZ+tmp_k;
                        AV[ind]=data;   AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 15:
                {
                    int tmp_i=i-1;  int tmp_k=k+1;
                    if(i==0){ tmp_i=(tmp_i+NX)%NX;}
                    if(k==NZ-1) break;
                    else{
                        m_c=tmp_i*NY*NZ+j*NZ+tmp_k;
                        AV[ind]=data;   AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 16:
                {
                    int tmp_i=i+1;   int tmp_k=k+1;
                    if(i==NX-1){ tmp_i=(tmp_i+NX)%NX;}
                    if(k==NZ-1) break;
                    else{
                        m_c=tmp_i*NY*NZ+j*NZ+tmp_k;
                        AV[ind]=data;    AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 17:
                {
                    int tmp_i=i;  int tmp_j=j-1; int tmp_k=k+1;
                    if(j==0){ tmp_i=(tmp_i+NX/2)%NX; tmp_j=j;}
                    if(k==NZ-1) break;
                    else{
                        m_c=tmp_i*NY*NZ+tmp_j*NZ+tmp_k;
                        AV[ind]=data;   AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                case 18:
                {
                    int tmp_i=i;  int tmp_j=j+1; int tmp_k=k+1;
                    if(j==NY-1){ tmp_i=(tmp_i+NX/2)%NX; tmp_j=j;}
                    if(k==NZ-1) break;
                    else{
                        m_c=tmp_i*NY*NZ+tmp_j*NZ+tmp_k;
                        AV[ind]=data;AJ[ind]=m_c;
                        ind++;
                    }
                    break;
                }
                default:
                    break;
            }
        }
 //       printf("----all data--->%d\n",AI[m_i]);

    }
    fclose(fp);
    FILE *fw=fopen("data.txt","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    int size=ind;
    fwrite(&size,4,1,fw);
    fwrite(AI,4,N,fw);
    fwrite(AJ,4,size,fw);
    fwrite(AV,8,size,fw);

    fclose(fw);
    free(AI);	free(AV); free(AJ);
}
void handleB()
{
    FILE *fp=fopen("b.txt","r");
    if(fp==NULL){
        printf("can not open b.txt file!\n");
        exit(0);
    }
    int i,j,k;
    double value;
    double *data=(double*)malloc(N*sizeof(double));
    int m_i ;
    for(m_i=0;m_i<N;m_i++) {
        fscanf(fp, "%d", &i);   fscanf(fp, "%d", &j);
        fscanf(fp, "%d", &k);   fscanf(fp, "%lf", &value);
        data[m_i]=value;
     //   printf("(%d, %d, %d)--->%.30lf\n",i,j,k,value);
    }
    fclose(fp);

    FILE *fw=fopen("v1_b.txt","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    fwrite(data,8,N,fw);
    fclose(fw);
    free(data);
}

void handleX()
{
    FILE *fp=fopen("x0.txt","r");
    if(fp==NULL){
        printf("can not open b.txt file!\n");
        exit(0);
    }
    int i,j,k;
    double value;
    double *data=(double*)malloc(N*sizeof(double));
    int m_i ;
    for(m_i=0;m_i<N;m_i++) {
        fscanf(fp, "%d", &i);   fscanf(fp, "%d", &j);
        fscanf(fp, "%d", &k);   fscanf(fp, "%lf", &value);
        data[m_i]=value;
        //printf("(%d, %d, %d)--->%.30lf\n",i,j,k,value);
    }
    fclose(fp);

    FILE *fw=fopen("v1_x0.txt","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    fwrite(data,8,N,fw);
    fclose(fw);
    free(data);
}

int main(int argc, char **argv) {
   printf("start A.txt\n"); 
   handleA();
   printf("end A.txt\n"); 
   printf("start b.txt\n");
   handleB();
   printf("end  b.txt\n"); 
   printf("start x0.txt\n"); 
   handleX();
   printf("end x0.txt\n"); 
    return 0;
}


































