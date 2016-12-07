/**
 * 预处理文件: data.txt  第1: int nSize(矩阵中非零个数)
 *N个int值表示每行中第一个数值的偏移 nSize个列 nSize个double类型的值
 *v1_b.txt 下面N个double类型的值
 *v1_x0.txt下面N个double类型的值
 *
**/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 360
#define NY 180
#define NZ 38
#define N NX*NY*NZ
double *AV,*b,*x0,*diag;
int *AI,*AJ;
int nSize;


int InitAll()
{
    FILE *fp1=fopen("data.txt","rb");
    if(fp1==NULL){
        printf("can't open read file\n");
        exit(0);
    }
    int size;
    fread(&size,4,1,fp1);
    AV=(double *)malloc(size* sizeof(double));
    AI=(int*)malloc(N*sizeof(int));
    AJ=(int*)malloc(size* sizeof(int));
    fread(AI,4,N,fp1);
    fread(AJ,4,size,fp1);
    fread(AV,8,size,fp1);
    fclose(fp1);

    FILE *fp2=fopen("v1_b.txt","rb");
    if(fp2==NULL){
        printf("can't open read file\n");
        exit(0);
    }
    b=(double *)malloc(N* sizeof(double));
    fread(b,8,N,fp2);
    fclose(fp2);

    FILE *fp3=fopen("v1_x0.txt","rb");
    if(fp3==NULL){
        printf("can't open read file\n");
        exit(0);
    }
    x0=(double *)malloc(N* sizeof(double));
    fread(x0,8,N,fp3);
    fclose(fp3);
    return size;
}

void getAmVector(double *result,double *data) //矩阵向量相乘
{
    int ind=0;
    int i,j;
    int size;
    double tmp;
    for(i=0;i<N-1;i++) //N-1
    {
        ind=AI[i];
        size=AI[i+1]-AI[i];
        tmp=0.0;
        for(j=0;j<size;j++){
            tmp=tmp+AV[ind+j]*data[AJ[ind+j]];
        }
        result[i]=tmp;
    }
    size=nSize-AI[N-1];
    tmp=0.0;
    ind=AI[N-1];
    for(j=0;j<size;j++){
        tmp=tmp+AV[ind+j]*data[AJ[ind+j]];
    }
    result[N-1]=tmp;
}

double checkSum(double *data)
{
    int ind=0;
    int i,j;
    int size;
    double tmp;
    double sum=0.0;
    for(i=0;i<N-1;i++) //N-1行
    {
        ind=AI[i];
        size=AI[i+1]-AI[i];
        tmp=0.0;
        for(j=0;j<size;j++){
            tmp=tmp+AV[ind+j]*data[AJ[ind+j]];
        }
        sum=sum+(tmp-b[i])*(tmp-b[i]);
    }

    size=nSize-AI[N-1];
    tmp=0.0;
    ind=AI[N-1];
    for(j=0;j<size;j++){
        tmp=tmp+AV[ind+j]*data[AJ[ind+j]];
    }
    sum=sum+(tmp-b[N-1])*(tmp-b[N-1]);

    return sum;

}
double VectorDianJi(double *dataA,double *dataB)
{
    int i;
    double sum=0.0;
    for(i=0;i<N;i++){
        sum=sum+dataA[i]*dataB[i];
    }
    return sum;
}


void JacobiCal(double *Rm,double *R) //矩阵Rm=M^(-1)*R Jacobi迭代
{
    int ind=0;
    int i,j,size;
    double tmp,m_rst;
    for(i=0;i<N-1;i++){
        ind=AI[i];
        tmp=0.0;    m_rst=0.0;
        size=AI[i+1]-AI[i];
        for(j=0;j<size;j++){
            if(AJ[ind+j]==i){  //对角线上的元素
                tmp=AV[ind+j];
            }else{
                m_rst=m_rst+AV[ind+j]*R[AJ[ind+j]];
            }
        }
        Rm[i]= (R[i]-m_rst)/tmp;
    }
    size=nSize-AI[N-1];
    ind=AI[N-1];
    tmp=0.0; m_rst=0.0;
    for(j=0;j<size;j++){
        if(AJ[ind+j]==(N-1)){  //对角线上的元素
            tmp=AV[ind+j];
        }else{
            m_rst=m_rst+AV[ind+j]*R[AJ[ind+j]];
        }
    }
    Rm[N-1]=(R[N-1]-m_rst)/tmp;
}

double getjk(int j,int k)//获取Ajk值
{
  
   int i,size,ind;
   ind=0; size=0;
   if(j==(N-1)) size=nSize-AI[j];
   else{
     size=AI[j+1]-AI[j];	
   }
   ind=AI[j];
   double reVal=0.0;
   for(i=0;i<size;i++)
   {
      if(AJ[ind+i]==k){
         reVal=AV[ind+i];
         break;     
       }      
   }
   return reVal;
}

void factorDLU(double *diag)
{
    int ind,size;
    int i,j;
    for(i=0;i<N;i++)  //get diag of matrix A
    {
       if(i==(N-1)){size=nSize-AI[i];}
       else{
          size=AI[i+1]-AI[i];
       }
       ind=AI[i];
       for(j=0;j<size;j++){
          if(AJ[ind+j]==i){
             diag[i]=AV[ind+j];  //对角元素
             break;
          }
       }
    }
    for(i=0;i<N;i++)
    {
       if(i==(N-1)){size=nSize-AI[i];}
       else{
          size=AI[i+1]-AI[i];
       }
       ind=AI[i];
       for(j=0;j<size;j++){
          if(AJ[ind+j]==i){
             diag[i]=1.0/diag[i];  //对角元素
             break;
          }
       }

       for(j=0;j<size;j++){
         if(AJ[ind+j]>i){
           double tmp=getjk(AJ[ind+j],i);
	   if(tmp!=0){
	     diag[AJ[ind+j]]=diag[AJ[ind+j]]-tmp*diag[i]*AV[ind+j];	  	
           }    
	}      
      }
   } 
}
void getDLUCal(double *Rm,double *R)
{
   int i,j;
   int ind,size;
   double *z=(double*)malloc(N*sizeof(double)); //用于存储对角元素 
   double tmp;
   for(i=0;i<N;i++)
   {
       if(i==(N-1)){size=nSize-AI[i];}
       else{
          size=AI[i+1]-AI[i];
       }
       ind=AI[i];
       tmp=0.0;
       for(j=0;j<size;j++){
         if(AJ[ind+j]<i){
            tmp=tmp+AV[ind+j]*z[AJ[ind+j]];
         }
       }
      z[i]=R[i]-tmp*diag[i];
   }
   for(i=N-1;i>=0;i--)
   {
       if(i==(N-1)){size=nSize-AI[i];}
       else{
          size=AI[i+1]-AI[i];
       }
       ind=AI[i];
       tmp=0.0;
       for(j=0;j<size;j++){
         if(AJ[ind+j]>i){
            tmp=tmp+AV[ind+j]*Rm[AJ[ind+j]];
         }
       }
       Rm[i]=z[i]-tmp*diag[i];
   }
    free(z); 
}




void calGCR()
{
    int k=5;
    int i,j;
    double *R=(double*)malloc(N*sizeof(double));
    double *Rm=(double*)malloc(N*sizeof(double));
    double *p=(double*)malloc(N*sizeof(double));
    double *Ap=(double*)malloc(N* sizeof(double));
    double *ARm=(double*)malloc(N* sizeof(double));
    double *Apj=(double*)malloc(N* sizeof(double));

    double *pj[5];
    double betaij[5];
    for(i=0;i<5;i++){
        pj[i]=(double*)malloc(N*sizeof(double));
        betaij[i]=0;
    }
    getAmVector(R,x0);
    for(i=0;i<N;i++){   //get the R=b-Ax0
        R[i]=b[i]-R[i];
    }

    getDLUCal(Rm,R);    
    for(i=0;i<N;i++){
	p[i]=Rm[i];
    }
    double ai,checkVal;
    int ind=0,m=0;
    i=1;
    do{
        printf("Jacobi times------>%d\n",i);
        getAmVector(Ap,p);   //get Api-1 向量 
        
        for(j=0;j<N;j++){
          pj[(i-1)%k][j]=p[j];      //存储5个子空间向量 
        }

        ai=VectorDianJi(R,Ap)/VectorDianJi(Ap,Ap);
        //printf("ai---->%.30lf\n",ai);
        for(j=0;j<N;j++){       //xi=x(i-1)+a(i-1)*p(i-1)
            x0[j]=x0[j]+ai*p[j];
            R[j]=R[j]-ai*Ap[j];
       }
 //       JacobiCal(Rm,R);        //预处理
        getDLUCal(Rm,R);    
        ind=0;
        for(j=((i-1)/k)*k;j<i;j++){
            getAmVector(ARm,Rm);
            getAmVector(Apj,pj[ind]);
            betaij[ind]=-VectorDianJi(ARm,Apj)/VectorDianJi(Apj,Apj);
            ind++;
        }

       for(j=0;j<N;j++){
          p[j]=Rm[j];
       }     
        ind=0;
        for(j=((i-1)/k)*k;j<i;j++){
            for(m=0;m<N;m++){
                p[m]=p[m]+betaij[ind]*pj[ind][m];
            }    
            ind++;
        }
       checkVal=checkSum(x0);
       printf("checkSum---->%.30lf\n",checkVal); 
       i++;
    }while(checkVal>0.0000000001);
    free(R);free(Rm);free(p);free(Ap);free(ARm); free(Apj);
    for(i=0;i<5;i++){
       free(pj[i]);
    }
}



double checkData()
{
    int ind=0;
    int i,j;
    int size;
    double tmp;
    double sum=0.0;
    for(i=0;i<N-1;i++) //N-1行
    {
        ind=AI[i];
        size=AI[i+1]-AI[i];
        tmp=0.0;
        for(j=0;j<size;j++){
            tmp=tmp+AV[ind+j]*x0[AJ[ind+j]];
        }
        sum=sum+(tmp-b[i])*(tmp-b[i]);
    }

    size=nSize-AI[N-1];
    //printf("%d\n",AI[N-1]);
    tmp=0.0;
    ind=AI[N-1];
    for(j=0;j<size;j++){
        tmp=tmp+AV[ind+j]*x0[AJ[ind+j]];
    }
    sum=sum+(tmp-b[N-1])*(tmp-b[N-1]);

    return sum;
}



int main(int argc, char **argv) {
    AI=NULL;  AJ=NULL; AV=NULL; b=NULL; x0=NULL;
    printf(" start Init all data\n");
    nSize=InitAll();
    printf(" end Init all data----N: %d   nSize: %d\n",N,nSize);
//    double checkValue=checkData();
 //   printf("checksum ---->%.30lf\n",checkValue);

     diag=(double*)malloc(N*sizeof(double)); //用于存储对角元素    
     factorDLU(diag);
     calGCR();

    free(AI); free(AJ);free(AV);free(b);free(x0);
    free(diag);
    return 0;
}

































