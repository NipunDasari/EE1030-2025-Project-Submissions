#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void normalize(double *v,int n){
	double s=0; for(int i=0;i<n;i++) s+=v[i]*v[i];
	s=sqrt(s); if(s<1e-12)s=1e-12;
	//to avoid dividing by zero
	for(int i=0;i<n;i++) v[i]/=s;}
int main(){
	char in[100], out[100], magic[2];
	int r,w,h,maxv;
	printf("Input PGM filename: "); scanf("%s",in);
	printf("Rank r: "); scanf("%d",&r);
	printf("Output PGM filename: "); scanf("%s",out);
	FILE *f=fopen(in,"rb");
	if(!f){printf("Cannot open input file.\n"); return 1;}
	fscanf(f,"%2s %d %d %d",magic,&w,&h,&maxv);
	fgetc(f);
	int m=h,n=w,N=m*n;
	unsigned char *img=malloc(N);+
	fread(img,1,N,f);
	fclose(f);//since sizeof(double) is 8
	double *A=malloc(8*N);
	for(int i=0;i<N;i++) A[i]=img[i];
	double *U=malloc(m*r*8), *S=malloc(r*8), *V=malloc(n*r*8), *v=malloc(n*8), *Av=malloc(m*8), *AtAv=malloc(n*8);
	for(int col=0;col<r;col++){
		for(int i=0;i<n;i++) v[i]=((double)rand()/RAND_MAX)-0.5;
		for(int it=0;it<25;it++){
			for(int i=0;i<m;i++){
				double s=0;//matrix multiplication
				for(int j=0;j<n;j++) s+=A[i*n+j]*v[j];
				Av[i]=s;}
			for(int j=0;j<n;j++){
				double s=0;
				for(int i=0;i<m;i++) s+=A[i*n+j]*Av[i];
				AtAv[j]=s;}
			for(int p=0;p<col;p++){
				double d=0;
				for(int i=0;i<n;i++) d+=AtAv[i]*V[i+p*n];
				for(int i=0;i<n;i++) AtAv[i]-=d*V[i+p*n];}
			normalize(AtAv,n);//a unit basis vector for V matrix
			for(int i=0;i<n;i++) v[i]=AtAv[i];}
		for(int i=0;i<n;i++) V[i+col*n]=v[i];
		double sigma_sq=0;
		for(int i=0;i<m;i++){
			double s=0;
			for(int j=0;j<n;j++) s+=A[i*n+j]*v[j];
			U[i+col*m]=s;
			sigma_sq+=s*s;}
		S[col]=sqrt(sigma_sq);
		for(int i=0;i<m;i++) U[i+col*m]/=S[col];}
	double *R=malloc(8*N);
	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++){
		double s=0;
		for(int t=0;t<r;t++) s+=U[i+t*m]*S[t]*V[j+t*n];
		R[i*n+j]=s;}
	for(int i=0;i<N;i++){
		double v=R[i]; if(v<0)v=0; if(v>255)v=255; 
		img[i]=(unsigned char)(v+0.5);}//to approximate to better values and avoid smaller values
	FILE *g=fopen(out,"wb");
	fprintf(g,"P5\n%d %d\n255\n",w,h);
	fwrite(img,1,N,g);
	fclose(g);
	double sum=0.0;
	for(int i=0;i<N;i++) 
      sum+=(A[i]-R[i])*(A[i]-R[i]); sum = sqrt(sum); printf("%lf\n", sum);
	return 0;}

