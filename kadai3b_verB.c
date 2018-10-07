#include <stdio.h>
#include <math.h>

#define epsilon 1e-8//逆行列の近似の許容誤差
#define delta 0.50

double sum_Mat(double X[][5]);//行列の行の和の最大値を返す
void mul_Mat(double[][5],double[][5]);//後ろの行列と前の行列の積を後ろの行列に代入する
void recur_Mat(double[][5],double[][5]);//課題2(e)の再帰式で定められた計算をする

int main(void){
  int n,i,j,loop;
  double mat_X[5][5],mat_Y[5][5],sum=0,gamma[5],nu=0,n_Eps,n_Chil,alpha=0;
  double mat_S[5][5] = {{0,1/2.0,0,1/2.0,0,},
			{0,0,1,0,0},
			{1/2.0,0,0,0,1/2.0},
			{0,0,0,0,1},
			{1/5.0,1/5.0,1/5.0,1/5.0,1/5.0}};//確率行列S

  for(loop=0;loop<9999;loop++){
    alpha += 0.0001;

    nu = 0;

    for(i=0;i<5;i++){
      for(j=0;j<5;j++){
	if(i==j){
	  mat_Y[i][j]=1;
	}else{
	  mat_Y[i][j]=0;
	}
	mat_X[i][j]= mat_S[i][j] * alpha;
      }
    }//Y=I,X=αSを代入

    while(1){
      if(sum_Mat(mat_Y) <= delta){
	break;
      }
      mul_Mat(mat_X,mat_Y);
      nu++;
    }//Yの行和がδを下回るまでYにXをかけて、その回数を記録する

    n_Eps = ceil(nu * log (epsilon * delta * (1-delta) / nu) / log(delta));

    n_Chil = ceil(log2(n_Eps + 1));

    for(i=0;i<5;i++){
      for(j=0;j<5;j++){
	if(i==j){
	  mat_Y[i][j]=1;
	}else{
	  mat_Y[i][j]=0;
	}
      }
    }//Y=Iを代入

    for(n=0;n<n_Chil;n++){
      recur_Mat(mat_X,mat_Y);
    }//V_N~(ε)を計算((I-F)^-1の近似)

    for(i=0;i<5;i++){
      gamma[i] = 0;
      for(j=0;j<5;j++){
	gamma[i] += mat_Y[j][i];
      }
      gamma[i] *= (1 - alpha) / 5.0;
    }//γの定義に従ってγを求める


  }

  return 0;
}

double sum_Mat(double X[][5]){
  int i,j;
  double max=0,sum;

  for(i=0;i<5;i++){
    sum=0;
    for(j=0;j<5;j++){
      sum+= X[i][j];
    }
    if(max < sum){
      max = sum;
    }
  }//Xの行和の最大値を求める

  return max;
}

void mul_Mat(double X[][5],double Y[][5]){
  int i,j,k;
  double memo[5][5];

  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      memo[i][j] = Y[i][j];
    }
  }//memoにYを記録

  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      Y[i][j]=0;
      for(k=0;k<5;k++){
	Y[i][j]+= X[i][k] * memo[k][j];
      }
    }
  }//YにXYを代入する
}

void recur_Mat(double X[][5],double Y[][5]){
  int i,j,k;
  double memo_X[5][5],memo_Y[5][5],i_Mat[5][5];

  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      memo_X[i][j] = X[i][j];
      memo_Y[i][j] = Y[i][j];
      if(i==j){
	i_Mat[i][j] = 1;
      }else{
	i_Mat[i][j] = 0;
      }
    }
  }//memo_XにX,memo_YにYを代入、単位行列の生成

  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      X[i][j] = 0;
      Y[i][j] = 0;
      for(k=0;k<5;k++){
	X[i][j] += memo_X[i][k] * memo_X[k][j];
	Y[i][j] += (i_Mat[i][k] + memo_X[i][k]) * memo_Y[k][j];
      }
    }
  }//Xにmemo_Xの二乗、Yに(I+memo_X)とmemo_Yの積を代入する
}
