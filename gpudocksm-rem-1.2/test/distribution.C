#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926
#define ABS(x)  ( ( (x) < 0) ? -(x) : (x) )


// // #define NORM(normed_x, scale)          1/((scale) * (sqrt(2 * (PI)))) * (exp(0.f - (normed_x)*(normed_x)/2))

// #define CAUCHY(normed_x, scale)        1/((PI) * (scale)) * (1/(1 + (normed_x) * (normed_x)))

// #define EXPONENTIAL(normed_x)          exp(0.f - (normed_x))
// #define LOGISTIC(normed_x, scale)      1/(scale) * ((EXPONENTIAL(normed_x)) / (pow(1 + (EXPONENTIAL(normed_x)), 2,0)))



float 
normPdf(float x, float loc, float scale){
  float norm_para, prob, pdf_val;

  norm_para = 1/(scale * sqrt(2 * PI));
  prob = exp (0.f - (x - loc) * (x - loc) / (2 * scale * scale));
  
  pdf_val = norm_para * prob;
  
  return pdf_val;
}

float
cauchyPdf(float x, float loc, float scale){
  float norm_para, prob, pdf_val;
  
  norm_para = 1/(PI * scale);
  prob = 1/(1 + ( (x-loc)/scale) * ((x-loc)/scale) );
  
  pdf_val = norm_para * prob;
  
  return pdf_val;
}

float
logisticPdf(float x, float loc, float scale){
  float norm_para, e_power, prob, pdf_val;
  
  norm_para = 1/scale;
  e_power = exp( -(x - loc) / scale);
  prob = e_power / pow(1+e_power, 2.0);
  
  pdf_val = norm_para * prob;
  
  return pdf_val;
}

float
waldPdf(float x, float loc, float scale){
  float norm_para, prob, pdf_val;
  
  float normed_x = (x - loc)/scale;

  norm_para = 1/(sqrt( 2 * PI * pow(normed_x, 3.0) ) * scale);
  prob = exp(-pow(normed_x-1, 2)/(2*normed_x));

  if (normed_x < 0) 
    pdf_val = 0.0f;
  else
    pdf_val = norm_para * prob;

  return pdf_val;
}

float 
laplacePdf(float x, float loc, float scale){
  float normed_x, pdf_val;

  normed_x = fabs(x-loc) / scale;

  pdf_val = (1/(2 * scale)) * exp(- normed_x);
  
  return pdf_val;
}

int main(void){
  float loc = -1.0041f;
  float scale = 0.11517f;
  float x = -1.3f;

  // float pdf_val = normPdf(x, loc, scale);
  // float pdf_val = cauchyPdf(x, loc, scale);
  // float pdf_val = logisticPdf(x, loc, scale);
  // float pdf_val = waldPdf(x, loc, scale);
  // float pdf_val = laplacePdf(x, loc, scale);
  
  // printf("pdf_val: %2.6f\n", pdf_val);
  
  return 0;
}
