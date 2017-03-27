#include <iostream>
#include <stdlib.h>

#include <visp/vpDebug.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageFilter.h>
#include <visp/vpDisplayX.h>


using namespace std ;


vpImage<float> winnerTakesAll(vpImage<unsigned char> left, vpImage<unsigned char> right)
{
  vpImage<float> disparity(288,384,0);

  int actual, read, temp;
  for (int i = 0; i < left.getRows(); ++i)
  {
    for (int j = 0; j < left.getCols(); ++j)     //boucle pixel
    {
      actual = 9000;
      for (int k = 0; k < right.getCols(); ++k)  //boucle ligne
      {
        read = abs(right[i][k] - left[i][j]);
        if(read < actual)
        {
          actual = read;
          temp = abs(j-k);
        }
      }
      disparity[i][j]=temp;
    }
  }
  return disparity;
}



vpImage<float> SSD(vpImage<unsigned char> left, vpImage<unsigned char> right, int size, double* kernel)
{
  vpImage<float> disparity(288,384,0);

  int actual, read, temp, temp_i, temp_j;
  for (int i = size/2; i < left.getRows()-size/2; ++i)       //boucle pixel
  {
    std::cout<< i << std::endl;
    for (int j = size/2; j < left.getCols()-size/2; ++j)     
    {
      actual = 9000;
      for (int k = size/2; k < right.getCols()-size/2; ++k)  //boucle ligne
      {
        read = 0;
        temp_i=0;
        temp_j=0;

        for(int ii = i-size/2; ii<i+size; ++ii)              //calcul de la SSD
        {
          for (int jj = j-size/2; jj < j+size; ++jj)
          {
            read += kernel[temp_i*size+temp_j]*(right[i][k] - left[i][j]) * (right[i][k] - left[i][j]);
            temp_j++;
          }
          temp_i++;
        }

        if(read < actual)
        {
          actual = read;
          temp = abs(j-k);
        }
      }
      disparity[i][j]=temp;
    }
  }
  return disparity;
}



int main()
{
  int size=3;
  vpImage<unsigned char> Ileft(288,384);
  vpImage<unsigned char> Iright(288,384);

  vpImage<float> disparity(288,384,0);
  vpImage<unsigned char> disparity_UC(288,384,0);

  double kernel[size];
 
  vpImageIo::read(Ileft,"../Data/scene_l.pgm") ;
  vpImageIo::read(Iright,"../Data/scene_r.pgm") ;
  
  disparity = winnerTakesAll(Ileft, Iright);

  
  // vpImageFilter::getGaussianKernel(kernel, size);
  // disparity = SSD(Ileft, Iright, size, kernel);
  
  vpImageConvert::convert(disparity,disparity_UC);

  vpDisplayX d1(disparity_UC);
  vpDisplay::display(disparity_UC) ;
  vpDisplay::flush(disparity_UC) ;
  vpDisplay::getClick(disparity_UC) ;

  vpImageIo::write(disparity_UC, "../result/disparity.pgm");


  return 0;
}
