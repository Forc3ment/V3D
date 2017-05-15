#include <iostream>
#include <stdlib.h>

#include <visp/vpDebug.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageFilter.h>
#include <visp/vpDisplayX.h>


using namespace std ;

void gaussianKernel(const int & taille, vpImage<double> kernel){

  int demi_taille=(taille+1)/2;

  double filter[demi_taille];

  vpImageFilter::getGaussianKernel(filter, taille);

  for(int i=0; i<taille; i++){
    for(int j=0; j<taille; j++){
      int ii=i-taille/2, jj=j-taille/2;
      kernel[i][j]=filter[abs(ii)]*filter[abs(jj)];
    }
  }
}

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
          //disparity[i][j]=k;
          disparity[i][j] = abs(j-k);
        }
      }
    }
  }
  return disparity;
}



vpImage<float> SSD(vpImage<unsigned char> left, vpImage<unsigned char> right, int size, vpImage<double> kernel)
{
  vpImage<float> disparity(288,384,0);

  double actual, read;
  int temp, temp_i, temp_j;
  for (int i = size/2; i < left.getRows()-size/2; ++i)       //boucle pixel
  {
    // std::cout<< i << std::endl;
    for (int j = size/2; j < left.getCols()-size/2; ++j)     
    {
      actual = 10000000;
      for (int k = size/2; k < right.getCols()-size/2; ++k)  //boucle ligne
      {

        read = 0;
        temp_i=0;
        temp_j=0;

        for(int ii = -size/2; ii<size/2; ++ii)              //calcul de la SSD
        {
          for (int jj = -size/2; jj < size/2; ++jj)
          {
            //std::cout << i << " " << j << std::endl;
            read += abs(kernel[temp_i][temp_j]*((right[i+ii][k+jj] - left[i+ii][j+jj]) * (right[i+ii][k+jj] - left[i+ii][j+jj])));
            //read += ((right[i+ii][k+jj] - left[i+ii][j+jj]) * (right[i+ii][k+jj] - left[i+ii][j+jj]));
            temp_j++;
          }
          temp_i++;
        }

        //std::cout<<read<<std::endl;
        if(read < actual)
        {
          actual = read;
          //temp = k;
          temp = abs(j-k);
        }
      }
      disparity[i][j]=temp;
      // std::cout << temp << std::endl;
    }
  }
  return disparity;
}



int main()
{
  int size=7;
  double sigma = 1;
  vpImage<unsigned char> Ileft(288,384);
  vpImage<unsigned char> Iright(288,384);

  vpImage<float> disparity(288,384,0);
  vpImage<unsigned char> disparity_UC(288,384,0);

  vpImage<double> kernel(size,size,0);
 
  vpImageIo::read(Ileft,"../Data/scene_l.pgm") ;
  vpImageIo::read(Iright,"../Data/scene_r.pgm") ;
  
  gaussianKernel(size, kernel);

//=======================================================================================
// Winner takes all
//=======================================================================================

  disparity = winnerTakesAll(Ileft, Iright);
  
  vpImageConvert::convert(disparity,disparity_UC);

  vpDisplayX d1(disparity_UC);
  vpDisplay::display(disparity_UC) ;
  vpDisplay::flush(disparity_UC) ;
  vpDisplay::getClick(disparity_UC) ;

  vpImageIo::write(disparity_UC, "../result/disparity_WTA.pgm");

//=======================================================================================
// SSD
//=======================================================================================

  disparity = SSD(Ileft, Iright, size, kernel);

  vpImageConvert::convert(disparity,disparity_UC);

  vpDisplayX d2(disparity_UC);
  vpDisplay::display(disparity_UC) ;
  vpDisplay::flush(disparity_UC) ;
  vpDisplay::getClick(disparity_UC) ;

  vpImageIo::write(disparity_UC, "../result/disparity_SSD.pgm");

  return 0;
}
