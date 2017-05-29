#include <iostream>

#include <visp/vpDebug.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageSimulator.h>
#include <visp/vpDisplayX.h>
#include <visp/vpDot.h>


using namespace std ;

// Calculer l'homographie aHb a partir des coordonnees des point p1 et p2
void DLT(unsigned int n,
	 vpImagePoint *p1,
	 vpImagePoint *p2,
	 vpMatrix &H12)
{

	#define NBPTMIN 5 
	if(n < NBPTMIN )
	{
		cout << "there must be at least " << NBPTMIN <<  " points in the both images\n" <<endl  ;
	    throw ;
	}  

	vpMatrix A(2*n,9);
	vpMatrix temp(2*n,9);
	vpMatrix V(9,9);
	vpColVector D(9);
	vpColVector h(9);
	H12.resize(3,3);

	float x1, y1, w1, x2, y2, w2;

	for (int i = 0; i < n; i++)
	{
		x2 = p2[i].get_u();
		y2 = p2[i].get_v();
		w2 = 1.0;
		x1 = p1[i].get_u();
		y1 = p1[i].get_v();
		w1 = 1.0;


		A[2*i][0]= 0. ;
        A[2*i][1]= 0. ;
        A[2*i][2]= 0. ;
        A[2*i][3]= - w1 * x2;
        A[2*i][4]= - w1 * y2;
        A[2*i][5]= - w1 * w2;
        A[2*i][6]= y1 * x2;
        A[2*i][7]= y1 * y2;
        A[2*i][8]= y1 * w2;
   
        A[2*i+1][0]= w1 * x2;
        A[2*i+1][1]= w1 * y2;
        A[2*i+1][2]= w1 * w2;
        A[2*i+1][3]= 0. ;
        A[2*i+1][4]= 0. ;
        A[2*i+1][5]= 0. ;
        A[2*i+1][6]= - x1 * x2;
        A[2*i+1][7]= - x1 * y2;
        A[2*i+1][8]= - x1 * w2;
	}

	for(unsigned int i =0;i<A.getRows();i++)
    {
        for(unsigned int j=0;j<A.getCols();j++)	temp[i][j]=A[i][j];
    }
	std::cout << A << std::endl;
	A.svd(D,V);

	double min = 9999999;
	int index = 0;
	for (int i = 0; i < D.getRows(); i++)
	{
		if(D[i] < min)
		{
			min = D[i];
			index = i;
		}
		// std:: cout << "index : " << i << std::endl;
		// std:: cout << "min : " << min << std::endl;
	}
	
	h=-V.getCol(index);

	// std::cout << "D" << std::endl;
	// std::cout << temp*h << std::endl;
	// std::cout  << "V" << std::endl;
	// std::cout << V << std::endl;
	std::cout  << "h" << std::endl;
	std::cout << h << std::endl;


    for(unsigned int i =0;i<3;i++)
    {
        for(unsigned int j=0;j<3;j++)	H12[i][j]=h[3*i+j];
    }

}

void transfer(const vpImage<unsigned char>& I1, const vpImage<unsigned char>& I2, const vpMatrix H12, vpImage<unsigned char>& result)
{

	int h = I1.getRows();
	int w = I1.getCols();
	int ii, jj;
	std::cout << "h: " << h << " w: " << w << std::endl; 

	// result.resize(2*h, 2*w, 0);
	

	// for (int i = 0; i < h; ++i)
	// {
	// 	for (int j = 0; j < w; ++j)
	// 	{
	// 	    vpColVector P2(3);
	// 	    vpColVector P1(3);

	// 	    P2[0] = i;
	// 	    P2[1] = j;
	// 	    P2[2] = 1;

	// 		P1 = H12*P2;

	// 		ii=P1[0]/P1[2];
	// 		jj=P1[1]/P1[2];
			
	// 		std::cout << "i: " << i << " j: " <<j << " ii: " << ii << " jj: " << jj ;
	// 		if( !(ii < -h/2 || ii > 3*h/2 || jj < -w/2 || jj > 3*w/2))
	// 		{ 
	// 			if(ii >= 0 && ii < h && ii >= 0 && jj < w)
	// 			{
	// 				std::cout << " Dans le if" << std::endl;
	// 				result[ii+h/2][jj+w/2] = I1[ii][jj]/*(I1[ii][jj] + I2[i][j])/2*/;
	// 			}
						
	// 			else	
	// 			{	
	// 				std::cout << " Dans le else" << std::endl;
	// 				//result[ii+h/2][jj+w/2] = I2[i][j];
	// 			}
	// 		}
	// 		else std::cout << std::endl;
	// 	}
	// }
	result.resize(h, w, 0);
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			result[i][j]=I1[i][j];
		}
	}

	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			vpColVector P2(3);
		    vpColVector P1(3);

		    P2[0] = j;
		    P2[1] = i;
		    P2[2] = 1;

			P1 = H12*P2;

			jj=P1[0]/P1[2];
			ii=P1[1]/P1[2];

			if( ii >= 0 && ii < h && jj >= 0 && jj < w)
				result[ii][jj] = (I2[i][j]+I1[ii][jj])/2;
		}
	}
}


int main()
{
  vpImage<unsigned char> I1;
  vpImage<unsigned char> I2;
  vpImage<unsigned char> Iresult;
  vpImage<vpRGBa> Iimage(876,1200);
  
 
  vpImageIo::read(I1,"../Data/I1.pgm") ;
  vpImageIo::read(I2,"../Data/I2.pgm") ;



  vpCameraParameters cam(800.0, 800.0, 200, 150);
  cam.printParameters() ;

  vpDisplayX d1(I1,10,10,"I1") ;
  vpDisplay::display(I1) ;
  vpDisplay::flush(I1) ;

  vpDisplayX d2(I2,450,10,"I2") ;
  vpDisplay::display(I2) ;
  vpDisplay::flush(I2) ;

  int nb = 5;
  vpImagePoint p1[nb], p2[nb];
  
  // clicker 5 point sur l'image I2 ; recuperer leur coordonnees
  for(unsigned int i=0; i<nb; i++)
    {
      vpDisplay::getClick(I1, p1[i]) ;
      vpDot d ;
      d.initTracking(I1,p1[i]) ;
      d.track(I1,p1[i]) ;
      char s[10] ;
      sprintf(s,"%d",i) ;
      vpDisplay::displayCross(I1,p1[i],10,vpColor::blue) ;
      vpDisplay::displayCharString(I1,p1[i],s,vpColor::red) ;
      vpDisplay::flush(I1) ;
    }
  
  // clicker 5 point sur l'image I1 ; recuperer leur coordonnees
  // faites attention a les clicker dans le meme ordre
  for(unsigned int i=0; i<nb; i++)
    {
      vpDisplay::getClick(I2, p2[i]) ;
      vpDot d ;
      d.initTracking(I2,p2[i]) ;
      d.track(I2,p2[i]) ;
      char s[10] ;
      sprintf(s,"%d",i) ;
      vpDisplay::displayCross(I2,p2[i],10,vpColor::green) ;
      vpDisplay::displayCharString(I2,p2[i],s,vpColor::red) ;
      vpDisplay::flush(I2) ;
    }
  

  // Calculer l'homographie
  vpMatrix H12;
  DLT(nb,p1, p2, H12) ;

  cout << "Homographie H12 : " <<endl ; 
  cout << H12 << endl ;
  cout << endl;

  //Verification 
  double residue =0 ;
  for (int i=0 ; i < nb ; i++) 
  {
  // Connaissant le formule permettant le transfert des points p2 dans p1
  // Calculer les coordonnées des point p1 connaissant p2 et dHg
    vpImagePoint p1_calcule  ;

    vpColVector P2(3);
    vpColVector P1(3);

    P2[0] = p2[i].get_u();
    P2[1] = p2[i].get_v();
    P2[2] = 1;

	P1 = H12*P2;

	P1[0]/=P1[2];
	P1[1]/=P1[2];
	P1[2]/=P1[2];

    p1_calcule.set_u(P1[0]);
    p1_calcule.set_v(P1[1]);

	// en deduire l'erreur sur commise sur chaque point et 
    // afficher un cercle de rayon 10 fois cette erreur
    double r ;
    r = vpImagePoint::distance(p1[i],p1_calcule) ;

    cout << "point " <<i << "  " << r <<endl ;;
    double rayon ;
    rayon = sqrt(r)*10 ; if (rayon < 10) rayon =10 ;
    vpDisplay::displayCircle(I1,p1_calcule,rayon,vpColor::green) ; ;
   }

   transfer(I1,I2,H12,Iresult);

  vpDisplayX dresult(Iresult,450,450,"Resultat fusion") ;
  vpDisplay::display(Iresult) ;
  vpDisplay::flush(Iresult) ;
  vpImageIo::write(Iresult,"resultat_fusion.jpg") ;

  vpDisplay::flush(I1) ;
  vpImage<vpRGBa> Ic ;
  vpDisplay::getImage(I1,Ic) ;
  vpImageIo::write(Ic,"resultat_distance.jpg") ;

  vpDisplay::getClick(I1) ;

  vpDisplay::close(I2) ;
  vpDisplay::close(I1) ;
  vpDisplay::close(Iresult) ;


  

  return 0;
}
