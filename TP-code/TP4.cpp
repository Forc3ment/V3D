#include <iostream>

#include <visp/vpDebug.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageSimulator.h>
#include <visp/vpDisplayX.h>
#include <visp/vpKeyPointSurf.h>
#include <visp/vpHomography.h>
#include <algorithm>



using namespace std ;

void AfficheAppariement(vpImage<unsigned char> &I1,  
			vpImage<unsigned char> &I2, 
			vpImage<unsigned char> &I,
      vpImagePoint* p1,
      vpImagePoint* p2,
      int nb
      )

{

  for (int i =0 ; i < I1.getRows() ; i++)
    for (int j = 0 ; j < I1.getCols() ; j++)
    {
      	I[i][j] = I1[i][j] ;
      	I[i][j+I1.getCols()] = I2[i][j] ;
    }

  vpDisplay::display(I) ;

  vpDisplay::flush(I) ;

  
  for(unsigned int i=0; i<nb; i++)
  {
    char s[10] ;
    sprintf(s,"%d",i) ;
    p2[i].set_u(p2[i].get_u()+I1.getCols()) ;
    vpDisplay::displayCharString(I,p1[i],s,vpColor::red) ;
    vpDisplay::displayCharString(I,p2[i],s,vpColor::red) ;
    vpDisplay::displayLine(I,p1[i],p2[i],vpColor::yellow) ;
  } 

  vpDisplay::flush(I) ;
  vpImage<vpRGBa> Ic ;
  vpDisplay::getImage(I,Ic) ;
  vpImageIo::write(Ic,"TP4_appariement.jpg") ;
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
  //  for (int j = 0; j < w; ++j)
  //  {
  //      vpColVector P2(3);
  //      vpColVector P1(3);

  //      P2[0] = i;
  //      P2[1] = j;
  //      P2[2] = 1;

  //    P1 = H12*P2;

  //    ii=P1[0]/P1[2];
  //    jj=P1[1]/P1[2];
      
  //    std::cout << "i: " << i << " j: " <<j << " ii: " << ii << " jj: " << jj ;
  //    if( !(ii < -h/2 || ii > 3*h/2 || jj < -w/2 || jj > 3*w/2))
  //    { 
  //      if(ii >= 0 && ii < h && ii >= 0 && jj < w)
  //      {
  //        std::cout << " Dans le if" << std::endl;
  //        result[ii+h/2][jj+w/2] = I1[ii][jj]/*(I1[ii][jj] + I2[i][j])/2*/;
  //      }
            
  //      else  
  //      { 
  //        std::cout << " Dans le else" << std::endl;
  //        //result[ii+h/2][jj+w/2] = I2[i][j];
  //      }
  //    }
  //    else std::cout << std::endl;
  //  }
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

void DLT(unsigned int n,
   const vpImagePoint *p1,
   const vpImagePoint *p2,
   vpMatrix &H12)
{
  std::cout << "--------- DLT ---------"<<endl;
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
  H12.resize(3,3,0);

  float x1, y1, w1, x2, y2, w2;

  for (int i = 0; i < n; i++)
  {
    x1 = p2[i].get_u();
    y1 = p2[i].get_v();
    w1 = 1.0;
    x2 = p1[i].get_u();
    y2 = p1[i].get_v();
    w2 = 1.0;


    A[2*i][0]= 0. ;
        A[2*i][1]= 0. ;
        A[2*i][2]= 0. ;
        A[2*i][3]= - w2 * x1;
        A[2*i][4]= - w2 * y1;
        A[2*i][5]= - w2 * w1;
        A[2*i][6]= y2 * x1;
        A[2*i][7]= y2 * y1;
        A[2*i][8]= y2 * w1;
   
        A[2*i+1][0]= w2 * x1;
        A[2*i+1][1]= w2 * y1;
        A[2*i+1][2]= w2 * w1;
        A[2*i+1][3]= 0. ;
        A[2*i+1][4]= 0. ;
        A[2*i+1][5]= 0. ;
        A[2*i+1][6]= - x2 * x1;
        A[2*i+1][7]= - x2 * y1;
        A[2*i+1][8]= - x2 * w1;
  }

  for(unsigned int i =0;i<A.getRows();i++)
    {
        for(unsigned int j=0;j<A.getCols();j++) temp[i][j]=A[i][j];
    }

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
  }
  
  h=-V.getCol(index);

  for(unsigned int i =0;i<3;i++)
  {
    for(unsigned int j=0;j<3;j++) H12[i][j]=h[3*i+j];
  }
  std::cout << "H12" << std::endl;
  std::cout << H12 << "\n" <<std::endl;
  std::cout << "--------- End DLT ---------"<<endl;

}

void getCorrespondances(vpImagePoint * p1, vpImagePoint * p2)
{

    p1[0].set_u(117.5130997);    p1[0].set_v(62.34123611);    p2[0].set_u(202.841095);     p2[0].set_v(36.29648209);

    p1[1].set_u(84.06044006);    p1[1].set_v(67.55551147);    p2[1].set_u(169.5350189);    p2[1].set_v(26.80556679);

    p1[2].set_u(80.27194214);    p1[2].set_v(111.0672302);    p2[2].set_u(147.9641113);    p2[2].set_v(64.5475769);

    p1[3].set_u(342.6855164);    p1[3].set_v(199.8661346);    p2[3].set_u(63.4621048);     p2[3].set_v(68.28819275);

    p1[4].set_u(302.6676636);    p1[4].set_v(226.6687317);    p2[4].set_u(300.4017639);    p2[4].set_v(263.6835022);

    p1[5].set_u(101.5870972);    p1[5].set_v(63.0242424);     p2[5].set_u(187.8421478);    p2[5].set_v(29.56011963);

    p1[6].set_u(153.4119415);    p1[6].set_v(91.05652618);    p2[6].set_u(222.968277);     p2[6].set_v(77.2434845);

    p1[7].set_u(190.6780548);    p1[7].set_v(110.7231598);    p2[7].set_u(247.8312683);    p2[7].set_v(110.4263763);

    p1[8].set_u(302.8087463);    p1[8].set_v(133.9337616);    p2[8].set_u(339.9194641);    p2[8].set_v(178.880661);

    p1[9].set_u(162.7279968);    p1[9].set_v(276.4970398);    p2[9].set_u(152.7050171);    p2[9].set_v(248.9367065);

    p1[10].set_u(151.0850067);   p1[10].set_v(36.12360764);   p2[10].set_u(244.672287);    p2[10].set_v(25.44586563);

    p1[11].set_u(171.7740173);   p1[11].set_v(53.67162704);   p2[11].set_u(256.0083618);   p2[11].set_v(49.99362183);

    p1[12].set_u(116.7895355);   p1[12].set_v(74.19098663);   p2[12].set_u(196.8202972);   p2[12].set_v(45.97808456);

    p1[13].set_u(104.2023163);   p1[13].set_v(83.85998535);   p2[13].set_u(181.4200439);   p2[13].set_v(50.26084518);

    p1[14].set_u(84.71365356);   p1[14].set_v(190.8507233);   p2[14].set_u(300.4017639);   p2[14].set_v(263.6835022);

    p1[15].set_u(138.8526764);   p1[15].set_v(273.5761719);   p2[15].set_u(131.6974182);   p2[15].set_v(236.8515778);

    p1[16].set_u(167.2081451);   p1[16].set_v(96.59983063);   p2[16].set_u(233.1238556);   p2[16].set_v(88.96112061);

}

void generateRandomSubset(int setSize, int subSetSize, std::vector<int>& subset)
{
  for (int i = 0; i < subSetSize; i++)
  {
    int rnd = rand() % setSize;
    if(std::find(subset.begin(), subset.end(), rnd) == subset.end()) 
    {
      subset.push_back(rnd);
    } 
    else 
    {
      i--;
    }
  }
}

double computeError(const std::vector<int> subset, const vpImagePoint* p1, const vpImagePoint* p2, const vpMatrix H12)
{
  double error = 0;
  std::cout << "--------- computeError ---------"<<endl;
  // std::cout << "H12" << std::endl;
  // std::cout << H12 << "\n" <<std::endl;

  for (int i=0 ; i < subset.size() ; i++) 
  {
  // Connaissant le formule permettant le transfert des points p2 dans p1
  // Calculer les coordonnées des point p1 connaissant p2 et dHg
    vpImagePoint p1_calcule  ;

    vpColVector P2(3);
    vpColVector P1(3);

    P2[0] = p2[subset[i]].get_u();
    P2[1] = p2[subset[i]].get_v();
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
    r = vpImagePoint::distance(p1[subset[i]],p1_calcule);
    cout << "point " <<i << "  " << sqrt(r) <<endl;
    error += sqrt(r);
  }
  cout << "Subset error : " << error/subset.size() << endl;
  std::cout << "--------- End computeError ---------"<<endl;
  return error/subset.size();
}

int checkError(const vpImagePoint* p1, const vpImagePoint* p2, const vpMatrix H12, const double seuil, std::vector<bool>& ok)
{
  int somme = 0;
  std::cout << "--------- checkError ---------"<<endl;
  // std::cout << "H12" << std::endl;
  // std::cout << H12 << "\n" <<std::endl;

  for (int i=0 ; i < ok.size() ; i++) 
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
    cout << "point " <<i << "  " << sqrt(r) <<endl;
    if(sqrt(r) < seuil) 
    {
      ok[i] = true;
      somme++ ;
    }
  }
  std::cout << "--------- End checkError ---------"<<endl;
  return somme;

}

std::vector<bool> ransac(const vpImagePoint* p1, const vpImagePoint* p2, vpMatrix & H12, int nbPoint=17, int minSize=5, int seuilSubsetError=50, int seuilIndiv=10, int nbIter=40)
{
  std::vector<bool> ok(nbPoint,false);
  
  H12.resize(3,3,0);

  double sub_error = 0;
  int temp = 0;
  int nbMax = 0;

  for (int i = 0; i < nbIter; ++i)
  {
    std::cout << "------------------- Iteration "<<i<<"-------------------" << std::endl;

    std::vector<int> subset;
    std::vector<bool> temp_ok(nbPoint,false);
    vpImagePoint temp1[minSize], temp2[minSize];

    generateRandomSubset(nbPoint, minSize, subset);
    
    int k = 0;
    for (int i = 0; i < minSize; ++i)
    {
      temp1[i]=p1[subset[i]];
      temp2[i]=p2[subset[i]];
    }



    DLT(minSize,temp1,temp2,H12);

    //std::cout << "H12" << std::endl;
    //std::cout << H12 << "\n" <<std::endl;

    sub_error = computeError(subset, p1, p2, H12);

    std::cout << "Erreur subset :"<< sub_error << std::endl;
    if(sub_error < seuilSubsetError)
    {
      temp = checkError(p1, p2, H12, seuilIndiv, temp_ok);
      std::cout << "Nb point passant inf au seuil :"<< temp << std::endl;
      if(temp>nbMax)
      {
        nbMax=temp;
        for (int i = 0; i < ok.size(); ++i)
        {
          ok[i]=temp_ok[i];

        } 
      }
    }
  }

  for (int i = 0; i < ok.size(); ++i)
    if(ok[i]) std::cout << i << std::endl;

  vpImagePoint res1[nbMax], res2[nbMax];
  int k = 0;
  for (int i = 0; i < ok.size(); ++i)
  {
    if(ok[i])
    {
      res1[k]=p1[i];
      res2[k]=p2[i];
      k++;
    }
  }

  DLT(nbMax,res1,res2,H12);
  return ok;

}

int main()
{
  vpImage<unsigned char> I1(300,400,0);
  vpImage<unsigned char> I2(300,400,0);
  vpImage<unsigned char> Iresult;
  vpImage<vpRGBa> Iimage(876,1200);

  std::vector<bool> ok;
  
  vpMatrix H12;
 
  vpImageIo::read(Iimage,"../Data/big-sleep.jpg") ;
  

  double L = 0.600 ;
  double l = 0.438;
  // Initialise the 3D coordinates of the Iimage corners
  vpColVector X[4];
  for (int i = 0; i < 4; i++) X[i].resize(3);
  // Top left corner
  X[0][0] = -L;
  X[0][1] = -l;
  X[0][2] = 0;
  
  // Top right corner
  X[1][0] = L;
  X[1][1] = -l;
  X[1][2] = 0;
  
  // Bottom right corner
  X[2][0] = L;
  X[2][1] = l;
  X[2][2] = 0;
  
  //Bottom left corner
  X[3][0] = -L;
  X[3][1] = l;
  X[3][2] = 0;
  


  vpImageSimulator sim;
  sim.init(Iimage, X);
  vpCameraParameters cam(800.0, 800.0, 200, 150);
  

  cam.printParameters() ;


  // I1g
  vpHomogeneousMatrix  c1gMo(0,0,2,  vpMath::rad(0),vpMath::rad(0),0) ;
  sim.setCameraPosition(c1gMo);
  sim.getImage(I1,cam);
  cout << "Image I1g " <<endl ;
  cout << c1gMo << endl ;

  // I1d
  vpHomogeneousMatrix c1dMo(0.1,0,2, 
			    vpMath::rad(0),vpMath::rad(0),vpMath::rad(25)) ; //0.1,0,2, vpMath::rad(0),vpMath::rad(0),0) ;
  sim.setCameraPosition(c1dMo);
  sim.getImage(I2,cam);  
  cout << "Image I1d " <<endl ;
  cout << c1dMo << endl ;
 
  vpHomogeneousMatrix cgMcd = c1gMo * c1dMo.inverse() ;
  vpMatrix K = cam.get_K() ;

  vpDisplayX dg(I1,10,10,"I1") ;
  vpDisplay::display(I1) ;
  vpDisplay::flush(I1) ;

  vpDisplayX dd(I2,450,10,"I2") ;
  vpDisplay::display(I2) ;
  vpDisplay::flush(I2) ;

  // Image resultat
  vpImage<unsigned char> I ;
  I.resize(I1.getRows(), I1.getCols()*2) ;
  vpDisplayX d(I,10,400,"I") ;
  
  int nb = 17;
  vpImagePoint p1[nb], p2[nb];
  getCorrespondances(p1, p2);

  
  AfficheAppariement(I1,  I2, I, p1, p2, nb) ;

  if(nb >= 0){ // ... add paired points to vectPts
    
    for(unsigned int i=0; i<nb; i++)
    {
    	char s[10] ;
    	sprintf(s,"%d",i) ;
    	// cout << i <<"  "  << p1[i].get_u() <<"  " << p1[i].get_v() <<"  " ;
    	// cout <<  p2[i].get_u() <<"  " << p2[i].get_v() << endl;
    	vpDisplay::displayCharString(I1,p1[i],s,vpColor::yellow) ;
    	vpDisplay::displayCharString(I2,p2[i],s,vpColor::yellow) ;
    }

  }

  vpDisplay::getClick(I);

  ok = ransac(p1, p2, H12);

  std::cout << "\nH12" << std::endl;
  std::cout << H12 << std::endl;

  for (int i=0 ; i < nb ; i++) 
  {
    if(ok[i])
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

      cout << "point " <<i << "  " << sqrt(r) <<endl;
      double rayon ;
      rayon = sqrt(r)*10 ; if (rayon < 10) rayon =10 ;
      vpDisplay::displayCircle(I1,p1_calcule,rayon,vpColor::green);
    }    
  }

  transfer(I1,I2,H12,Iresult);

  vpDisplayX dresult(Iresult,450,450,"Resultat fusion") ;
  vpDisplay::display(Iresult) ;
  vpDisplay::flush(Iresult) ;
  vpImageIo::write(Iresult,"TP4_resultat_fusion.jpg") ;

  vpDisplay::flush(I1) ;
  vpImage<vpRGBa> Ic ;
  vpDisplay::getImage(I1,Ic) ;
  vpImageIo::write(Ic,"TP4_resultat_distance.jpg") ;

  vpDisplay::getClick(I1) ;

  vpDisplay::close(I2) ;
  vpDisplay::close(I1) ;
  vpDisplay::close(Iresult) ;
  
  vpDisplay::getClick(I);

  return 0;
}
